# ─────────────────────────────────────────────────────────────────────────────
# 0) Load required packages
# ─────────────────────────────────────────────────────────────────────────────
library(vroom)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(broom)      # for tidy(t.test)
library(forcats)
library(patchwork)  # for laying out the final figure

# ─────────────────────────────────────────────────────────────────────────────
# 1) Read in & wrangle the two Z-scored datasets
# ─────────────────────────────────────────────────────────────────────────────
control  <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Control_all_lipids_final_non_normalized.csv")
lowinput <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv")

valid_classes <- c(
  "TG", "DG", "MG",
  "PC", "PE", "PI",
  "DGDG", "MGDG",
  "SQDG","SM","AEG","LPC","LPE","PG","PA"
)
class_pattern <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

process_long <- function(df, condition_label) {
  df %>%
    pivot_longer(
      cols      = -Compound_Name,
      names_to  = "Lipid",
      values_to = "Intensity"
    ) %>%
    mutate(
      Class     = str_extract(Lipid, class_pattern),
      Condition = condition_label
    ) %>%
    filter(!is.na(Class), Class %in% valid_classes) %>%
    # **WE DO NOT back-transform**; we only have Z-scores, so we sum those.
    group_by(Compound_Name, Condition, Class) %>%
    summarise(
      sum_intensity = sum(Intensity, na.rm = TRUE),
      .groups       = "drop"
    ) %>%
    rename(Sample = Compound_Name)
}


control_long  <- process_long(control,  "Control")

lowinput_long <- process_long(lowinput, "LowInput")

per_sample_sums <- bind_rows(control_long, lowinput_long)

per_sample_wide <- per_sample_sums %>%
  select(Sample, Condition, Class, sum_intensity) %>%
  pivot_wider(
    names_from  = Class,
    values_from = sum_intensity,
    values_fill = list(sum_intensity = 0)
  )
# Now per_sample_wide columns are:
#   Sample | Condition | DG | DGDG | MG | MGDG | PC | PE | PI | SQDG | TG
# Where each column (DG, DGDG, etc.) is “sum of Z-scores within that class.”

# ─────────────────────────────────────────────────────────────────────────────
# 2) Compute the 10 “important” ratios per sample (no LPC/LPE here)
# ─────────────────────────────────────────────────────────────────────────────
per_sample_ratios <- per_sample_wide %>%
  mutate(
    ## ─ Low-P (6 ratios)
    DGDG_PC    = ifelse(PC   > 0, DGDG      / PC,     NA_real_),
    DGDG_PE    = ifelse(PE   > 0, DGDG      / PE,     NA_real_),
    MGDG_PC    = ifelse(PC   > 0, MGDG      / PC,     NA_real_),
    MGDG_PE    = ifelse(PE   > 0, MGDG      / PE,     NA_real_),
    SQDG_PC    = ifelse(PC   > 0, SQDG      / PC,     NA_real_),
    SQDG_PE    = ifelse(PE   > 0, SQDG      / PE,     NA_real_),
    PG_PC          = ifelse(PC  > 0, PG  / PC,  NA_real_),
    LPC_PC         = ifelse(PC  > 0, LPC / PC, NA_real_),
    LPE_PE         = ifelse(PE  > 0, LPE / PE, NA_real_),
    nonP_P_total   = (MGDG + DGDG + SQDG) /
      (PC + PE + PG + PA + LPC + LPE),
    
    
    ## ─ Low-N (3 ratios)
    TG_PC      = ifelse(PC   > 0, TG        / PC,     NA_real_),
    DG_PC      = ifelse(PC   > 0, DG        / PC,     NA_real_),
    TG_DG      = ifelse(DG   > 0, TG        / DG,     NA_real_),
    stor_vs_photo  = (TG + DG) / (MGDG + PC),
    DG_Phospho     = ifelse((PC + PE) > 0, DG / (PC + PE), NA_real_),
    
    ## ─ Cold acclimation (3 ratios)
    MGDG_DGDG  = ifelse(DGDG > 0, MGDG      / DGDG,   NA_real_),
    SQDG_DGDG  = ifelse(DGDG > 0, SQDG      / DGDG,   NA_real_),
    SQDG_MGDG  = ifelse(MGDG > 0, SQDG      / MGDG,   NA_real_),
    DGDG_total     = ifelse((MGDG + SQDG) > 0,
                            DGDG / (MGDG + SQDG), NA_real_),
    
    ## ─ Membrane integrity (1 ratio)
    PC_PE      = ifelse(PE   > 0, PC        / PE,     NA_real_),
    Lyso_ratio     = (LPC + LPE) / (PC + PE + PG),
  )

# ─────────────────────────────────────────────────────────────────────────────
# 3) Pivot those 10 ratio columns into “long” form
# ─────────────────────────────────────────────────────────────────────────────
ratio_long <- per_sample_ratios %>%
  select(
    Sample, Condition,
    DGDG_PC, DGDG_PE, MGDG_PC, MGDG_PE, SQDG_PC, SQDG_PE, PG_PC, LPC_PC, LPE_PE, nonP_P_total,
    TG_PC,   DG_PC,   TG_DG,stor_vs_photo, DG_Phospho,
    MGDG_DGDG, SQDG_DGDG, SQDG_MGDG,DGDG_total,
    PC_PE,Lyso_ratio
  ) %>%
  pivot_longer(
    cols      = -c(Sample, Condition),
    names_to  = "Ratio",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# ─────────────────────────────────────────────────────────────────────────────
# 4) Run t-tests per Ratio → extract p-values & stars, compute y_pos
# ─────────────────────────────────────────────────────────────────────────────
stat_df <- ratio_long %>%
  group_by(Ratio) %>%
  summarise(
    p_value = t.test(Value ~ Condition)$p.value,
    y_max   = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    star = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    y_pos = y_max * 1.10   # place stars 10% above max
  ) %>%
  filter(star != "ns")

# ─────────────────────────────────────────────────────────────────────────────
# 5) Define your “nature_theme” (remove all panel borders/background)
# ─────────────────────────────────────────────────────────────────────────────
nature_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(
      size   = 14,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title       = element_text(size = 12, face = "bold"),
    axis.text        = element_text(color = "black"),
    axis.text.y      = element_text(face = "bold", hjust = 1),
    axis.line        = element_line(color = "black"),
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    panel.background = element_blank(),
    legend.position  = "none",
    plot.margin      = margin(12, 12, 12, 12)
  )

# 6) Revised helper to build a boxplot WITHOUT dots, on a linear y‐axis

condition_colors <- c(Control = "#440154FF", LowInput = "#FDE725FF")

make_category_plot <- function(ratio_names, ncol_facet) {
  df_sub <- filter(ratio_long, Ratio %in% ratio_names)
  
  ggplot(df_sub, aes(x = Condition, y = Value, fill = Condition, color = Condition)) +
    # show the full distribution
    geom_violin(
      position = position_dodge(width = 0.8),
      width    = 0.9,
      trim     = TRUE,
      alpha    = 0.3,
      size     = 0
    ) +
    # then overlay boxplots (no outliers)
    geom_boxplot(
      outlier.shape = NA,
      width         = 0.6,
      colour        = "black",
      alpha         = 0.8,
      position      = position_dodge(width = 0.8)
    ) +
    facet_wrap(~ Ratio, ncol = ncol_facet, scales = "free_y", drop = FALSE) +
    scale_fill_manual(values = condition_colors) +
    scale_color_manual(values = condition_colors) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
    labs(x = NULL, y = "Lipid-class ratio") +
    nature_theme +
    theme(
      strip.text            = element_text(face = "bold", size = 10),
      panel.spacing         = unit(0.5, "lines"),
      # add light grey grid lines on both axes
      panel.grid.major.x     = element_line(color = "grey90"),
      panel.grid.minor.x     = element_line(color = "grey95"),
      panel.grid.major.y     = element_line(color = "grey90"),
      panel.grid.minor.y     = element_line(color = "grey95")
    ) +
    # significance stars
    geom_text(
      data = filter(stat_df, Ratio %in% ratio_names),
      aes(x = 1.5, y = y_pos, label = star),
      inherit.aes = FALSE,
      size        = 5,
      colour      = "black"
    )
}


# example: redraw the Low-P panel
p_lowP <- make_category_plot(ratios_lowP, ncol_facet = 2)
quartz()
print(p_lowP)

p_lowN <- make_category_plot(ratios_lowN, ncol_facet = 2)
quartz()
print(p_lowN)

p_cold <- make_category_plot(ratios_cold, ncol_facet = 2)
quartz()
print(p_cold)

p_membrane <- make_category_plot(ratios_membrane, ncol_facet = 2)
quartz()
print(p_membrane)

# ─────────────────────────────────────────────────────────────────────────────
# 7) Build one ggplot for each of the four categories (exactly as before)
# ─────────────────────────────────────────────────────────────────────────────
ratios_lowP     <- c("DGDG_PC","DGDG_PE","MGDG_PC","MGDG_PE","SQDG_PC","SQDG_PE",
                     "PG_PC","LPC_PC","LPE_PE","nonP_P_total")
p_lowP          <- make_category_plot(ratios_lowP,     ncol_facet = 3)

ratios_membrane <- c("PC_PE", "Lyso_ratio")
p_membrane      <- make_category_plot(ratios_membrane, ncol_facet = 1)

ratios_lowN     <- c("TG_PC","DG_PC","TG_DG", "stor_vs_photo", "DG_Phospho")
p_lowN          <- make_category_plot(ratios_lowN,     ncol_facet = 3)

ratios_cold     <- c("MGDG_DGDG","SQDG_DGDG","SQDG_MGDG","DGDG_total")
p_cold          <- make_category_plot(ratios_cold,     ncol_facet = 3)

# ─────────────────────────────────────────────────────────────────────────────
# 8) Build the “heading” grobs (unchanged)
# ─────────────────────────────────────────────────────────────────────────────
make_group_label <- function(label_text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = label_text,
             size = 6,
             fontface = "bold") +
    theme_void()
}
p_label_lowP     <- make_group_label("Low P")
p_label_membrane <- make_group_label("Membrane integrity")
p_label_lowN     <- make_group_label("Low N")
p_label_cold     <- make_group_label("Cold acclimation")

# ─────────────────────────────────────────────────────────────────────────────
# 9) Stack each heading above its plot with patchwork
# ─────────────────────────────────────────────────────────────────────────────
col1 <- p_label_lowP     / p_lowP     + plot_layout(heights = c(1, 4))
col2 <- p_label_membrane / p_membrane + plot_layout(heights = c(1, 4))
col3 <- p_label_lowN     / p_lowN     + plot_layout(heights = c(1, 4))
col4 <- p_label_cold     / p_cold     + plot_layout(heights = c(1, 4))

quartz()
col1
# ─────────────────────────────────────────────────────────────────────────────
# 10) Arrange all four columns in a 2×2 grid (no panel borders)
# ─────────────────────────────────────────────────────────────────────────────
final_figure <- (col1 | col2) /
  (col3 | col4) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  ) &
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_blank(),
    panel.border     = element_blank(),
    plot.margin      = margin(0, 0, 0, 0)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 11) Display and/or save
# ─────────────────────────────────────────────────────────────────────────────
quartz()
print(final_figure)

# If you want to write out to a PNG:
ggsave("Fig2_lipid_ratio_linear.png", final_figure,
       width = 12, height = 10, dpi = 300)
