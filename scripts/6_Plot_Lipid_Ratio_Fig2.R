# ╔══════════════════════════════════════════════════════════════════╗
# ║ 0)  PACKAGES                                                    ║
# ╚══════════════════════════════════════════════════════════════════╝
library(vroom)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(broom)      # tidy(t.test)
library(patchwork)  # arranging plots
library(scales)

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 1)  READ Z‑SCORED PEAK TABLES                                    ║
# ╚══════════════════════════════════════════════════════════════════╝
ctrl_file <- "/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Control_all_lipids_final_non_normalized.csv"
lowp_file <- "/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv"

control  <- vroom(ctrl_file,  show_col_types = FALSE)
lowinput <- vroom(lowp_file,  show_col_types = FALSE)

valid_classes <- c("TG","DG","MG",
                   "PC","PE","PI",
                   "DGDG","MGDG",
                   "SQDG","SM","AEG",
                   "LPC","LPE","PG","PA")
class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

long_prep <- function(df, label){
  df %>% 
    pivot_longer(-Compound_Name,
                 names_to  = "Lipid",
                 values_to = "Z") %>%          # intensities are already Z‑scores
    mutate(Class     = str_extract(Lipid, class_pat),
           Condition = label) %>% 
    filter(!is.na(Class)) %>% 
    group_by(Sample = Compound_Name, Condition, Class) %>%
    summarise(class_z = mean(Z, na.rm = TRUE), .groups = "drop") # mean Z per class
}

long_ctrl <- long_prep(control,  "Control")
long_lowp <- long_prep(lowinput, "LowInput")

wide_z <- bind_rows(long_ctrl, long_lowp) %>% 
  pivot_wider(names_from  = Class,
              values_from = class_z,
              values_fill = NA_real_)

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 2)  COMPUTE CLASS‑PAIR  Z‑DIFFERENCES (  zA – zB  )              ║
# ╚══════════════════════════════════════════════════════════════════╝
z_diff <- wide_z %>% 
  mutate(
    ## — Low‑P contrasts
    DGDG_PC        = DGDG - PC,
    DGDG_PE        = DGDG - PE,
    MGDG_PC        = MGDG - PC,
    MGDG_PE        = MGDG - PE,

    PG_retention   = PG - (PC + PE)/2,
    
    PC_PE          = PC    - PE,
  
    Lyso_activity  = (LPC + LPE) - (PC + PE),
    nonP_P_total   = (MGDG + DGDG) - (PC + PE),
    
    # SUPPLEMENTARY FIGURE
    #SQDG_Spared    = SQDG - (DGDG + MGDG + PG)/3,
    #LPC_PC         = LPC   - PC,
    #LPE_PE         = LPE   - PE,
    
    
    ## — Low‑N contrasts
    TG_PC          = TG - PC,
    DG_PC          = DG - PC,
    TG_DG          = TG - DG,
    DG_Phospho     = DG - (PC + PE)/2,
    TG_Phospho     = TG - (PC + PE)/2,
    stor_vs_photo  = (TG + DG)/2 - (MGDG + PC)/2,
    
    ## — Cold contrasts
    MGDG_DGDG      = MGDG - DGDG,
    SQDG_DGDG      = SQDG - DGDG,
    SQDG_MGDG      = SQDG - MGDG,
    DGDG_total     = DGDG - ((MGDG + SQDG)/2),
    MGDG_total     = MGDG - ((DGDG + SQDG)/2),
    SQDG_total     = SQDG - ((DGDG + MGDG)/2),
    
    ## — Membrane integrity (2)
    PC_PE          = PC - PE,
    PC_PG          = PC - PG,
    PE_PG          = PE - PG,
    PC_total       = (PC + PE + PG)/3,
    PE_total       = (PC + PE + PG)/3,
    PG_total       = (PC + PE + PG)/3,
    Lyso_ratio     = (LPC + LPE)/2 - (PC + PE + PG)/3
  )

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 3)  LONG FORMAT + t‑TESTS + p‑value STARS                        ║
# ╚══════════════════════════════════════════════════════════════════╝
ratio_long <- z_diff %>% 
  select(Sample, Condition,
         # Low P
         DGDG_PC,DGDG_PE,
         MGDG_PC,MGDG_PE,
         PG_retention,PC_PE,
         Lyso_activity,nonP_P_total,
         
         # Low P supplementary 
         #SQDG_Spared,LPC_PC,LPE_PE,
         
         TG_PC,DG_PC,TG_DG,
         DG_Phospho,TG_Phospho,stor_vs_photo,
         MGDG_DGDG,SQDG_DGDG,SQDG_MGDG,
         DGDG_total,MGDG_total,SQDG_total,
         PC_PE,PC_PG,PE_PG,Lyso_ratio,
         PC_total,PE_total,PG_total) %>% 
  pivot_longer(-c(Sample, Condition),
               names_to  = "Ratio",
               values_to = "Value") %>% 
  filter(!is.na(Value))

stat_df <- ratio_long %>% 
  group_by(Ratio) %>% 
  summarise(
    p_value = t.test(Value ~ Condition)$p.value,
    y_pos   = max(Value, na.rm = TRUE) + 0.25,
    .groups = "drop"
  ) %>% 
  mutate(
    star = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>% 
  filter(star != "ns")

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 4)  PLOTTING HELPER                                              ║
# ╚══════════════════════════════════════════════════════════════════╝
condition_cols <- c(Control = "#440154FF", LowInput = "#FDE725FF")

nature_theme <- theme_minimal(12) +
  theme(axis.line        = element_line(colour = "black"),
        axis.text        = element_text(colour = "black"),
        panel.grid.major = element_line(colour = "grey92"),
        panel.grid.minor = element_line(colour = "grey95"),
        legend.position  = "none")

make_plot <- function(ratio_vec, ncol_facets){
  df_sub <- filter(ratio_long, Ratio %in% ratio_vec) %>%
    mutate(Ratio = factor(Ratio, levels = ratio_vec))  # set custom facet order
  stat_sub <- filter(stat_df, Ratio %in% ratio_vec) %>%
    mutate(Ratio = factor(Ratio, levels = ratio_vec))  # ensure matching factor levels
  
  ggplot(df_sub, aes(x = Condition, y = Value, fill = Condition, colour = Condition)) +
    geom_violin(width = 0.9, alpha = 0.3, trim = TRUE, size = 0) +
    geom_boxplot(width = 0.6, alpha = 0.8, colour = "black", outlier.shape = NA) +
    facet_wrap(~ Ratio, ncol = ncol_facets, scales = "free_y") +
    geom_segment(
      data = stat_sub,
      inherit.aes = FALSE,
      aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
      colour = "black", size = 0.3
    ) +
    geom_text(
      data = stat_sub,
      inherit.aes = FALSE,
      aes(x = 1.5, y = y_pos * 1.02, label = star),
      size = 5, colour = "black"
    ) +
    scale_x_discrete(
      name = "Condition",
      labels = c("Control", "LowInput"),
      expand = expansion(add = 0.5)
    ) +
    scale_fill_manual(values = condition_cols) +
    scale_colour_manual(values = condition_cols) +
    labs(
      x = "Condition",
      y = "Δ Z-scores"
    ) +
    nature_theme +
    theme(
      plot.title         = element_text(face = "bold", hjust = 0.5),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold", size = 10),
      panel.spacing      = unit(0.5, "lines"),
      panel.grid.major.y = element_line(color = "grey90", size = 0.3),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.line.x.bottom = element_line(color = "black"),
      axis.line.x.top    = element_blank(),
      axis.line.y.left   = element_line(color = "black"),
      axis.line.y.right  = element_blank(),
      axis.ticks.x       = element_line(color = "black"),
      axis.ticks.y       = element_line(color = "black"),
      axis.text.x        = element_text(color = "black"),
      axis.text.y        = element_text(color = "black"),
      panel.border       = element_blank(),
      plot.margin        = margin(10, 10, 10, 10)
    )
}


# groupings
ratios_lowP <- c("DGDG_PC","DGDG_PE",
                 "MGDG_PC","MGDG_PE",
                 "PG_retention","PC_PE",
                 "Lyso_activity", "nonP_P_total")

# Low P supplementary
#ratios_lowP <- c("SQDG_Spared","LPC_PC","LPE_PE")

ratios_lowN      <- c("TG_PC","DG_PC","TG_DG",
                      "DG_Phospho","TG_Phospho","stor_vs_photo")
ratios_cold      <- c("MGDG_DGDG","SQDG_DGDG","SQDG_MGDG",
                      "DGDG_total","MGDG_total","SQDG_total")
ratios_membrane  <- c("PC_PE","PC_PG","PE_PG",
                      "PC_total","PE_total","PG_total",
                      "Lyso_ratio")

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 5)  DRAW & SAVE  FIGURES                                         ║
# ╚══════════════════════════════════════════════════════════════════╝
p_lowP <- make_plot(ratios_lowP,     2)
quartz()
p_lowP
ggsave("Fig2a_lipid_ratio_linear_lowP.png",  p_lowP,  width = 6, height = 12, dpi = 300, bg = "white")
#ggsave("SupFig_lipid_ratio_linear_lowP.png",  p_lowP,  width = 9, height = 3, dpi = 300, bg = "white")


p_lowN <- make_plot(ratios_lowN,     3)
quartz()
p_lowN
ggsave("Fig2b_lipid_ratio_linear_lowN.png",  p_lowN,  width = 12, height = 6,  dpi = 300, bg = "white")

p_cold <- make_plot(ratios_cold,     3)
quartz()
p_cold
ggsave("Fig2c_lipid_ratio_linear_cold.png",  p_cold, width = 12, height = 6,  dpi = 300, bg = "white")

p_memb <- make_plot(ratios_membrane, 3)
quartz()
p_memb
ggsave("Fig2d_lipid_ratio_linear_membrane.png", p_memb, width = 12, height = 9, dpi = 300, bg = "white")



library(dplyr)
library(tidyr)

# 1a) Compute median ΔZ per Ratio and Condition
median_df <- ratio_long %>%
  group_by(Ratio, Condition) %>%
  summarise(
    medianZ = median(Value, na.rm = TRUE),
    .groups = "drop"
  )

# 1b) Pivot so Control and LowInput are side by side, then compute the difference
median_diff <- median_df %>%
  pivot_wider(names_from = Condition, values_from = medianZ) %>%
  # Adjust the names “Control” and “LowInput” if your labels differ
  mutate(
    delta_median = LowInput - Control
  )

# 1c) Inspect only the Low-P ratios (replace ratios_lowP with your vector of names)
median_diff_lowP <- median_diff %>%
  filter(Ratio %in% ratios_lowP) %>%
  arrange(desc(delta_median))

print(median_diff_lowP, n = Inf)
