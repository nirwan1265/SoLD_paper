# ╔══════════════════════════════════════════════════════════════════╗
# ║ 0)  PACKAGES                                                    ║
# ╚══════════════════════════════════════════════════════════════════╝
library(vroom)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(broom)      
library(patchwork)  
library(scales)

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 1)  READ RAW INTENSITY TABLES                                   ║
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

# ╔══════════════════════════════════════════════════════════════════╗
# ║ Modified long_prep: TIC normalisation + log10 per lipid (NO Z)  ║
# ╚══════════════════════════════════════════════════════════════════╝
long_prep_global <- function(df) {
  df_long <- df %>%
    pivot_longer(-c(Compound_Name, Condition),
                 names_to = "Lipid", values_to = "Intensity") %>%
    rename(Sample = Compound_Name)
  
  # 1) TIC per sample → relative abundance
  df_long <- df_long %>%
    group_by(Sample) %>%
    mutate(TIC       = sum(Intensity, na.rm = TRUE),
           rel_abund = Intensity / TIC) %>%
    ungroup()
  
  # 2) small pseudo-count per sample → log10
  df_long <- df_long %>%
    group_by(Sample) %>%
    mutate(minpos = ifelse(
      all(is.na(rel_abund) | rel_abund <= 0),
      NA_real_, min(rel_abund[rel_abund > 0], na.rm = TRUE)),
      eps    = ifelse(is.na(minpos), 0, minpos * 0.5)) %>%
    ungroup() %>%
    mutate(log_rel = log10(rel_abund + eps))
  
  # 3) summarise mean log10(relative abundance) per lipid class
  df_long <- df_long %>%
    mutate(Class = str_extract(Lipid, class_pat)) %>%
    filter(!is.na(Class)) %>%
    group_by(Sample, Condition, Class) %>%
    summarise(class_log = mean(log_rel, na.rm = TRUE), .groups = "drop")
  
  return(df_long)
}

# combine Control + LowInput, then run prep
combined  <- bind_rows(control  %>% mutate(Condition = "Control"),
                       lowinput %>% mutate(Condition = "LowInput"))
long_all  <- long_prep_global(combined)

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 2)  PIVOT WIDE                                                  ║
# ╚══════════════════════════════════════════════════════════════════╝
wide_log <- long_all %>%
  pivot_wider(names_from = Class, values_from = class_log,
              values_fill = NA_real_)

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 3)  COMPUTE CLASS-PAIR  log10-RATIOS                            ║
# ╚══════════════════════════════════════════════════════════════════╝
ratio_tbl <- wide_log %>%
  mutate(
    ## — Low-P contrasts 
    DGDG_PC        = DGDG - PC,
    DGDG_PE        = DGDG - PE,
    MGDG_PC        = MGDG - PC,
    MGDG_PE        = MGDG - PE,

    PG_retention   = PG - (PC + PE)/2,

    PC_PE          = PC - PE,

    Lyso_activity  = (LPC + LPE) - (PC + PE),
    nonP_P_total   = (MGDG + DGDG) - (PC + PE),
    
    
    # Low P supplementary
    # SQDG_Spared    = SQDG - (MGDG + DGDG)/2,
    # LPC_PC         = LPC - PC,
    # LPE_PE         = LPE - PE,
    # 
    
    
    ## — Low-N contrasts
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
    
    ## — Membrane integrity
    PC_PE          = PC - PE,
    PC_PG          = PC - PG,
    PE_PG          = PE - PG,
    PC_total       = (PC + PE + PG)/3,
    PE_total       = (PC + PE + PG)/3,
    PG_total       = (PC + PE + PG)/3,
    Lyso_ratio     = (LPC + LPE)/2 - (PC + PE + PG)/3
  )

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 4)  LONG FORMAT + Welch t-tests                                 ║
# ╚══════════════════════════════════════════════════════════════════╝
ratio_long <- ratio_tbl %>%
  select(Sample, Condition,
         # choose the ratios you want – same as before
         DGDG_PC,DGDG_PE,MGDG_PC,MGDG_PE,
         PG_retention,PC_PE,Lyso_activity,nonP_P_total,

         # P supplementary
         #SQDG_Spared,LPC_PC,LPE_PE,
         
         TG_PC,DG_PC,TG_DG,
         DG_Phospho,TG_Phospho,stor_vs_photo,
         MGDG_DGDG,SQDG_DGDG,SQDG_MGDG,
         DGDG_total,MGDG_total,SQDG_total,
         PC_PE,PC_PG,PE_PG,Lyso_ratio,
         PC_total,PE_total,PG_total) %>%
  pivot_longer(-c(Sample, Condition),
               names_to = "Ratio", values_to = "Value") %>%
  filter(!is.na(Value))

stat_df <- ratio_long %>%
  group_by(Ratio) %>%
  summarise(p_value = t.test(Value ~ Condition)$p.value,
            y_pos   = max(Value, na.rm = TRUE) + 0.25,
            .groups = "drop") %>%
  mutate(star = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ "ns")) %>%
  filter(star != "ns")

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 5)  PLOTTING FUNCTION                                           ║
# ╚══════════════════════════════════════════════════════════════════╝
condition_cols <- c(Control = "#440154FF", LowInput = "#FDE725FF")

nature_theme <- theme_minimal(12) +
  theme(axis.line        = element_line(colour = "black"),
        axis.text        = element_text(colour = "black"),
        panel.grid.major = element_line(colour = "grey92"),
        panel.grid.minor = element_line(colour = "grey95"),
        legend.position  = "none")

make_plot <- function(ratio_vec, ncol_facets){
  df_sub  <- filter(ratio_long, Ratio %in% ratio_vec) %>%
    mutate(Ratio = factor(Ratio, levels = ratio_vec))
  stat_sub <- filter(stat_df, Ratio %in% ratio_vec) %>%
    mutate(Ratio = factor(Ratio, levels = ratio_vec))
  
  ggplot(df_sub, aes(Condition, Value, fill = Condition, colour = Condition)) +
    geom_violin(width = 0.9, alpha = 0.3, trim = TRUE, size = 0) +
    geom_boxplot(width = 0.6, alpha = 0.8, colour = "black", outlier.shape = NA) +
    facet_wrap(~ Ratio, ncol = ncol_facets, scales = "free_y") +
    geom_segment(data = stat_sub, inherit.aes = FALSE,
                 aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
                 colour = "black", size = 0.3) +
    geom_text(data = stat_sub, inherit.aes = FALSE,
              aes(x = 1.5, y = y_pos * 1.02, label = star),
              size = 5, colour = "black") +
    scale_x_discrete(labels = c("Control", "LowInput"),
                     expand  = expansion(add = 0.5)) +
    scale_fill_manual(values = condition_cols) +
    scale_colour_manual(values = condition_cols) +
    labs(x = NULL, y = "Z-scores") +   # <-- new y-axis label
    nature_theme +
    theme(strip.background = element_blank(),
          strip.text       = element_text(face = "bold", size = 10),
          panel.spacing    = unit(0.5, "lines"))
}

# groupings
ratios_lowP <- c("DGDG_PC","DGDG_PE","MGDG_PC","MGDG_PE",
                "PG_retention","PC_PE","Lyso_activity","nonP_P_total")


# Low P supplementary
#ratios_lowP <- c("SQDG_Spared","LPC_PC","LPE_PE")


ratios_lowN <- c("TG_PC","DG_PC","TG_DG","DG_Phospho","TG_Phospho","stor_vs_photo")
ratios_cold <- c("MGDG_DGDG","SQDG_DGDG","SQDG_MGDG",
                 "DGDG_total","MGDG_total","SQDG_total")
ratios_membrane <- c("PC_PE","PC_PG","PE_PG",
                     "PC_total","PE_total","PG_total","Lyso_ratio")

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 6)  DRAW & SAVE  FIGURES                                         ║
# ╚══════════════════════════════════════════════════════════════════╝
p_lowP <- make_plot(ratios_lowP, 4)
quartz();  p_lowP
ggsave("Fig2_lipid_ratio_log_lowP.png",  p_lowP, width = 12, height = 6,
       dpi = 300, bg = "white")

p_lowN <- make_plot(ratios_lowN, 3)
quartz();  p_lowN
ggsave("Fig3_lipid_ratio_log_lowN.png", p_lowN, width = 12, height = 6,
       dpi = 300, bg = "white")

p_cold <- make_plot(ratios_cold, 3)
quartz();  p_cold
ggsave("Fig4_lipid_ratio_log_cold.png", p_cold, width = 12, height = 6,
       dpi = 300, bg = "white")

p_memb <- make_plot(ratios_membrane, 3)
quartz();  p_memb
ggsave("SuppFig_lipid_ratio_log_membrane.png", p_memb, width = 12, height = 9,
       dpi = 300, bg = "white")




# ╔═══════════════════════════════════════════════════════════════════╗
# ║ 6)   MEDIAN ΔZ PER RATIO AND CONDITION                            ║
# ╚═══════════════════════════════════════════════════════════════════╝



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



]version
