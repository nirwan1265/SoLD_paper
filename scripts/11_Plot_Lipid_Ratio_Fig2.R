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
library(cowplot)
library(ggh4x)


# ╔══════════════════════════════════════════════════════════════════╗
# ║ 1) READ RAW INTENSITY TABLES AND CLASS                           ║
# ╚══════════════════════════════════════════════════════════════════╝

# Read the raw files 
control  <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
colnames(control)[1] <- "Compound_Name"  


lowinput  <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
colnames(lowinput)[1] <- "Compound_Name"  


# Valid classes
valid_classes <- c("TG","DG","MG",
                   "PC","PE","PI",
                   "DGDG","MGDG",
                   "SQDG","SM","AEG",
                   "LPC","LPE","PG","PA","PS")
class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 2) TIC normalisation + log10 per lipid                           ║
# ╚══════════════════════════════════════════════════════════════════╝

# Long format
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
# ║ 3)  PIVOT WIDE                                                   ║
# ╚══════════════════════════════════════════════════════════════════╝

# Wide format
wide_log <- long_all %>%
  pivot_wider(names_from = Class, values_from = class_log,
              values_fill = NA_real_)

#wide_log <- wide_log[,-3]
# ╔══════════════════════════════════════════════════════════════════╗
# ║ 4)  COMPUTE CLASS-PAIR                                           ║
# ╚══════════════════════════════════════════════════════════════════╝

# Ratio table
ratio_tbl <- wide_log %>%
  mutate(
    
    ### Membrane Lipid remodeling
    
    # ## 1. Sulfolipid Adjustments
    # #SQDG is the most vulnerable lipid:
    # SQDG_PC    = SQDG - PC,
    # SQDG_PG    = SQDG - PG,
    # SQDG_DGDG  = SQDG - DGDG, # Anionic lipid trade-offs
    # SQDG_MGDG  = SQDG - MGDG, 
    # SQDG_total = SQDG - ((DGDG + MGDG) / 2), # Sulfur allocation balance
    
    # VIP
    
    MG_SQDG    = MG - SQDG, # Monogalactosyl vs. sulfolipid balance
    PG_SQDG    = PG - SQDG, # Phosphatidylglycerol vs. sulfolipid balance
    PE_SQDG    = PE - SQDG, # Phosphatidylethanolamine vs. sulfolipid balance
    DGDG_SQDG  = DGDG - SQDG, # Digalactosyl vs. sulfolipid balance
    PC_SQDG    = PC - SQDG, # Phosphatidylcholine vs. sulfolipid balance
    SQDG_TG    = SQDG - TG, # Sulfolipid vs. storage lipid balance
    
    
    
    # ## 2. Galactolipid Dynamics
    # #MGDG is preferentially degraded:
    # MGDG_DGDG  = MGDG - DGDG, # Monogalactosyl vs. digalactosyl balance (bilayer vs. non-bilayer)
    # MGDG_PC    = MGDG - PC,
    # MGDG_total = MGDG - ((DGDG + SQDG) / 2), # Absolute monogalactolipid shift
    
    # VIP
    
    MG_MGDG  = MGDG - MG, # Monogalactosyl vs. monogalactolipid balance
    DGDG_MG  = DGDG - MGDG, # Digalactosyl vs. monogalactolipid balance
    
    
    # # DGDG is partially spared but outcompeted by phospholipids:
    # DGDG_total = DGDG - ((MGDG + SQDG) / 2), # Absolute digalactolipid shift
    # DGDG_PC    = DGDG - PC, #  Tests phospholipid replacement capacity under P stress
    # DGDG_PE    = DGDG - PE, #  Tests phospholipid replacement capacity under P stress
    # DGDG_PG    = DGDG - PG, #  Assess photosynthetic membrane integrity
    
    
    
    # ## 3. Phospholipid Homeostasis
    # PC_PE        = PC - PE, # Major bilayer asymmetry
    # PC_PG        = PC - PG,
    # PC_PS        = PC - PS, # Phosphatidylserine balance
    # PE_PS        = PE - PS, # Phosphatidylserine balance
    # PG_retention = PG - ((PC + PE) / 2), # Photosynthetic membrane priority
    # NonP_Phospho = (MGDG + DGDG + SQDG) - ((PC + PE) / 2), # Non-phospholipid balance
    # PC_total    = PC - ((PE + PG + PS) / 3), # Phospholipid balance
    # PE_total    = PE - ((PC + PG + PS) / 3), # Phospholipid balance
    # PS_total    = PS - ((PC + PE + PG) / 3), # Phospholipid balance
    # 
    
    #VIP
    
    PG_PS        = PG - PS, # Phosphatidylglycerol vs. phosphatidylserine balance
    PC_PS        = PC - PS, # Phosphatidylcholine vs. phosphatidylserine balance
    PE_PS        = PE - PS, # Phosphatidylethanolamine vs. phosphatidylserine balance
    PS_TG        = PS - TG, # Phosphatidylserine vs. storage lipid balance
    PA_PS        = PA - PS, # Phosphatidic acid vs. phosphatidylserine balance
    PG_retention = PG - ((PC + PE) / 2), # Photosynthetic membrane priority
    
    
    
    ### Lipid Turnover and Signaling
    
    ## 1. Diacylglycerol Flux Hub
    # DG Accumulates from Degraded Lipids but is Excluded from Membranes
    DG_DGDG    = DG   - DGDG,
    DG_MGDG    = DG   - MGDG, # Precursor flux to galactolipids
    
    DG_PC        = DG  - PC,
    DG_PE        = DG  - PE,
    DG_Phospho   = DG  - ((PC + PE) / 2), # Phospholipid-derived precursor pool
    DG_SQDG    = DG   - SQDG,
   
    
    
    ## 2. Lysophospholipids Remodeling
    LPC_PC       = LPC - PC,
    LPE_PE       = LPE - PE,
    Lyso_activity = (LPC + LPE) - (PC + PE), # Global phospholipase activity
    
    
    
    ### Carbon sink and allocation
    
    ## 1. Storage Lipid Reallocation
    
    # Triacylglycerol (TG) Synthesis Dominates Under Stress
    TG_DG        = TG  - DG, # DAG→TAG conversion efficiency
    TG_Galacto   = TG - (MGDG + DGDG)/2,  # Storage vs. chloroplast lipids
    
    
    
    TG_Phospho   = TG  - ((PC + PE) / 2), # Storage vs. membrane lipids
    DG_to_Gala   = DG - ((MGDG + DGDG)/2), # Precursor flux to galactolipids
    
    
    ## 2. Metabolic Trade-offs
    # Metabolic Trade-offs Favor Survival Over Growth
    Storage_vs_Photo    = (TG + DG)/2 - (MGDG + DGDG)/2,    # Storage vs. photosynthetic machinery
    Storage_vs_Membrane = (TG + DG)/2 - (PC + PE)/2,        # Carbon storage vs. membranes
    #Global_storage_vs_Membrane_balance = TG - ((PC + PE + PG) / 3)
    
  )


# ╔══════════════════════════════════════════════════════════════════╗
# ║ 5)  LONG FORMAT FOR RATIOS + Welch t-tests                       ║
# ╚══════════════════════════════════════════════════════════════════╝

# Long format ratio table
ratio_long <- ratio_tbl %>%
  dplyr::select(Sample,Condition,
                
                # # Sulfolipid Adjustments
                # SQDG_PC, SQDG_PG, SQDG_DGDG, SQDG_MGDG, SQDG_total,
                # 
                # # Galactolipid Dynamics
                # MGDG_DGDG,  MGDG_PC, MGDG_total, 
                # DGDG_total,  DGDG_PC, DGDG_PE, DGDG_PG,
                # 
                # # Phospholipid Homeostasis
                # PC_PE, PC_PG, PC_PS, PE_PS, PG_retention, NonP_Phospho, PC_total, PE_total, PS_total,
                
                #VIP
                MG_SQDG, PG_SQDG, PE_SQDG, DGDG_SQDG, PC_SQDG, SQDG_TG,
                # Galactolipid Dynamics
                MG_MGDG, DGDG_MG,
                # DGDG_total, DGDG_PC, DGDG_PE, DGDG_PG,
                # Phospholipid Homeostasis
                PG_PS, PC_PS, PE_PS, PS_TG, PA_PS, PG_retention,
                
                
                
                # Diacylglycerol Flux Hub
                DG_DGDG, DG_MGDG, 
                DG_PC, DG_PE,  DG_Phospho, DG_SQDG,
                
                # Lysophospholipids Remodeling 
                LPC_PC, LPE_PE, Lyso_activity,
                
                # Storage Lipid Reallocation
                TG_DG, TG_Galacto, 
                TG_Phospho,DG_to_Gala,
                
                # Metabolic Trade-offs
                Storage_vs_Photo, Storage_vs_Membrane
               

                
                                
  ) %>%
  pivot_longer(-c(Sample, Condition),
               names_to = "Ratio", values_to = "Value") %>%
  dplyr::filter(!is.na(Value) & !is.infinite(Value))


# T-test helper
stat_df <- ratio_long %>%
  dplyr::group_by(Ratio) %>%
  dplyr::summarise(p_value = t.test(Value ~ Condition)$p.value,
                   y_pos   = max(Value, na.rm = TRUE) + 0.25,
                   .groups = "drop") %>%
  dplyr::mutate(star = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ "ns")) %>%
  dplyr::filter(star != "ns")



# ╔══════════════════════════════════════════════════════════════════╗
# ║ 6)  RATIOS FOR EACH GROUP                                        ║
# ╚══════════════════════════════════════════════════════════════════╝

### Membrane Lipid remodeling
# mem_sulfolipid     <- c("SQDG_PC","SQDG_PG","SQDG_DGDG","SQDG_MGDG","SQDG_total")
# mem_galactolipid   <- c("MGDG_DGDG",  "MGDG_PC", "MGDG_total","DGDG_total",  "DGDG_PC", "DGDG_PE", "DGDG_PG" )
# mem_phospholipid   <- c("PC_PE","PC_PG","PC_PS", "PE_PS","PG_retention","NonP_Phospho","PC_total", "PE_total", "PS_total")

# VIP lipid ratios
mem_sulfolipid <- c("MG_SQDG","PG_SQDG","PE_SQDG","DGDG_SQDG","PC_SQDG","SQDG_TG")
mem_galactolipid <- c("MG_MGDG","DGDG_MG")
mem_phospholipid <- c("PG_PS","PC_PS","PE_PS","PS_TG","PA_PS")

### Lipid Turnover and Signaling
# turn_diacylglycerol<- c("DG_DGDG", "DG_MGDG","DG_PC","DG_PE","DG_Phospho","DG_SQDG")
# turn_lysophospho   <- c("LPC_PC","LPE_PE","Lyso_activity")

# VIP lipid ratios
turn_diacylglycerol <- c("DG_MG","DG_PG","DG_PE","DG_SQDG","DG_TG","LPC_TG")
turn_lysophospho <- c("LPC_MG","LPE_SQDG","LPC_LPE","LPE_MGDG","DG_LPE","LPC_PE","LPC_PG","LPC_TG")

### Carbon sink and allocation
# carbon_storage     <- c("TG_DG","TG_Galacto", "TG_Phospho","DG_to_Gala","")
# carbon_photo       <- c( "Storage_vs_Photo","Storage_vs_Membrane")

# VIP lipid ratios
carbon_storage     <- c("TG_DG","TG_Galacto", "TG_Phospho","DG_to_Gala")
carbon_photo       <- c( "Storage_vs_Photo","Storage_vs_Membrane")


# ╔══════════════════════════════════════════════════════════════════╗
# ║ 7)  PLOTTING FUNCTION                                            ║
# ╚══════════════════════════════════════════════════════════════════╝

# Colors
condition_cols <- c(Control = "#440154FF", LowInput = "#FDE725FF")

# Pre-used theme
nature_theme <- theme_minimal(12) +
  theme(axis.line        = element_line(colour = "black"),
        axis.text        = element_text(colour = "black"),
        panel.grid.major = element_line(colour = "grey92"),
        panel.grid.minor = element_line(colour = "grey95"),
        legend.position  = "none")

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 8)  PLOTTING                                                     ║
# ╚══════════════════════════════════════════════════════════════════╝

# Define a lookup, e.g. via a named vector or case_when:
ratio_to_group <- function(r) {
  if (r %in% mem_sulfolipid) {
    "Sulfolipid Adjustments"
    
  } else if (r %in% mem_galactolipid) {
    "Galactolipid Dynamics"
    
  } else if (r %in% mem_phospholipid) {
    "Phospholipid Homeostasis"
  } else {
    NA_character_
  }
}



##### Membrane Lipid Remodeling

desired_order <- c(mem_sulfolipid,mem_galactolipid,mem_phospholipid)
# Apply the function to create a new column in ratio_long
ratio_long2 <- ratio_long %>%
  mutate(Group = vapply(Ratio, ratio_to_group, character(1))) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Ratio = factor(Ratio, levels = desired_order),
    Group = factor(Group, levels = c(
    "Sulfolipid Adjustments",
    "Galactolipid Dynamics",
    "Phospholipid Homeostasis"
  )))

# Apply the function to create a new column in ratio_long
stat_df2 <- stat_df %>%
  mutate(Group = vapply(Ratio, ratio_to_group, character(1))) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Ratio = factor(Ratio, levels = desired_order),
    Group = factor(Group, levels = c(
    "Sulfolipid Adjustments",
    "Galactolipid Dynamics",
    "Phospholipid Homeostasis"
  )))

# Plotting
p_nested <- ggplot(ratio_long2, aes(Condition, Value, fill = Condition, colour = Condition)) +
  geom_violin(width = 0.9, alpha = 0.3, trim = TRUE) +
  geom_boxplot(width = 0.6, alpha = 0.8, colour = "black", outlier.shape = NA) +
  geom_segment(data = stat_df2, inherit.aes = FALSE,
               aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
               colour = "black", size = 0.3) +
  geom_text(data = stat_df2, inherit.aes = FALSE,
            aes(x = 1.5, y = y_pos * 1.02, label = star),
            size = 5, colour = "black") +
  scale_x_discrete(labels = c("Control", "LowInput"),
                   expand  = expansion(add = 0.5)) +
  scale_fill_manual(values = condition_cols) +
  scale_colour_manual(values = condition_cols) +
  labs(x = NULL, y = "Z-scores") +
  nature_theme +
  facet_grid2(
    rows = NULL,
    cols = vars(Group, Ratio),
    scales = "free_y",
    space = "free_x",
    strip = strip_nested(
      # Both levels use white background if desired, or NA for transparent:
      background_x = elem_list_rect(fill = c("white", "white")),
      # Both Group-level and Ratio-level text use plain face and same size:
      text_x = list(
        element_text(face = "plain", size = 10),  # Group-level
        element_text(face = "plain", size = 10)   # Ratio-level
      )
    )
  ) +
  theme(
    panel.spacing      = unit(0.5, "lines"),
    axis.line          = element_line(colour = "black"),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.major.x = element_blank(),
    legend.position    = "none",
    plot.title         = element_text(hjust = 0.5, face = "plain", size = 14)
  ) +
  ggtitle("Membrane Lipid Remodeling")

quartz(width = 12, height = 6)
print(p_nested)

ggsave("nested_lipid_remodeling_analysis2.png", p_nested,
       width = 16, height = 6, units = "in", dpi = 300, bg = "white")



##### Lipid Turnover and Signaling

# Define a lookup, e.g. via a named vector or case_when:
ratio_to_group2 <- function(r) {
  if (r %in% turn_diacylglycerol) {
    "Diacylglycerol Flux Hub"
    
  } else if (r %in% turn_lysophospho) {
    "Lysophospholipids Remodeling"
  } else {
    NA_character_
  }
}

desired_order <- c(turn_diacylglycerol,turn_lysophospho)

# Apply the function to create a new column in ratio_long
ratio_long3 <- ratio_long %>%
  mutate(Group = vapply(Ratio, ratio_to_group2, character(1))) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Ratio = factor(Ratio, levels = desired_order),
    Group = factor(Group, levels = c(
    "Diacylglycerol Flux Hub",
    "Lysophospholipids Remodeling"
  )))

# Apply the function to create a new column in stat_df
stat_df3 <- stat_df %>%
  mutate(Group = vapply(Ratio, ratio_to_group2, character(1))) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Ratio = factor(Ratio, levels = desired_order),
    Group = factor(Group, levels = c(
    "Diacylglycerol Flux Hub",
    "Lysophospholipids Remodeling"
  )))

# Plotting
p_nested2 <- ggplot(ratio_long3, aes(Condition, Value, fill = Condition, colour = Condition)) +
  geom_violin(width = 0.9, alpha = 0.3, trim = TRUE) +
  geom_boxplot(width = 0.6, alpha = 0.8, colour = "black", outlier.shape = NA) +
  geom_segment(data = stat_df3, inherit.aes = FALSE,
               aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
               colour = "black", size = 0.3) +
  geom_text(data = stat_df3, inherit.aes = FALSE,
            aes(x = 1.5, y = y_pos * 1.02, label = star),
            size = 5, colour = "black") +
  scale_x_discrete(labels = c("Control", "LowInput"),
                   expand  = expansion(add = 0.5)) +
  scale_fill_manual(values = condition_cols) +
  scale_colour_manual(values = condition_cols) +
  labs(x = NULL, y = "Z-scores") +
  nature_theme +
  facet_grid2(
    rows = NULL,
    cols   = vars(Group, Ratio),
    scales = "free_y",
    space  = "free_x",
    strip  = strip_nested(
      background_x   = elem_list_rect(fill   = c("white", "white")),
      text_x         = list(
        element_text(face   = "plain", size   = 10), # Group-level
        element_text(face   = "plain", size   = 10) # Ratio-level
      )
    )
  ) +
  theme(
    panel.spacing      = unit(0.5, "lines"),
    axis.line          = element_line(colour   = "black"),
    panel.grid.major.y = element_line(colour   = "grey90"),
    panel.grid.major.x = element_blank(),
    legend.position    = "none",
    plot.title         = element_text(hjust   = 0.5, face   ="plain", size   =14)
  ) +
  ggtitle("Lipid Turnover and Signaling")


quartz(width = 12, height = 6)
print(p_nested2)
ggsave("nested_lipid_turnover_analysis2.png", p_nested2,
       width = 12, height = 6, units = "in", dpi = 300, bg = "white")



##### Carbon sink and allocation


# Define a lookup, e.g. via a named vector or case_when:
ratio_to_group3 <- function(r) {
  if (r %in% carbon_storage) {
    "Storage Lipid Reallocation"
  } else if (r %in% carbon_photo) {
    "Metabolic Trade-offs"
  } else {
    NA_character_
  }
}

desired_order <- c(carbon_storage, carbon_photo)

# Apply the function to create a new column in ratio_long
ratio_long4 <- ratio_long %>%
  mutate(Group = vapply(Ratio, ratio_to_group3, character(1))) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Ratio = factor(Ratio, levels = desired_order),
    Group = factor(Group, levels = c(
    "Storage Lipid Reallocation",
    "Metabolic Trade-offs"
  )))
# Apply the function to create a new column in stat_df
stat_df4 <- stat_df %>%
  mutate(Group = vapply(Ratio, ratio_to_group3, character(1))) %>%
  filter(!is.na(Group)) %>%
  mutate(
    Ratio = factor(Ratio, levels = desired_order),
    Group = factor(Group, levels = c(
    "Storage Lipid Reallocation",
    "Metabolic Trade-offs"
  )))
# Plotting
p_nested3 <- ggplot(ratio_long4, aes(Condition, Value, fill = Condition, colour = Condition)) +
  geom_violin(width = 0.9, alpha = 0.3, trim = TRUE) +
  geom_boxplot(width = 0.6, alpha = 0.8, colour = "black", outlier.shape = NA) +
  geom_segment(data = stat_df4, inherit.aes = FALSE,
               aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
               colour = "black", size = 0.3) +
  geom_text(data = stat_df4, inherit.aes = FALSE,
            aes(x = 1.5, y = y_pos * 1.02, label = star),
            size = 5, colour = "black") +
  scale_x_discrete(labels = c("Control", "LowInput"),
                   expand  = expansion(add = 0.5)) +
  scale_fill_manual(values = condition_cols) +
  scale_colour_manual(values = condition_cols) +
  labs(x = NULL, y = "Z-scores") +
  nature_theme +
  facet_grid2(
    rows   = NULL,
    cols   = vars(Group, Ratio),
    scales = "free_y",
    space  = "free_x",
    strip  = strip_nested(
      background_x   = elem_list_rect(fill   = c("white", "white")),
      text_x         = list(
        element_text(face   ="plain", size   =10), # Group-level
        element_text(face   ="plain", size   =10) # Ratio-level
      )
    )
  ) +
  theme(
    panel.spacing      = unit(0.5, "lines"),
    axis.line          = element_line(colour   ="black"),
    panel.grid.major.y = element_line(colour   ="grey90"),
    panel.grid.major.x = element_blank(),
    legend.position    = "none",
    plot.title         = element_text(hjust   =0.5, face   ="plain", size   =14)
  ) +
  ggtitle("Carbon sink and allocation")

quartz(width = 12, height = 6)
print(p_nested3)
ggsave("nested_lipid_carbon_analysis2.png", p_nested3,
       width = 10, height = 6, units = "in", dpi = 300, bg = "white")








# OLD PLOT BY STRESS

# # ╔══════════════════════════════════════════════════════════════════╗
# # ║ 0)  PACKAGES                                                    ║
# # ╚══════════════════════════════════════════════════════════════════╝
# library(vroom)
# library(dplyr)
# library(tidyr)
# library(stringr)
# library(ggplot2)
# library(broom)      
# library(patchwork)  
# library(scales)
# library(cowplot)
# # ╔══════════════════════════════════════════════════════════════════╗
# # ║ 1)  READ RAW INTENSITY TABLES                                   ║
# # ╚══════════════════════════════════════════════════════════════════╝
# #ctrl_file <- "/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Control_all_lipids_final_non_normalized.csv"
# #lowp_file <- "/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv"
# 
# 
# ctrl_file <- "/Users/nirwantandukar/Documents/Github/SoLD_paper/results/spats_correction/control/control_all_lipids_fitted_phenotype_non_normalized.csv"
# 
# lowp_file <- "/Users/nirwantandukar/Documents/Github/SoLD_paper/results/spats_correction/lowinput/lowinput_all_lipids_fitted_phenotype_non_normalized.csv"
# 
# 
# 
# control  <- vroom(ctrl_file,  show_col_types = FALSE)
# lowinput <- vroom(lowp_file,  show_col_types = FALSE)
# 
# valid_classes <- c("TG","DG","MG",
#                    "PC","PE","PI",
#                    "DGDG","MGDG",
#                    "SQDG","SM","AEG",
#                    "LPC","LPE","PG","PA")
# class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")
# 
# # ╔══════════════════════════════════════════════════════════════════╗
# # ║ Modified long_prep: TIC normalisation + log10 per lipid (NO Z)  ║
# # ╚══════════════════════════════════════════════════════════════════╝
# long_prep_global <- function(df) {
#   df_long <- df %>%
#     pivot_longer(-c(Compound_Name, Condition),
#                  names_to = "Lipid", values_to = "Intensity") %>%
#     rename(Sample = Compound_Name)
#   
#   # 1) TIC per sample → relative abundance
#   df_long <- df_long %>%
#     group_by(Sample) %>%
#     mutate(TIC       = sum(Intensity, na.rm = TRUE),
#            rel_abund = Intensity / TIC) %>%
#     ungroup()
#   
#   # 2) small pseudo-count per sample → log10
#   df_long <- df_long %>%
#     group_by(Sample) %>%
#     mutate(minpos = ifelse(
#       all(is.na(rel_abund) | rel_abund <= 0),
#       NA_real_, min(rel_abund[rel_abund > 0], na.rm = TRUE)),
#       eps    = ifelse(is.na(minpos), 0, minpos * 0.5)) %>%
#     ungroup() %>%
#     mutate(log_rel = log10(rel_abund + eps))
#   
#   # 3) summarise mean log10(relative abundance) per lipid class
#   df_long <- df_long %>%
#     mutate(Class = str_extract(Lipid, class_pat)) %>%
#     filter(!is.na(Class)) %>%
#     group_by(Sample, Condition, Class) %>%
#     summarise(class_log = mean(log_rel, na.rm = TRUE), .groups = "drop")
#   
#   return(df_long)
# }
# 
# # combine Control + LowInput, then run prep
# combined  <- bind_rows(control  %>% mutate(Condition = "Control"),
#                        lowinput %>% mutate(Condition = "LowInput"))
# long_all  <- long_prep_global(combined)
# 
# # ╔══════════════════════════════════════════════════════════════════╗
# # ║ 2)  PIVOT WIDE                                                  ║
# # ╚══════════════════════════════════════════════════════════════════╝
# wide_log <- long_all %>%
#   pivot_wider(names_from = Class, values_from = class_log,
#               values_fill = NA_real_)
# 
# # ╔══════════════════════════════════════════════════════════════════╗
# # ║ 3)  COMPUTE CLASS-PAIR  log10-RATIOS                            ║
# # ╚══════════════════════════════════════════════════════════════════╝
# # ratio_tbl <- wide_log %>%
# #   mutate(
# #     ## — Low-P contrasts 
# #     DGDG_PC        = DGDG - PC,
# #     DGDG_PE        = DGDG - PE,
# #     MGDG_PC        = MGDG - PC,
# #     MGDG_PE        = MGDG - PE,
# # 
# #     PG_retention   = PG - (PC + PE)/2,
# # 
# #     PC_PE          = PC - PE,
# # 
# #     Lyso_activity  = (LPC + LPE) - (PC + PE),
# #     nonP_P_total   = (MGDG + DGDG) - (PC + PE),
# #     
# #     
# #     # Low P supplementary
# #     # SQDG_Spared    = SQDG - (MGDG + DGDG)/2,
# #     # LPC_PC         = LPC - PC,
# #     # LPE_PE         = LPE - PE,
# #     # 
# #     
# #     
# #     ## — Low-N contrasts
# #     TG_PC          = TG - PC,
# #     DG_PC          = DG - PC,
# #     TG_DG          = TG - DG,
# #     DG_Phospho     = DG - (PC + PE)/2,
# #     TG_Phospho     = TG - (PC + PE)/2,
# #     stor_vs_photo  = (TG + DG)/2 - (MGDG + PC)/2,
# #     
# #     ## — Cold contrasts
# #     MGDG_DGDG      = MGDG - DGDG,
# #     SQDG_DGDG      = SQDG - DGDG,
# #     SQDG_MGDG      = SQDG - MGDG,
# #     DGDG_total     = DGDG - ((MGDG + SQDG)/2),
# #     MGDG_total     = MGDG - ((DGDG + SQDG)/2),
# #     SQDG_total     = SQDG - ((DGDG + MGDG)/2),
# #     
# #     ## — Membrane integrity
# #     PC_PE          = PC - PE,
# #     PC_PG          = PC - PG,
# #     PE_PG          = PE - PG,
# #     PC_total       = (PC + PE + PG)/3,
# #     PE_total       = (PC + PE + PG)/3,
# #     PG_total       = (PC + PE + PG)/3,
# #     Lyso_ratio     = (LPC + LPE)/2 - (PC + PE + PG)/3
# #   )
# 
# 
# ratio_tbl <- wide_log %>%
#   mutate(
#     ### Membrane Lipid remodeling
#     
#     ## 1. Galactolipid Dynamics
#     DGDG_PC    = DGDG - PC, #  Tests phospholipid replacement capacity under P stress
#     DGDG_PE    = DGDG - PE, #  Tests phospholipid replacement capacity under P stress
#     MGDG_PC    = MGDG - PC,
#     #MGDG_PE    = MGDG - PE,
#     MGDG_DGDG  = MGDG - DGDG, # Monogalactosyl vs. digalactosyl balance (bilayer vs. non-bilayer)
#     DGDG_total = DGDG - ((MGDG + SQDG) / 2), # Absolute digalactolipid shift
#     MGDG_total = MGDG - ((DGDG + SQDG) / 2), # Absolute monogalactolipid shift
#     
#     ## 2. Phospholipid Homeostasis
#     PC_PE        = PC - PE, # Major bilayer asymmetry
#     PC_PG        = PC - PG,
#     PG_retention = PG - ((PC + PE) / 2), # Photosynthetic membrane priority
#     
#     ## 3. Sulfolipid Adjustments
#     SQDG_PC    = SQDG - PC,
#     SQDG_PG    = SQDG - PG,
#     SQDG_DGDG  = SQDG - DGDG, # Anionic lipid trade-offs
#     SQDG_MGDG  = SQDG - MGDG, 
#     SQDG_total = SQDG - ((DGDG + MGDG) / 2), # Sulfur allocation balance
#     
#     ### Lipid Turnover and Signaling
#     
#     ## 1. Lysophospholipids
#     LPC_PC       = LPC - PC,
#     LPE_PE       = LPE - PE,
#     Lyso_activity = (LPC + LPE) - (PC + PE), # Global phospholipase activity
#     
#     ## 2. Diacylglycerol_Flux Hub
#     DG_MGDG    = DG   - MGDG, # Precursor flux to galactolipids
#     DG_DGDG    = DG   - DGDG,
#     DG_SQDG    = DG   - SQDG,
#     DG_Phospho   = DG  - ((PC + PE) / 2), # Phospholipid-derived precursor pool
#     DG_PC        = DG  - PC,
#     DG_PE        = DG  - PE,
#     
#     
#     ### Carbon sink and allocation
#     
#     ## 1. Storage Lipid Reallocation
#     TG_DG        = TG  - DG, # DAG→TAG conversion efficiency
#     TG_Phospho   = TG  - ((PC + PE + PG) / 3),
#     
#     ## 2. Photosynthesis VS Storage
#     stor_vs_photo = ((TG + DG) / 2) - ((MGDG + PC) / 2) # Carbon sink vs structural investment
#     
#   )
# 
# 
# # ╔══════════════════════════════════════════════════════════════════╗
# # ║ 4)  LONG FORMAT + Welch t-tests                                 ║
# # ╚══════════════════════════════════════════════════════════════════╝
# # ratio_long <- ratio_tbl %>%
# #   dplyr::select(Sample, Condition,
# #          # choose the ratios you want – same as before
# #          DGDG_PC,DGDG_PE,MGDG_PC,MGDG_PE,
# #          PG_retention,PC_PE,Lyso_activity,nonP_P_total,
# # 
# #          # P supplementary
# #          #SQDG_Spared,LPC_PC,LPE_PE,
# #          
# #          TG_PC,DG_PC,TG_DG,
# #          DG_Phospho,TG_Phospho,stor_vs_photo,
# #          MGDG_DGDG,SQDG_DGDG,SQDG_MGDG,
# #          DGDG_total,MGDG_total,SQDG_total,
# #          PC_PE,PC_PG,PE_PG,Lyso_ratio,
# #          PC_total,PE_total,PG_total) %>%
# #   pivot_longer(-c(Sample, Condition),
# #                names_to = "Ratio", values_to = "Value") %>%
# #   dplyr::filter(!is.na(Value))
# 
# 
# ratio_long <- ratio_tbl %>%
#   dplyr::select(Sample,Condition,
#           # Galactolipid Dynamics
#           DGDG_PC, DGDG_PE, MGDG_PC, MGDG_DGDG,
#           DGDG_total, MGDG_total,
#           # Phospholipid Homeostasis
#           PC_PE, PC_PG, PG_retention,
#           # Sulfolipid Adjustments
#           SQDG_PC, SQDG_PG, SQDG_DGDG, SQDG_MGDG, SQDG_total,
#           # Lysophospholipids
#           LPC_PC, LPE_PE, Lyso_activity,
#           # Diacylglycerol_Flux Hub
#           DG_MGDG, DG_DGDG, DG_SQDG, DG_Phospho, DG_PC, DG_PE,
#           # Storage Lipid Reallocation
#           TG_DG, TG_Phospho,
#           # Carbon sink and allocation
#           stor_vs_photo
#           
#          
#                 ) %>%
#   pivot_longer(-c(Sample, Condition),
#                names_to = "Ratio", values_to = "Value") %>%
#   dplyr::filter(!is.na(Value) & !is.infinite(Value))
#                 
#     
# stat_df <- ratio_long %>%
#   dplyr::group_by(Ratio) %>%
#   dplyr::summarise(p_value = t.test(Value ~ Condition)$p.value,
#             y_pos   = max(Value, na.rm = TRUE) + 0.25,
#             .groups = "drop") %>%
#   dplyr::mutate(star = case_when(
#     p_value < 0.001 ~ "***",
#     p_value < 0.01  ~ "**",
#     p_value < 0.05  ~ "*",
#     TRUE            ~ "ns")) %>%
#   dplyr::filter(star != "ns")
# 
# # ╔══════════════════════════════════════════════════════════════════╗
# # ║ 5)  PLOTTING FUNCTION                                           ║
# # ╚══════════════════════════════════════════════════════════════════╝
# condition_cols <- c(Control = "#440154FF", LowInput = "#FDE725FF")
# 
# nature_theme <- theme_minimal(12) +
#   theme(axis.line        = element_line(colour = "black"),
#         axis.text        = element_text(colour = "black"),
#         panel.grid.major = element_line(colour = "grey92"),
#         panel.grid.minor = element_line(colour = "grey95"),
#         legend.position  = "none")
# 
# make_plot <- function(ratio_vec, ncol_facets){
#   df_sub  <- filter(ratio_long, Ratio %in% ratio_vec) %>%
#     mutate(Ratio = factor(Ratio, levels = ratio_vec))
#   stat_sub <- filter(stat_df, Ratio %in% ratio_vec) %>%
#     mutate(Ratio = factor(Ratio, levels = ratio_vec))
# 
#   ggplot(df_sub, aes(Condition, Value, fill = Condition, colour = Condition)) +
#     geom_violin(width = 0.9, alpha = 0.3, trim = TRUE, size = 0) +
#     geom_boxplot(width = 0.6, alpha = 0.8, colour = "black", outlier.shape = NA) +
#     facet_wrap(~ Ratio, ncol = ncol_facets, scales = "free_y") +
#     geom_segment(data = stat_sub, inherit.aes = FALSE,
#                  aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
#                  colour = "black", size = 0.3) +
#     geom_text(data = stat_sub, inherit.aes = FALSE,
#               aes(x = 1.5, y = y_pos * 1.02, label = star),
#               size = 5, colour = "black") +
#     scale_x_discrete(labels = c("Control", "LowInput"),
#                      expand  = expansion(add = 0.5)) +
#     scale_fill_manual(values = condition_cols) +
#     scale_colour_manual(values = condition_cols) +
#     labs(x = NULL, y = "Z-scores") +   # <-- new y-axis label
#     nature_theme +
#     theme(strip.background = element_blank(),
#           strip.text       = element_text(face = "bold", size = 10),
#           panel.spacing    = unit(0.5, "lines"))
# }
# 
# # groupings
# #ratios_lowP <- c("DGDG_PC","DGDG_PE","MGDG_PC","MGDG_PE",
#                 #"PG_retention","PC_PE","Lyso_activity","nonP_P_total")
# 
# 
# # Low P supplementary
# #ratios_lowP <- c("SQDG_Spared","LPC_PC","LPE_PE")
# 
# 
# #ratios_lowN <- c("TG_PC","DG_PC","TG_DG","DG_Phospho","TG_Phospho","stor_vs_photo")
# #ratios_cold <- c("MGDG_DGDG","SQDG_DGDG","SQDG_MGDG",
# #                 "DGDG_total","MGDG_total","SQDG_total")
# #ratios_membrane <- c("PC_PE","PC_PG","PE_PG",
# #                     "PC_total","PE_total","PG_total","Lyso_ratio")
# 
# 
# 
# # 1) Define ratio vectors for each subgroup
# 
# ### Membrane Lipid remodeling
# mem_galactolipid   <- c("DGDG_PC","DGDG_PE","MGDG_PC","MGDG_DGDG","DGDG_total","MGDG_total")
# mem_phospholipid   <- c("PC_PE","PC_PG","PG_retention")
# mem_sulfolipid     <- c("SQDG_PC","SQDG_PG","SQDG_DGDG","SQDG_MGDG","SQDG_total")
# 
# 
# ### Lipid Turnover and Signaling
# turn_lysophospho   <- c("LPC_PC","LPE_PE","Lyso_activity")
# turn_diacylglycerol<- c("DG_MGDG","DG_DGDG","DG_SQDG","DG_Phospho","DG_PC","DG_PE")
# 
# 
# ### Carbon sink and allocation
# carbon_storage     <- c("TG_DG","TG_Phospho")
# carbon_photo       <- c("stor_vs_photo")
# 
# 
# 
# # ╔══════════════════════════════════════════════════════════════════╗
# # ║ 6)  DRAW & SAVE  FIGURES                                         ║
# # ╚══════════════════════════════════════════════════════════════════╝
# 
# ### Membrane Lipid remodeling
# 
# # Galactolipid homeostasis
# galactolipid_homeostasis_plot <- make_plot(mem_galactolipid, 3)
# quartz();  galactolipid_homeostasis_plot
# ggsave("Fig2_lipid_ratio_log_galactolipid_homeostasis.png", 
#        galactolipid_homeostasis_plot, width = 12, height = 6,
#        dpi = 300, bg = "white")
# 
# # Phospholipid dynamics
# phospholipid_dynamics_plot <- make_plot(mem_phospholipid, 3)
# quartz();  phospholipid_dynamics_plot
# ggsave("Fig3_lipid_ratio_log_phospholipid_dynamics.png", 
#        phospholipid_dynamics_plot, width = 12, height = 6,
#        dpi = 300, bg = "white")
# 
# # Sulfolipid adjustments
# sulfolipid_adjustments_plot <- make_plot(mem_sulfolipid, 3)
# quartz();  sulfolipid_adjustments_plot
# ggsave("Fig5_lipid_ratio_log_sulfolipid_adjustments.png", 
#        sulfolipid_adjustments_plot, width = 12, height = 6,
#        dpi = 300, bg = "white")
# 
# 
# ### Lipid Turnover and Signaling
# 
# # Membrane Hyrolysis
# turn_lysophospho_plot <- make_plot(turn_lysophospho, 3)
# quartz();  turn_lysophospho_plot
# ggsave("Fig3_lipid_ratio_log_turn_lysophospho.png", 
#        turn_lysophospho_plot, width = 12, height = 6,
#        dpi = 300, bg = "white")
# 
# # Diacylglycerol Flux Hub
# turn_diacylglycerol_plot <- make_plot(turn_diacylglycerol, 3)
# quartz();  turn_diacylglycerol_plot
# 
# ggsave("Fig4_lipid_ratio_log_turn_diacylglycerol.png", 
#        turn_diacylglycerol_plot, width = 12, height = 6,
#        dpi = 300, bg = "white")
# 
# 
# ### Carbon sink and allocation
# 
# # Carbon storage reallocation
# carbon_storage_plot <- make_plot(carbon_storage, 3)
# 
# quartz();  carbon_storage_plot
# ggsave("Fig5_lipid_ratio_log_carbon_storage.png", 
#        carbon_storage_plot, width = 12, height = 6,
#        dpi = 300, bg = "white")
# 
# # Carbon sink vs photosynthesis
# carbon_photo_plot <- make_plot(carbon_photo, 3)
# quartz();  carbon_photo_plot
# ggsave("Fig6_lipid_ratio_log_carbon_photo.png", 
#        carbon_photo_plot, width = 12, height = 6,
#        dpi = 300, bg = "white")
# 
# 
# 
# 
# 
# 
# 
# library(dplyr)
# library(ggplot2)
# 
# # Assume ratio_long is your tibble: Sample, Condition, Ratio, Value
# # Assume stat_df is your tibble: Ratio, p_value, y_pos, star
# # Assume condition_cols exists, e.g.
# condition_cols <- c(Control = "#440154FF", LowInput = "#FDE725FF")
# 
# # 1. Define subclass mapping:
# ratio_long2 <- ratio_long %>%
#   mutate(
#     Subclass = case_when(
#       Ratio %in% mem_galactolipid  ~ "1. Galactolipid Dynamics",
#       Ratio %in% mem_phospholipid  ~ "2. Phospholipid Homeostasis",
#       Ratio %in% mem_sulfolipid    ~ "3. Sulfolipid Adjustments",
#       TRUE                         ~ NA_character_
#     ),
#     Condition = factor(Condition, levels = c("Control","LowInput"))
#   ) %>%
#   filter(!is.na(Subclass))  # only keep the three groups
# 
# stat_df2 <- stat_df %>%
#   filter(Ratio %in% c(mem_galactolipid, mem_phospholipid, mem_sulfolipid)) %>%
#   mutate(
#     Subclass = case_when(
#       Ratio %in% mem_galactolipid  ~ "1. Galactolipid Dynamics",
#       Ratio %in% mem_phospholipid  ~ "2. Phospholipid Homeostasis",
#       Ratio %in% mem_sulfolipid    ~ "3. Sulfolipid Adjustments",
#       TRUE                         ~ NA_character_
#     )
#   ) %>%
#   filter(!is.na(Subclass))
# 
# # 2. Ensure Ratio is a factor in the desired order within each subclass
# # We build a single vector in the desired order:
# ratio_order <- c(mem_galactolipid, mem_phospholipid, mem_sulfolipid)
# ratio_long2 <- ratio_long2 %>%
#   mutate(Ratio = factor(Ratio, levels = ratio_order),
#          Subclass = factor(Subclass, levels = c(
#            "1. Galactolipid Dynamics",
#            "2. Phospholipid Homeostasis",
#            "3. Sulfolipid Adjustments"
#          )))
# 
# stat_df2 <- stat_df2 %>%
#   mutate(Ratio = factor(Ratio, levels = ratio_order),
#          Subclass = factor(Subclass, levels = c(
#            "1. Galactolipid Dynamics",
#            "2. Phospholipid Homeostasis",
#            "3. Sulfolipid Adjustments"
#          )))
# 
# # 3. Build the plot
# p <- ggplot(ratio_long2, aes(x = Condition, y = Value, fill = Condition, colour = Condition)) +
#   geom_violin(width = 0.9, alpha = 0.3, trim = TRUE) +
#   geom_boxplot(width = 0.6, alpha = 0.8, colour = "black", outlier.shape = NA) +
#   # add significance segments & stars:
#   geom_segment(data = stat_df2,
#                inherit.aes = FALSE,
#                aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
#                colour = "black", size = 0.3) +
#   geom_text(data = stat_df2,
#             inherit.aes = FALSE,
#             aes(x = 1.5, y = y_pos * 1.02, label = star),
#             size = 4, colour = "black") +
#   # facet_grid: rows = Subclass, cols = Ratio
#   facet_grid(Subclass ~ Ratio,
#              scales = "free_y",
#              space = "free_x",
#              switch = "y") +
#   scale_x_discrete(labels = c("Control","LowInput"), expand = expansion(add = 0.5)) +
#   scale_fill_manual(values = condition_cols) +
#   scale_colour_manual(values = condition_cols) +
#   labs(x = NULL, y = "Z-scores") +
#   theme_minimal(base_size = 11) +
#   theme(
#     # Move the row strips (Subclass) to the left outside the panels:
#     strip.placement = "outside",
#     # For column strips (Ratio), keep at top:
#     strip.text.x = element_text(face = "bold", size = 10),
#     # For row strips, rotate text horizontally:
#     strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
#     # Grey background for the row strips:
#     strip.background = element_rect(fill = "grey80", colour = "black"),
#     # Panel spacing
#     panel.spacing = unit(0.5, "lines"),
#     # Keep axis lines
#     panel.border = element_blank(),
#     axis.line = element_line(colour = "black"),
#     # Remove individual panel border boxes (we rely on the strip grouping)
#     panel.grid.major.y = element_line(colour = "grey90"),
#     panel.grid.major.x = element_blank(),
#     legend.position = "none"
#   )
# 
# # 4. Add an overall title above:
# #    We can draw the title in cowplot::plot_grid by stacking a blank with title
# library(cowplot)
# title_grob <- ggdraw() + 
#   draw_label("Membrane Lipid Remodeling", fontface="bold", size=14)
# 
# final <- plot_grid(
#   title_grob,
#   p + theme(plot.margin = margin(t = 0, r = 5, b = 5, l = 5)),
#   ncol = 1,
#   rel_heights = c(0.05, 1)
# )
# 
# # 5. Display and save
# quartz()
# print(final)
# ggsave("Fig_lipid_ratios_by_subclass.png", final, width=12, height=8, dpi=300, bg="white")
# 
# 
# 
# 
# 
# 
# # ╔══════════════════════════════════════════════════════════════════╗
# # ║ 6)  DRAW & SAVE  FIGURES                                         ║
# # ╚══════════════════════════════════════════════════════════════════╝
# # p_lowP <- make_plot(ratios_lowP, 4)
# # quartz();  p_lowP
# # ggsave("Fig2_lipid_ratio_log_lowP.png",  p_lowP, width = 12, height = 6,
# #        dpi = 300, bg = "white")
# # 
# # p_lowN <- make_plot(ratios_lowN, 3)
# # quartz();  p_lowN
# # ggsave("Fig3_lipid_ratio_log_lowN.png", p_lowN, width = 12, height = 6,
# #        dpi = 300, bg = "white")
# # 
# # p_cold <- make_plot(ratios_cold, 3)
# # quartz();  p_cold
# # ggsave("Fig4_lipid_ratio_log_cold.png", p_cold, width = 12, height = 6,
# #        dpi = 300, bg = "white")
# # 
# # p_memb <- make_plot(ratios_membrane, 3)
# # quartz();  p_memb
# # ggsave("SuppFig_lipid_ratio_log_membrane.png", p_memb, width = 12, height = 9,
# #        dpi = 300, bg = "white")
# # 
# 
# galactolipid_homeostasis_plot <- make_plot(galactolipid_homeostasis, 3)
# quartz();  galactolipid_homeostasis_plot
# 
# ggsave("Fig2_lipid_ratio_log_galactolipid_homeostasis.png", 
#        galactolipid_homeostasis_plot, width = 12, height = 6,
#        dpi = 300, bg = "white")
# 
# phospholipid_dynamics_plot <- make_plot(phospholipid_dynamics, 3)
# quartz();  phospholipid_dynamics_plot
# 
# ggsave("Fig3_lipid_ratio_log_phospholipid_dynamics.png", 
#        phospholipid_dynamics_plot, width = 12, height = 6,
#        dpi = 300, bg = "white")
# 
# storage_lipid_reallocation_plot <- make_plot(storage_lipid_reallocation, 3)
# quartz();  storage_lipid_reallocation_plot
# 
# ggsave("Fig4_lipid_ratio_log_storage_lipid_reallocation.png", 
#        storage_lipid_reallocation_plot, width = 12, height = 6,
#        dpi = 300, bg = "white")
# 
# sulfolipid_adjustments_plot <- make_plot(sulfolipid_adjustments, 3)
# quartz();  sulfolipid_adjustments_plot
# 
# ggsave("Fig5_lipid_ratio_log_sulfolipid_adjustments.png", 
#        sulfolipid_adjustments_plot, width = 12, height = 6,
#        dpi = 300, bg = "white")
# 
# # ╔═══════════════════════════════════════════════════════════════════╗
# # ║ 6)   MEDIAN ΔZ PER RATIO AND CONDITION                            ║
# # ╚═══════════════════════════════════════════════════════════════════╝
# 
# 
# 
# library(dplyr)
# library(tidyr)
# 
# # 1a) Compute median ΔZ per Ratio and Condition
# median_df <- ratio_long %>%
#   group_by(Ratio, Condition) %>%
#   summarise(
#     medianZ = median(Value, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # 1b) Pivot so Control and LowInput are side by side, then compute the difference
# median_diff <- median_df %>%
#   pivot_wider(names_from = Condition, values_from = medianZ) %>%
#   # Adjust the names “Control” and “LowInput” if your labels differ
#   mutate(
#     delta_median = LowInput - Control
#   )
# 
# # 1c) Inspect only the Low-P ratios (replace ratios_lowP with your vector of names)
# median_diff_lowP <- median_diff %>%
#   filter(Ratio %in% ratios_lowP) %>%
#   arrange(desc(delta_median))
# 
# print(median_diff_lowP, n = Inf)
# 
