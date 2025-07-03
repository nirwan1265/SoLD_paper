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
library(ropls)

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 1) READ RAW INTENSITY TABLES AND CLASS                           ║
# ╚══════════════════════════════════════════════════════════════════╝

# Read the raw files 
control  <- vroom("data/SPATS_fitted/non_normalized_intensities/control_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
colnames(control)[1] <- "Compound_Name"  

lowinput  <- vroom("data/SPATS_fitted/non_normalized_intensities/lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
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

# Removing  AEG and PG 
wide_log <- wide_log[,-c(3,13)]
wide_log <- wide_log[,-c(3)]




# 1) pivot into “long” so we have one row per Sample/Condition/Class
long_classes <- wide_log %>%
  pivot_longer(
    cols      = -c(Sample, Condition),
    names_to  = "Class",
    values_to = "Abundance"
  )

# 2) self-join on Sample+Condition, then keep only one direction of each pair
ratios_long <- long_classes %>%
  inner_join(long_classes, 
             by     = c("Sample","Condition"),
             suffix = c(".num", ".den")) %>%
  # only keep each unordered pair once
  filter(Class.num < Class.den) %>%
  transmute(
    Sample,
    Condition,
    RatioName = paste(Class.num, Class.den, sep = "/"),
    Ratio     = Abundance.num / Abundance.den
  )

# 3) inspect    
ratios_long

# Remove rows with string count less than 3 from Sample
ratios_long <- ratios_long %>% 
  filter(nchar(Sample) >= 5)

table(ratios_long$Sample)


# find the “full” count
full_n <- ratios_long %>% 
  group_by(Sample) %>% 
  tally() %>% 
  pull(n) %>% 
  max()

# keep only those with the full complement
ratios_clean <- ratios_long %>% 
  group_by(Sample) %>% 
  filter(n() == full_n) %>% 
  ungroup()


# 1) Pivot your ratios into a (n_samples × n_features) matrix
wide_ratios <- ratios_clean %>% 
  # make one row per Sample×Condition
  pivot_wider(
    id_cols   = c(Sample, Condition),
    names_from = RatioName,
    values_from = Ratio
  ) %>% 
  arrange(Sample, Condition)

# keep the metadata
meta <- wide_ratios %>% dplyr::select(Sample, Condition)

# pull out the numeric feature matrix
lipid_mat <- wide_ratios %>% 
  dplyr::select(-Sample, -Condition) %>% 
  as.matrix()

# 2) build your class vector: a factor of your two Conditions
class_vec <- factor(meta$Condition)

# sanity check
stopifnot(nrow(lipid_mat) == length(class_vec))

# 3) fit your OPLS
set.seed(123)  # for reproducibility of any CV‐based choices
quartz()
opls_mod <- opls(
  lipid_mat, 
  class_vec,
  predI  = 1,       # one predictive component
  orthoI = NA,      # have opls choose orthogonal comps by CV
  scaleC = "standard"  # mean‐center & unit‐variance scale
)

# 4) take a look
print(opls_mod)
plot(opls_mod)       # scores, loadings, CV‐metrics, etc.



# xtract VIP scores
vip_vec <- slot(opls_mod, "vipVn")     # numeric vector with names = lipid IDs

vip_df <- tibble(
  Lipid = names(vip_vec),
  VIP   = as.numeric(vip_vec)
)

# 2.  Keep VIP > 1.3  (and sort descending)
vip_hits <- vip_df %>%
  filter(VIP > 1) %>%
  arrange(desc(VIP))

# How many?
cat("Number of discriminatory lipids (VIP > 1.3):", nrow(vip_hits), "\n")
print(vip_hits, n = Inf)

pg_cv <- raw_peak_tbl %>%     # your raw species-level table
  filter(Class == "PG") %>% 
  group_by(SampleGroup) %>%   # QC or pooled injections
  summarise(CV = sd(Intensity) / mean(Intensity)) 


quartz()
ggplot(vip_hits, aes(x = reorder(Lipid, VIP), y = VIP)) +
  geom_col(
    width    = 0.7,
    fill     = "grey50",
    colour   = "black",
    linewidth = 0.2
  ) +
  coord_flip() +
  geom_hline(
    yintercept = 1.3,
    linetype   = "dashed",
    colour     = "red"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x     = NULL,
    y     = "VIP score",
    title = "Lipids driving Control vs Low-input separation\n(OPLS-DA, VIP > 1.3)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.line          = element_line(colour = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90"),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    plot.title         = element_text(hjust = 0.5, face = "plain", size = 10)
  )



# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ################################################################################
# #### SOrghum Lipidomics Database (SOLD)
# ################################################################################
# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 
# 
# ###############################################################################
# ## Supplementary Figures – Overlap of Control vs Low‑input lipid lists
# ## Three panels:
# ##   S1 – All detected lipids               (Venn)
# ##   S2 – Traditional classes only          (Venn)
# ##   S3 – Curated “Class” annotation set    (table + bar *optional*)
# ###############################################################################
# 
# # ------------------------------------------------------------------------------
# # 1. Load libraries
# # ------------------------------------------------------------------------------
# library(vroom)      # fast CSV import
# library(dplyr)      # data wrangling
# library(tidyr)      # pivoting
# library(ggplot2)    # plotting
# library(stringr)    # string manipulation
# library(forcats)    # factor reordering
# library(tidyverse)
# library(ggpubr)      # For stat_compare_means
# library(scales) 
# library(ggVennDiagram)   # nice, ggplot‑based 2‑set Venn
# library(readr) 
# library(sf)
# library(ggVennDiagram) # For Venn diagrams)
# library(viridis)
# library(ropls)    # OPLS‑DA
# library(limma)
# library(knitr)
# library(kableExtra)
# 
# # ------------------------------------------------------------------------------
# # 2. Read in the data
# # ------------------------------------------------------------------------------
# #control  <- vroom("data/SPATS_fitted/BLUP_GWAS_phenotype/control_all_lipids_BLUPs_rename.csv")
# control  <- vroom("data/SPATS_fitted/non_normalized_intensities/control_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
# colnames(control)[1] <- "Compound_Name"  
# 
# # Select only partial string match to columns with strings TG, DG, PC, PE, PA, PS, DGDG, SQDG, MGDG, MG
# control  <- control %>% dplyr::select(Compound_Name, contains("TG("), contains("DG("), contains("PC("),
#                                       contains("PE("), contains("PA("), contains("PS("),
#                                       contains("DGDG("), contains("SQDG("), contains("MGDG("), contains("MG("))
# 
# # Remove rows with strin gcount less than 5 from Compound_Name
# control <- control %>% 
#   filter(nchar(Compound_Name) >= 5)
# 
# # Lowinput
# #lowinput  <- vroom("data/SPATS_fitted/BLUP_GWAS_phenotype/lowinput_all_lipids_BLUPs_rename.csv")
# lowinput  <- vroom("data/SPATS_fitted/non_normalized_intensities/lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
# colnames(lowinput)[1] <- "Compound_Name"  
# 
# 
# # Select only partial string match to columns with strings TG, DG, PC, PE, PA, PS, DGDG, SQDG, MGDG, MG
# lowinput  <- lowinput %>% dplyr::select(Compound_Name, contains("TG("), contains("DG("), contains("PC("),
#                                           contains("PE("), contains("PA("), contains("PS("),
#                                           contains("DGDG("), contains("SQDG("), contains("MGDG("), contains("MG("))
# 
# 
# # Remove rows with strin gcount less than 5 from Compound_Name
# lowinput <- lowinput %>% 
#   filter(nchar(Compound_Name) >= 5)
# 
# 
# 
# #### CHANGE THE NAME TO COMMONNAME
# # Lipid class
# # lipid_class_info <- vroom::vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/data/lipid_class.csv", show_col_types = FALSE) %>%
# #   dplyr::filter(!is.na(Class))
# lipid_class_info <- vroom::vroom("data/lipid_class/lipid_classes.csv", show_col_types = FALSE) %>% dplyr::distinct(Lipids, .keep_all = TRUE)  # keep unique Lipids with their CommonName
# 
# # Create a named vector of replacements: names = original, values = new names
# name_map <- lipid_class_info %>%
#   dplyr::filter(!is.na(CommonName)) %>%
#   dplyr::select(Lipids, CommonName) %>%
#   deframe()  # turns it into a named vector: Lipids -> CommonName
# 
# # Function to rename lipid columns in any dataset
# rename_lipid_columns <- function(df, name_map) {
#   original_names <- colnames(df)
#   
#   renamed <- original_names %>%
#     # Replace only if a match exists
#     sapply(function(x) if (x %in% names(name_map)) name_map[[x]] else x)
#   
#   # Apply renamed vector back
#   colnames(df) <- renamed
#   return(df)
# }
# 
# # Apply to both control and lowinput
# control  <- rename_lipid_columns(control,  name_map)
# lowinput <- rename_lipid_columns(lowinput, name_map)
# 
# 
# 
# # 0.  Ensure identical lipid columns in both data frames
# shared_lipids <- intersect(names(control)[-1],   # exclude Compound_Name
#                            names(lowinput)[-1])
# 
# ctrl  <- control  %>% dplyr::select(Compound_Name, all_of(shared_lipids))
# lowin <- lowinput %>% dplyr::select(Compound_Name, all_of(shared_lipids))
# 
# 
# # 1.  Assemble one long table with a Condition flag
# full_tbl <- bind_rows(
#   ctrl  %>% mutate(Condition = "Control"),
#   lowin %>% mutate(Condition = "LowInput")
# )
# 
# 
# # 2.  Build lipid_mat  (samples × lipids)  and  class_vec
# 
# lipid_mat <- full_tbl %>%
#   dplyr::select(-Compound_Name, -Condition) %>%   # leave only numeric columns
#   as.matrix()
# 
# class_vec <- factor(full_tbl$Condition, levels = c("Control", "LowInput"))
# 
# 
# # 3.  (Optional) Pareto‑scale or autoscale before OPLS‑DA
# #  -- Pareto scaling: mean‑centre, divide by sqrt(sd)
# lipid_mat <- scale(lipid_mat, center = TRUE, scale = sqrt(apply(lipid_mat, 2, sd)))
# 
# # Check dimensions
# dim(lipid_mat)   #  ➜  samples × lipids
# table(class_vec) #  ➜  sample counts per condition
# lipid_mat[1:5,1:5]
# str(lipid_mat)
# 
# # 1. Pull out only the numeric matrix
# #    (drop any non-numeric columns if you have e.g. lipid names in the first column)
# lipid_numeric <- as.matrix(lipid_mat)  
# lipid_numeric <- t(lipid_numeric)
# dim(lipid_numeric)
# # if lipid_mat has a "Compound_Name" column, do:
# # lipid_numeric <- as.matrix(lipid_mat %>% select(-Compound_Name))
# 
# 
# # 2. Build a simple factor/vector of batch labels
# batch_vec <- c(
#   rep("Control", 394),
#   rep("LowInput", 363)
# )
# 
# # check it matches the number of columns
# stopifnot(length(batch_vec) == ncol(lipid_numeric))
# 
# # 2.  Run OPLS‑DA  -----------------------------------------------------
# set.seed(123)
# quartz()
# opls_mod <- opls(lipid_mat,
#                  class_vec,
#                  predI    = 1,       # 1 predictive component (default)
#                  orthoI   = NA,      # ropls chooses # orthogonal comps via CV
#                  scaleC   = "standard")  # autoscale (mean‑center, unit variance)
# 
# # 3.  Model diagnostics  ----------------------------------------------
# summary(opls_mod)      # R2X, R2Y, Q2, # components
# 
# quartz()
# plot(opls_mod, type = "overview")      # score + loading plots
# quartz()
# plot(opls_mod, type = "correlation")   # VIP vs correlation
# 
# # 4.  Permutation test (over‑fit check) -------------------------------
# perm_res <- opls(lipid_mat, class_vec,
#                  permI = 200,         # 200 permuted models
#                  predI = 1)
# 
# # look at p‑value in perm_res@summaryDF
# perm_res@summaryDF
# 
# 
# # 1.  Extract VIP scores
# 
# vip_vec <- slot(opls_mod, "vipVn")     # numeric vector with names = lipid IDs
# 
# vip_df <- tibble(
#   Lipid = names(vip_vec),
#   VIP   = as.numeric(vip_vec)
# )
# 
# # 2.  Keep VIP > 1.3  (and sort descending)
# 
# vip_hits <- vip_df %>%
#   filter(VIP > 1.2) %>%
#   arrange(desc(VIP))
# 
# # How many?
# cat("Number of discriminatory lipids (VIP > 1.3):", nrow(vip_hits), "\n")
# 
# 
# # 3.  Annotate with lipid class
# lookup_tbl <-
#   lipid_class_info %>%
#   dplyr::select(Lipid = Lipids, SubClass) %>%
#   mutate(Lipid = str_trim(Lipid)) %>%
#   distinct(Lipid, .keep_all = TRUE)
# 
# 
# # Compare with vip hits
# vip_hits <- vip_hits %>%
#   left_join(lookup_tbl, by = "Lipid") %>%          # add Class column
#   mutate(SubClass = ifelse(is.na(SubClass), "Unk", SubClass))
# 
# 
# # 4A.  Save as Supplementary Table (CSV)
# write.csv(vip_hits,
#           "VIP_lipids_OPLSDA_Control_vs_LowInput.csv",
#           row.names = FALSE)
# 
# # 4B.  Publication bar plot  (top 30 VIPs)
# vip_top30 <- vip_hits %>% slice_head(n = 30)
# 
# quartz()
# ggplot(vip_top30,
#        aes(x = reorder(Lipid, VIP),
#            y = VIP,
#            fill = SubClass)) +
#   geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
#   coord_flip() +
#   geom_hline(yintercept = 1.3, linetype = "dashed", colour = "red") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
#   scale_fill_viridis_d(option = "D", direction = 1) +
#   labs(
#     x     = NULL,
#     y     = "VIP score",
#     title = "Top 30 lipids driving Control vs Low-input separation\n(OPLS-DA, VIP > 1.3)"
#   ) +
#   theme_bw(base_size = 9) +
#   theme(
#     axis.line          = element_line(colour = "black"),
#     panel.grid.major.y = element_blank(),
#     panel.grid.major.x = element_line(colour = "grey90"),
#     panel.grid.minor   = element_blank(),
#     legend.position    = "right",
#     legend.title       = element_blank(),
#     plot.title         = element_text(hjust = 0.5, face = "plain", size = 10)
#   )










