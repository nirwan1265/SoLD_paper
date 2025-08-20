#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
#### SOrghum Lipidomics Database (SOLD)
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

################################################################################
###############################################################################
# 1.  LIBRARIES
###############################################################################

suppressPackageStartupMessages({
  library(vroom);   library(dplyr);   library(tidyr);  library(stringr)
  library(ggplot2); library(viridis); library(gridExtra); library(grid); library(purrr)
})


################################################################################
# 1.  LOAD DATA  (Control & LowInput)
################################################################################
#control  <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Control_all_lipids_final_non_normalized.csv")
#lowinput <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv")
# control <- vroom("/Users/nirwantandukar/Documents/Github/SoLD_paper/results/spats_correction/control/control_all_lipids_fitted_phenotype_non_normalized.csv")

control  <- vroom("data/SPATS_fitted/BLUP_GWAS_phenotype/control_all_lipids_BLUPs.csv") %>% dplyr::select(-c(2,3,4))
colnames(control)[1] <- "Compound_Name"  

lowinput  <- vroom("data/SPATS_fitted/BLUP_GWAS_phenotype/lowinput_all_lipids_BLUPs.csv") %>% dplyr::select(-c(2,3,4))
colnames(lowinput)[1] <- "Compound_Name"  


################################################################################
# 2.  SELECT THE CLASSES 
################################################################################

valid_classes <- c("TG","DG","MG","PC","PE","PG","PI",
                   "LPC","LPE","DGDG","MGDG","Cer","GalCer","SM","FA","SQDG","AEG","PA","PS")

#valid_classes <- c("LPE","LPC")
#valid_classes <- c("PC","PE","PA","PG","PS")
#valid_classes <- c("TG","DG","MG","PC","PE",
#                   "DGDG","MGDG","SQDG")

class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

reshape_plate <- function(df, label){
  df %>%
    pivot_longer(-Compound_Name, names_to = "Lipid", values_to = "Intensity") %>%
    mutate(Class = str_extract(Lipid, class_pat)) %>%
    filter(!is.na(Class)) %>%
    group_by(Sample = Compound_Name, Class) %>%
    summarise(sum_int = sum(Intensity, na.rm = TRUE), .groups = "drop") %>%
    mutate(Condition = label)
}

ctrl_long <- reshape_plate(control,  "Control")
low_long  <- reshape_plate(lowinput, "LowInput")

# Remove checks
ctrl_long <- ctrl_long %>%
  filter(!str_detect(Sample, "CHECK"))
low_long <- low_long %>%
  filter(!str_detect(Sample, "Check"))

# Transfrom to wide format
ctrl_wide <- ctrl_long %>%
  pivot_wider(names_from = Class, values_from = sum_int, values_fill = 0) %>%
  dplyr::select(-2)

low_wide <- low_long %>%
  pivot_wider(names_from = Class, values_from = sum_int, values_fill = 0) %>%
  dplyr::select(-2)

# Add "Sum_" to every string in colnames
colnames(ctrl_wide)[-1] <- paste0("Sum_", colnames(ctrl_wide)[-1])
colnames(low_wide)[-1] <- paste0("Sum_", colnames(low_wide)[-1])




# 1) grab just the lipid columns
lipid_cols_ctrl <- names(ctrl_wide)[-1] 
lipid_cols_low  <- names(low_wide)[-1]

# 2) all combinations of 2
pairs_ctrl <- combn(lipid_cols_ctrl, 2, simplify = FALSE)
pairs_low  <- combn(lipid_cols_low, 2, simplify = FALSE)

# 3) compute each ratio and name it "Sum_X/Sum_Y"
ratio_df_ctrl <- map_dfc(pairs_ctrl, ~ {
  num <- .x[1]; den <- .x[2]
  ctrl_wide[[num]] / ctrl_wide[[den]]
}) %>% 
  set_names(map_chr(pairs_ctrl, ~ paste(.x, collapse = "/")))

ratio_df_low <- map_dfc(pairs_low, ~ {
  num <- .x[1]; den <- .x[2]
  low_wide[[num]] / low_wide[[den]]
}) %>% 
  set_names(map_chr(pairs_low, ~ paste(.x, collapse = "/")))


# 4) bind the Sample column back on front
ctrl_wide_ratios <- bind_cols(
  ctrl_wide %>% dplyr::select(Sample),
  ratio_df
)
low_wide_ratios <- bind_cols(
  low_wide %>% dplyr::select(Sample),
  ratio_df_low
)

# result: ctrl_wide_ratios has 394 rows × (1 + choose(17,2)=137) cols
ctrl_wide_ratios
# result: low_wide_ratios has 363 rows × (1 + choose(17,2)=137) cols
low_wide_ratios

# Combine the sums as well
ctrl_wide_ratios <- ctrl_wide_ratios %>%
  left_join(ctrl_wide, by = "Sample")
low_wide_ratios <- low_wide_ratios %>%
  left_join(low_wide, by = "Sample")


# Save the results
write.csv(ctrl_wide_ratios, "data/SPATS_fitted/BLUP_GWAS_phenotype/control_all_lipids_BLUPs_sum_ratios.csv", row.names = FALSE, quote = FALSE)
write.csv(low_wide_ratios, "data/SPATS_fitted/BLUP_GWAS_phenotype/lowinput_all_lipids_BLUPs_sum_ratios.csv", row.names = FALSE, quote = FALSE)

