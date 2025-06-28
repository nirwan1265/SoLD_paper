#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
#### SOrghum Lipidomics Database (SOLD)
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


###############################################################################
## Supplementary Figures – Overlap of Control vs Low‑input lipid lists
## Three panels:
##   S1 – All detected lipids               (Venn)
##   S2 – Traditional classes only          (Venn)
##   S3 – Curated “Class” annotation set    (table + bar *optional*)
###############################################################################

# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------
library(vroom)      # fast CSV import
library(dplyr)      # data wrangling
library(tidyr)      # pivoting
library(ggplot2)    # plotting
library(stringr)    # string manipulation
library(forcats)    # factor reordering
library(tidyverse)
library(ggpubr)      # For stat_compare_means
library(scales) 
library(ggVennDiagram)   # nice, ggplot‑based 2‑set Venn
library(readr) 
library(sf)
library(ggVennDiagram) # For Venn diagrams)
library(viridis)

# ------------------------------------------------------------------------------
# 2. Read in the data
# ------------------------------------------------------------------------------
control  <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Control_all_lipids_final_non_normalized.csv")

lowinput <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv")
colnames(control)
quartz()
hist(control$`.alpha.-Ionone`)

# Remove PC(17:0) from control and lowinput (internal standard) from the columns
control  <- control %>% dplyr::select(-`PC(17:0)`)
lowinput <- lowinput %>% dplyr::select(-`PC(17:0)`)

#### CHANGE THE NAME TO COMMONNAME
# Lipid class
lipid_class_info <- vroom::vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/data/lipid_class.csv", show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Class))


# Create a named vector of replacements: names = original, values = new names
name_map <- lipid_class_info %>%
  filter(!is.na(CommonName)) %>%
  select(Lipids, CommonName) %>%
  deframe()  # turns it into a named vector: Lipids -> CommonName

# Function to rename lipid columns in any dataset
rename_lipid_columns <- function(df, name_map) {
  original_names <- colnames(df)
  
  renamed <- original_names %>%
    # Replace only if a match exists
    sapply(function(x) if (x %in% names(name_map)) name_map[[x]] else x)
  
  # Apply renamed vector back
  colnames(df) <- renamed
  return(df)
}

# Apply to both control and lowinput
control  <- rename_lipid_columns(control,  name_map)
lowinput <- rename_lipid_columns(lowinput, name_map)



# Lipid class
lipid_class_info <- vroom::vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/data/lipid_class.csv", show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Class))
colnames(control)
# Lipid class Traditional lipids
traditional_lipid_classes <- c("Glycerolipid", "Glycerophospholipid", "Glycoglycerolipid",
                               "Sphingolipid", "Sterol", "Betaine lipid", "Fatty acid and derivative")
lipid_class_traditional <- lipid_class_info %>% 
  dplyr::filter(Class %in% traditional_lipid_classes)

# Lipid class Traditional lipids
nontraditional_lipid_class <- c("N-acylethanolamine","Terpenoid","Prenol","Tetrapyrrole","Vitamin")

lipid_class_nontraditional <- lipid_class_info %>% 
  dplyr::filter(Class %in% nontraditional_lipid_class)


# Use this updated lipid_class_info
lipid_class_info <- lipid_class_info %>%
  mutate(final_name = ifelse(is.na(CommonName), Lipids, CommonName))

# Redefine traditional and nontraditional sets
lipid_class_traditional <- lipid_class_info %>%
  filter(Class %in% traditional_lipid_classes)

lipid_class_nontraditional <- lipid_class_info %>%
  filter(Class %in% nontraditional_lipid_class)

# Use final_name column now
traditional_lipids     <- lipid_class_traditional$final_name
nontraditional_lipids  <- lipid_class_nontraditional$final_name

# Assuming you renamed control/lowinput already
trad_cols_control     <- intersect(traditional_lipids, colnames(control))
nontrad_cols_control  <- intersect(nontraditional_lipids, colnames(control))

trad_cols_lowinput    <- intersect(traditional_lipids, colnames(lowinput))
nontrad_cols_lowinput <- intersect(nontraditional_lipids, colnames(lowinput))

# Subset
control_traditional     <- control  %>% select(Compound_Name, all_of(trad_cols_control))
control_nontraditional  <- control  %>% select(Compound_Name, all_of(nontrad_cols_control))
lowinput_traditional    <- lowinput %>% select(Compound_Name, all_of(trad_cols_lowinput))
lowinput_nontraditional <- lowinput %>% select(Compound_Name, all_of(nontrad_cols_lowinput))




# ────────────────────────────────────────────────────────────────
# 1.  Helper: merge two matrices that share `Compound_Name`
#     and then log-transform all numeric columns
# ────────────────────────────────────────────────────────────────
make_log10_matrix <- function(df_trad, df_non, outfile, pseudo = NULL) {
  
  # a) merge (keeps every lipid column, no duplication of Compound_Name)
  merged <- full_join(df_trad, df_non, by = "Compound_Name")
  
  # b) pick a pseudocount: half of the smallest positive value if not supplied
  if (is.null(pseudo)) {
    min_pos <- merged %>%
      select(-Compound_Name) %>%
      unlist(use.names = FALSE) %>%
      .[. > 0] %>%
      min(na.rm = TRUE)
    
    pseudo <- min_pos / 2
  }
  
  # c) log10-transform every numeric column
  log_mat <- merged %>%
    mutate(across(-Compound_Name, ~ log10(.x + pseudo)))
  
  # d) save
  write.csv(log_mat, outfile, row.names = FALSE)
  message("✓ wrote ", outfile, "  ( pseudocount = ", signif(pseudo, 3), " )")
  
  invisible(log_mat)
}

# ────────────────────────────────────────────────────────────────
# 2.  Apply to your four data frames
#     (these were created earlier in your script)
# ────────────────────────────────────────────────────────────────
log_control  <- make_log10_matrix(
  control_traditional,
  control_nontraditional,
  "control_allLipids_log10.csv"        # <- output path
)

log_lowinput <- make_log10_matrix(
  lowinput_traditional,
  lowinput_nontraditional,
  "lowinput_allLipids_log10.csv"
)


# Save the data
write.csv(log_control, "Control_allLipids_log10_GWAS.csv", row.names = FALSE)
write.csv(log_lowinput, "Lowinput_allLipids_log10_GWAS.csv", row.names = FALSE)


