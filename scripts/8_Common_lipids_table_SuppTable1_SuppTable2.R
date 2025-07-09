#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
#### SOrghum Lipidomics Database (SOLD)
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

###############################################################################
## Supplementary Tables – Overlap of Control vs Low‑input lipid lists
## Three panels:
##   S1 – All detected lipids               (Venn)
##   S2 – Traditional classes only          (Venn)
##   S3 – Curated “Class” annotation set    (table + bar *optional*)
###############################################################################



# ──────────────────────────────────────────────────────────────
#  0)  PACKAGES
# ──────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(stringr)
library(vroom)

# ──────────────────────────────────────────────────────────────
#  1)  READ & COMBINE RAW INTENSITIES
# ──────────────────────────────────────────────────────────────

# Control
control <- vroom("data/SPATS_fitted/non_normalized_intensities/control_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  select(-c(2,3,4)) %>%
  rename(Compound_Name = 1) %>%
  mutate(Condition = "Control") #%>%
  #select(-`PC(17:0)`)

# Lowinput
lowinput <- vroom("data/SPATS_fitted/non_normalized_intensities/lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  select(-c(2,3,4)) %>%
  rename(Compound_Name = 1) %>%
  mutate(Condition = "LowInput") #%>%
  #select(-`PC(17:0)`)

raw_all <- bind_rows(control, lowinput)

# ──────────────────────────────────────────────────────────────
#  2)  TIDY LONG
# ──────────────────────────────────────────────────────────────

valid_classes <- c("TG","DG","MG",
                   "PC","PE","PI",
                   "DGDG","MGDG",
                   "SQDG","SM","AEG",
                   "LPC","LPE","PG","PA")
class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

df_long <- raw_all %>%
  pivot_longer(
    cols      = -c(Compound_Name, Condition),
    names_to  = "Lipid",
    values_to = "Intensity"
  ) %>%
  mutate(
    Class  = str_extract(Lipid, class_pat),
    Sample = Compound_Name
  ) %>%
  filter(!is.na(Class))

# ──────────────────────────────────────────────────────────────
#  3)  COMPUTE WITHIN-CLASS PERCENTAGES (no lumping)
# ──────────────────────────────────────────────────────────────

compute_pct <- function(df, cond) {
  df %>%
    filter(Condition == cond) %>%
    group_by(Sample, Class) %>%
    mutate(class_total = sum(Intensity, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(rel_within = if_else(class_total > 0, Intensity / class_total, 0)) %>%
    group_by(Class, Lipid) %>%
    summarise(mean_prop = mean(rel_within, na.rm = TRUE), .groups = "drop") %>%
    mutate(mean_pct = round(mean_prop * 100, 1)) %>%   # one decimal place
    arrange(Class, desc(mean_pct))
}

table_control  <- compute_pct(df_long, "Control")
table_lowinput <- compute_pct(df_long, "LowInput")

# ──────────────────────────────────────────────────────────────
#  4)  MERGE & FORMAT FOR COMPARISON
# ──────────────────────────────────────────────────────────────

ctrl <- table_control  %>% rename(Pct_Control    = mean_pct, Species = Lipid)
lwp  <- table_lowinput %>% rename(Pct_LowInput  = mean_pct, Species = Lipid)

combined <- full_join(ctrl, lwp, by = c("Class", "Species")) %>%
  arrange(Class, desc(ifelse(is.na(Pct_Control), 0, Pct_Control) +
                        ifelse(is.na(Pct_LowInput),  0, Pct_LowInput)),
          Species) %>%
  dplyr::select(-mean_prop.x) %>%
  dplyr::select(-mean_prop.y) 
    

# clean column names
colnames(combined)[3:4] <- c("Control_Percentage", "LowInput_Percentage")

# drop Class column if you just want Species & the two percentages
# combined <- combined[ , -1 ]

print(combined, n = Inf)

# Save:
write.csv(combined, "SuppTable1_overlap_control_lowinput_lipid_species_percentages.csv", row.names = FALSE)


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
# control  <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Control_all_lipids_final_non_normalized.csv")
# 
# lowinput <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv")
# colnames(control)

control  <- vroom("data/SPATS_fitted/non_normalized_intensities/control_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4)) 
colnames(control)[1] <- "Compound_Name"  

lowinput  <- vroom("data/SPATS_fitted/non_normalized_intensities/lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
colnames(lowinput)[1] <- "Compound_Name"  



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



# 1) First, make sure you have a lookup from final_name → Class:
lookup <- lipid_class_info %>%
  mutate(final_name = ifelse(is.na(CommonName), Lipids, CommonName)) %>%
  select(final_name, Class)

# 2) Grab the lipid‐column names (minus the sample ID) for each condition/set
ctr_trad_cols   <- setdiff(names(control_traditional),    "Compound_Name")
low_trad_cols   <- setdiff(names(lowinput_traditional),   "Compound_Name")
ctr_non_cols    <- setdiff(names(control_nontraditional), "Compound_Name")
low_non_cols    <- setdiff(names(lowinput_nontraditional),"Compound_Name")

# 3) A helper that, for a vector of classes, builds the summary table
make_summary <- function(classes, ctrl_cols, low_cols){
  lapply(classes, function(cl){
    # which lipids in this class?
    members <- lookup %>% filter(Class == cl) %>% pull(final_name)
    in_ctrl <- intersect(ctrl_cols, members)
    in_low  <- intersect(low_cols,  members)
    tibble(
      Class            = cl,
      Shared           = paste0(intersect(in_ctrl, in_low), collapse = "; "),
      Unique_Control   = paste0(setdiff(in_ctrl, in_low), collapse = "; "),
      Unique_Lowinput  = paste0(setdiff(in_low,  in_ctrl), collapse = "; ")
    )
  }) %>% bind_rows()
}

# 4) Build both tables
trad_summary <- make_summary(
  traditional_lipid_classes,
  ctr_trad_cols, low_trad_cols
)

nontrad_summary <- make_summary(
  nontraditional_lipid_class,
  ctr_non_cols, low_non_cols
)

# 5) Inspect
print(trad_summary)
print(nontrad_summary)

# 6) (Optional) write out
write.csv(trad_summary,    "SuppTable_2_traditional_lipid_class_overlap.csv",    row.names = FALSE)
write.csv(nontrad_summary, "SuppTable_3_nontraditional_lipid_class_overlap.csv", row.names = FALSE)
