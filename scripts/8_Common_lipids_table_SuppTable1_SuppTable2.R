#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
#### SOrghum Lipidomics Database (SOLD)
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

###############################################################################
## Supplementary Tables – Overlap of Control vs Low‑input lipid lists
## Three panels:
##   S1 – All detected lipids               
##   S2 – Traditional classes only          
##   S3 – Curated “Class” annotation set    
###############################################################################



###############################################################################
### Suuplmentary Table 1
###############################################################################

### Load the packages
library(vroom)      
library(dplyr)      
library(tidyr)      
library(ggplot2)    
library(stringr)    
library(forcats)    
library(tidyverse)
library(ggpubr)     
library(scales) 
library(ggVennDiagram)  
library(readr) 
library(sf)
library(ggVennDiagram)
library(viridis)


### Read and combine the raw intensities for control and lowinput
# Control
control <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  select(-c(2,3,4)) %>%
  rename(Compound_Name = 1) %>%
  mutate(Condition = "Control")

# Lowinput
lowinput <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  select(-c(2,3,4)) %>%
  rename(Compound_Name = 1) %>%
  mutate(Condition = "LowInput")

# Combine
raw_all <- bind_rows(control, lowinput)


### Sub the different lipid species
valid_classes <- c("TG","DG","MG",
                   "PC","PE","PI",
                   "DGDG","MGDG",
                   "SQDG","SM","AEG",
                   "LPC","LPE","PG","PA","Cer","GalCer","FA")
class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

# Long format
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


### Compute within lipid species percentages
compute_pct <- function(df, cond) {
  df %>%
    filter(Condition == cond) %>%
    group_by(Sample, Class) %>%
    mutate(class_total = sum(Intensity, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(rel_within = if_else(class_total > 0, Intensity / class_total, 0)) %>%
    group_by(Class, Lipid) %>%
    summarise(mean_prop = mean(rel_within, na.rm = TRUE), .groups = "drop") %>%
    # Now format percentages, showing "<0.1" if it's truly nonzero but <0.1%
    mutate(
      pct_raw = mean_prop * 100,
      mean_pct = case_when(
        pct_raw == 0            ~ "0.0",           # truly zero
        pct_raw < 0.1           ~ "<0.1",          # nonzero but under 0.1%
        TRUE                    ~ sprintf("%.1f", pct_raw)
      )
    ) %>%
    select(-pct_raw, -mean_prop) %>%
    arrange(Class, desc(as.numeric(sub("<","",mean_pct))))
}


table_control  <- compute_pct(df_long, "Control")
table_lowinput <- compute_pct(df_long, "LowInput")

### Merge for comparison
ctrl <- table_control  %>% rename(Pct_Control    = mean_pct, Species = Lipid)
lwp  <- table_lowinput %>% rename(Pct_LowInput  = mean_pct, Species = Lipid)

combined <- full_join(ctrl, lwp, by = c("Class", "Species")) %>%
  # sort by Class, then by the *sum* of Control+LowInput percentages
  arrange(
    Class,
    desc(
      # for Control: if it starts with "<", treat as 0, else parse number
      (if_else(str_detect(Pct_Control, "^<"), 0, as.numeric(Pct_Control))) +
        # same for LowInput
        (if_else(str_detect(Pct_LowInput,  "^<"), 0, as.numeric(Pct_LowInput)))
    ),
    Species
  ) %>%
  # drop any helper mean_prop columns if they remain
  select(-matches("^mean_prop"))

    

# Clean column names
colnames(combined)[3:4] <- c("Control_Percentage", "LowInput_Percentage")

# SANITY CHECK
print(combined, n = Inf)

# Save:
write.csv(combined, "SuppTable1_overlap_control_lowinput_lipid_species_percentages.csv", row.names = FALSE)


###############################################################################
### Suuplmentary Table 2
###############################################################################


#### Get the lipid class
lipid_class_info <- vroom::vroom("data/lipid_class/final_lipid_classes.csv", show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Class))

# Create a named vector of replacements: names = original, values = new names
name_map <- lipid_class_info %>%
  filter(!is.na(CommonName)) %>%
  dplyr::select(Lipids, CommonName) %>%
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


table(lipid_class_info$Class)

# Lipid class Traditional lipids
lipid_classes <- c("Glycerolipid", "Glycerophospholipid",
                   "Sphingolipid", "Sterol", "Betaine lipid", "Fatty acid","Ether lipid","Terpenoid")



# 1) First, make sure you have a lookup from final_name → Class:
lookup <- lipid_class_info %>%
  mutate(final_name = ifelse(is.na(CommonName), Lipids, CommonName)) %>%
  dplyr::select(final_name, Class)

# 2) Grab the lipid‐column names (minus the sample ID) for each condition/set
ctr   <- setdiff(names(control),    "Compound_Name")
low   <- setdiff(names(lowinput),   "Compound_Name")

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
lipid_summary <- make_summary(
  lipid_classes,
  ctr, low
)


# 5) Inspect
print(lipid_summary)

# 6) (Optional) write out
write.csv(lipid_summary,    "SuppTable_2_lipid_class_overlap.csv",    row.names = FALSE)



