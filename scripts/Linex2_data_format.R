# Read the raw files 
control <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  select(-c(2,3,4)) %>%
  rename(Compound_Name = 1) %>%
  mutate(Compound_Name = paste0(Compound_Name, "_control"))

lowinput <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  select(-c(2,3,4)) %>%
  rename(Compound_Name = 1) %>%
  mutate(Compound_Name = paste0(Compound_Name, "_lowinput"))


# Valid classes
valid_classes <- c("DGDG","MGDG",
                   "TG","DG","MG",
                   "PC","PE",
                   "SQDG",
                   "LPC","LPE","PG","PA","PS"
)

class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

control_sub <- control %>%
  select(
    Compound_Name,
    matches(class_pat)
  )

lowinput_sub <- lowinput %>%
  select(
    Compound_Name,
    matches(class_pat)
  )

# inspect:
colnames(control_sub)
colnames(lowinput_sub)

# Get the common set of compounds
common_compounds <- intersect(colnames(control_sub), colnames(lowinput_sub))
# Remove Compound_Name
common_compounds <- common_compounds[common_compounds != "Compound_Name"]


 
# 1) pick only those columns (plus your sample ID & condition)
ctrl_trim <- control_sub %>%
  select(Compound_Name, all_of(common_compounds))

low_trim <- lowinput_sub %>%
  select(Compound_Name, all_of(common_compounds))

# 2) stack them into one df
shared_wide <- bind_rows(ctrl_trim, low_trim)

# check
glimpse(shared_wide)


# Save 
write.csv(shared_wide,"Linex2_data.csv",row.names = F)


library(tidyr)
formatted <- shared_wide %>% 
  # 1) break Compound_Name into sample ID + Condition
  separate(Compound_Name, into = c("Sample","Condition"), sep = "_", extra = "merge") %>% 
  
  # 2) within each Condition, number the replicates in order
  group_by(Condition) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  
  # 3) pivot to long so we have one row per Sample×Lipid
  pivot_longer(-c(Sample,Condition,replicate),
               names_to  = "Lipid",
               values_to = "Intensity") %>%
  
  # 4) pivot back to wide, naming columns like control_1, lowinput_1, control_2…
  unite(col = "cond_rep", Condition, replicate, sep = "_") %>%
  pivot_wider(names_from  = cond_rep,
              values_from = Intensity) %>%
  
  # 5) move Lipid to rownames if you really want that, or keep as a column
  column_to_rownames("Lipid")

formatted



