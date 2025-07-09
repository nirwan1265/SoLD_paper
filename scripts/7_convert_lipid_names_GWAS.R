library(dplyr)
library(stringr)
library(tibble)



control <- vroom("/Users/nirwantandukar/Documents/Github/SoLD_paper/data/SPATS_fitted/BLUP_GWAS_phenotype/control_all_lipids_BLUPs.csv")

# 1. Store original column names
orig_names <- colnames(control)

# 2. Convert non-alphanumeric characters to underscores
converted_names <- str_replace_all(orig_names, "[^A-Za-z0-9]", "_")

# 3. Apply to your data frame
colnames(control) <- converted_names

# 4. Make mapping data frame
name_mapping <- tibble(
  RealName = orig_names,
  ConvertedName = converted_names
)
name_mapping

colnames(control) <- name_mapping$ConvertedName


# save 
write.csv(control, "/Users/nirwantandukar/Documents/Github/SoLD_paper/data/SPATS_fitted/BLUP_GWAS_phenotype/control_all_lipids_BLUPs_rename.csv")
write.csv(name_mapping, 
          "/Users/nirwantandukar/Documents/Github/SoLD_paper/data/SPATS_fitted/BLUP_GWAS_phenotype/control_lipid_GWAS_map.csv",
          row.names = FALSE)

lowinput <- vroom("/Users/nirwantandukar/Documents/Github/SoLD_paper/data/SPATS_fitted/BLUP_GWAS_phenotype/lowinput_all_lipids_BLUPs.csv")
# 1. Store original column names
orig_names_lowinput <- colnames(lowinput)
# 2. Convert non-alphanumeric characters to underscores
converted_names_lowinput <- str_replace_all(orig_names_lowinput, "[^A-Za-z0-9]", "_")
# 3. Apply to your data frame
colnames(lowinput) <- converted_names_lowinput
# 4. Make mapping data frame
name_mapping_lowinput <- tibble(
  RealName = orig_names_lowinput,
  ConvertedName = converted_names_lowinput
)
name_mapping_lowinput
colnames(lowinput) <- name_mapping_lowinput$ConvertedName

# save
write_csv(lowinput, "/Users/nirwantandukar/Documents/Github/SoLD_paper/data/SPATS_fitted/BLUP_GWAS_phenotype/lowinput_all_lipids_BLUPs_rename.csv")
write.csv(name_mapping_lowinput, 
          "/Users/nirwantandukar/Documents/Github/SoLD_paper/data/SPATS_fitted/BLUP_GWAS_phenotype/lowinput_lipid_GWAS_map.csv",
          row.names = FALSE)
