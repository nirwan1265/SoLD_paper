# Load necessary library
library(dplyr)
library(stringr)
library(vroom)
library(purrr)

# List all .txt files in the directory

# Lowinput
# dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/lowinput/annotation/New Folder With Items/"
# #dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/lowinput/annotation/"
# file_list <- list.files(path = dir,pattern = "*.txt")
# file_list <- file_list[29]

# Select only that start with TG in file_list vector
#file_list <- file_list[str_detect(file_list, "^PG")]


# Control
# dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/control/annotation/remaining_lipids/"
# file_list <- list.files(path = dir,pattern = "*.txt")
# file_list <- file_list[5]



### Sums and Ratios
# Lowinput
dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/sum_ratio/lowinput/annotation/"
file_list <- list.files(path = dir,pattern = "*.txt")
file_list <- file_list[-c(1:70,133:186,331:336)]

file_list <- file_list[c(36,52,94,106,117,127,136,144,151,157,162,169,171)]



# Empty vectors
results <- tibble()          
p_value_threshold <- 1e-7
-log10(p_value_threshold)

# Loop through the files
for (fp in paste0(dir, file_list)) {
  
  ## ---------- read file ---------------------------------------------------
  dat <- vroom(fp, show_col_types = FALSE)
  
  ## ---------- phenotype from file name ------------------------------------
  phen <- sub("\\.txt$", "", basename(fp))
  
  ## ---------- keep only significant rows ----------------------------------
  sig <- dat %>%
    dplyr::filter(MinPValue < p_value_threshold)                         %>%
    dplyr::select(Chromosome, GeneID, SNPs, MinPValue)                   %>%
    dplyr::mutate(
      Phenotype = phen,
      n_snps    = str_count(SNPs, ",") + 1,                       # count SNPs
      best_snp  = str_split(SNPs, ",") |>                         # first SNP
        map_chr(1),
      best_pval = MinPValue                                       # already lowest
    )
  
  results <- bind_rows(results, sig)
}

## -------------------------------------------------------------------------
##  FINAL COLLAPSE: one line per Chromosome × gene
## -------------------------------------------------------------------------
combined_results <- results %>%
  group_by(Chromosome, gene = GeneID) %>%
  summarise(
    phenotype        = paste(Phenotype,  collapse = ";"),         # order = loop order
    n_snps           = paste(n_snps,     collapse = ";"),
    best_snp         = paste(best_snp,   collapse = ";"),
    best_pvalue      = paste(best_pval,  collapse = ";"),
    phenotype_number = n_distinct(Phenotype),
    .groups          = "drop"
  )

head(combined_results)


# directory prefix you want to strip (escape trailing slash if present)
dir_prefix <- paste0("^", dir)    # dir is your path string

combined_results <- combined_results %>%
  mutate(
    phenotype = map_chr(phenotype, function(x) {
      # split at commas → clean each piece → paste back
      str_split(x, ";")[[1]]                   |>
        sub(dir_prefix, "", x = _)             |>  # remove dir prefix
        str_remove("_mod_sub_.*$")             |>  # trim suffix
        unique()                               |>  # de‑duplicate if needed
        paste(collapse = ";")
    }),
    phenotype_number = str_count(phenotype, ";") + 1
  )

head(combined_results)

# arrange by phenotype number
combined_results <- combined_results %>%
  dplyr::arrange(desc(phenotype_number))

head(combined_results)

# Save the file
write.table(combined_results, "table/GWAS_results/All_sum_ratio_Lowinput_common_genes_individual_logpvalue_7_spat_fitted_BLUP.txt", row.names=F, quote=F, sep="\t")
getwd()
# write.table(combined_results, "table/GWAS_results/11_14_17_Eicosatrienoic_acid__Z_Z_Z___common_genes_lowinput_individual_logpvalue_5_spat_fitted_BLUP.txt", row.names=F, quote=F, sep="\t")
# write.table(combined_results, "table/GWAS_results/gibberellic_acid_common_genes_lowinput_individual_logpvalue_4_spat_fitted_BLUP.txt", row.names=F, quote=F, sep="\t")

