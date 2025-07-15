# Load necessary library
library(dplyr)
library(stringr)
library(vroom)



# List all .txt files in the directory

# All
# dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/annotation/lowinput/all/"
# file_list <- list.files(path = dir,pattern = "*.txt")

# Traditional lipids
dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/annotation/lowinput/Traditional_lipids/"
file_list <- list.files(path = dir,pattern = "*.txt")

# Non-Traditional lipids
dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/annotation/lowinput/Nontraditional_lipids/"
file_list <- list.files(path = dir,pattern = "*.txt")
file_list <- file_list[60]

# Sums and Ratios
dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/sum_ratio_BLUP/annotation/"
file_list <- list.files(path = dir,pattern = "*.txt")




# Empty vectors
results <- tibble()          
p_value_threshold <- 1e-5
-log10(p_value_threshold)

# Loop through the files
for (fp in paste0(dir, file_list)) {
  
  ## ---------- read file ---------------------------------------------------
  dat <- vroom(fp, show_col_types = FALSE)
  
  ## ---------- phenotype from file name ------------------------------------
  phen <- sub("\\.txt$", "", basename(fp))
  
  ## ---------- keep only significant rows ----------------------------------
  sig <- dat %>%
    filter(MinPValue < p_value_threshold)                         %>%
    select(Chromosome, GeneID, SNPs, MinPValue)                   %>%
    mutate(
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
    phenotype        = paste(Phenotype,  collapse = ","),         # order = loop order
    n_snps           = paste(n_snps,     collapse = ","),
    best_snp         = paste(best_snp,   collapse = ","),
    best_pvalue      = paste(best_pval,  collapse = ","),
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
      str_split(x, ",")[[1]]                   |>
        sub(dir_prefix, "", x = _)             |>  # remove dir prefix
        str_remove("_mod_sub_.*$")             |>  # trim suffix
        unique()                               |>  # de‑duplicate if needed
        paste(collapse = ",")
    }),
    phenotype_number = str_count(phenotype, ",") + 1
  )

head(combined_results)

# arrange by phenotype number
combined_results <- combined_results %>%
  dplyr::arrange(desc(phenotype_number))

head(combined_results)

# Save the file
# write.table(combined_results, "table/GWAS_results/NonTraditional_lipids_common_genes_lowinput_individual_logpvalue_5_spat_fitted_BLUP.txt", row.names=F, quote=F, sep="\t")
getwd()
# write.table(combined_results, "table/GWAS_results/SQDG_32_0_common_genes_lowinput_individual_logpvalue_5_spat_fitted_BLUP.txt", row.names=F, quote=F, sep="\t")
write.table(combined_results, "table/GWAS_results/gibberellic_acid_common_genes_lowinput_individual_logpvalue_4_spat_fitted_BLUP.txt", row.names=F, quote=F, sep="\t")

