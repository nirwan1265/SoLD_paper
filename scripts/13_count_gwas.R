# Load necessary library
library(dplyr)
library(stringr)
library(vroom)
library(purrr)

# List all .txt files in the directory

# Control
# dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/control/annotation/all_lipids/"
# file_list <- list.files(path = dir,pattern = "*.txt")

# Lowinput
# dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/lowinput/annotation/New Folder With Items/"
# dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/lowinput/annotation/"
# 
# # All
dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/lowinput/annotation/all_lipids/"

file_list <- list.files(path = dir,pattern = "*.txt")


-log10(0.05/234264)


# ## MG and DG
# file_list <- file_list[c(9:24,38:42)]
# 
# # PE and DG
# file_list <- file_list[c(90:105)]
# 
# # SQDG
# file_list <- file_list[c(112:117)]
# 
# # Neutral lipids
# file_list <- file_list[c(38:42,9:24)]
# 
# 
# # Galactolipids
# file_list <- file_list[c(43:47,25:31)]
# 
# # Phospholipids
file_list <- file_list[c(74:75,77:135,137)]

# Select only that start with TG in file_list vector
#file_list <- file_list[str_detect(file_list, "^PG")]


# Control
# dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/control/annotation/remaining_lipids/"
# file_list <- list.files(path = dir,pattern = "*.txt")
# file_list <- file_list[5]



### Sums and Ratios

# Control
dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/sum_ratio/lowinput/annotation/"
file_list <- list.files(path = dir,pattern = "*.txt")
file_list <- file_list[c(164,160,155,149,65,
                         118,132,57,131,130,58,129,
                         147,153,158,140,
                         42,51,143,47,46,
                         96,102,95,97,113,100,101,
                         170,165,135,
                         125,48,115,63,112,120,119, 121,41,109,50,105)]

file_list <- file_list[c(164,160,155,149,65)]

#Lowinput
# dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/sum_ratio/lowinput/annotation/"
# file_list <- list.files(path = dir,pattern = "*.txt")
# 
# # PE DG and TG LPE
# file_list <- file_list[c(46,116)]
# 
# 
# # Sum MG
# file_list <- file_list[str_detect(file_list, "Sum_MG")]
# file_list <- file_list[17]
# 
# 
# 
# # Sum SQDG
# file_list <- file_list[str_detect(file_list, "Sum_PG")]
# file_list <- file_list[14]
# 
# 
# # Sum PG
# file_list <- file_list[str_detect(file_list, "Sum_SQDG")]
# file_list <- file_list[17]
# 
# # Sum and Ratio Galactolipids
# # Select all files with MG or MGDG in the name
# pattern <- "(?<=^Sum_)(MGDG|DGDG)(?=_)|(?<=_Sum_)(MGDG|DGDG)(?=_)"
# file_list <- file_list[str_detect(file_list, pattern)]
# # Remove all files with FA or Cer or GalCer, SM, or AEG
# file_list <- file_list[!str_detect(file_list, "FA|Cer|GalCer|SM|AEG")]
# 
# 
# # Sum and Ratio Sulfolipids
# # Select all files with MG or MGDG in the name
# pattern <- "(?<=^Sum_)(SQDG)(?=_)|(?<=_Sum_)(SQDG)"
# file_list <- file_list[str_detect(file_list, pattern)]
# # Remove all files with FA or Cer or GalCer, SM, or AEG
# file_list <- file_list[!str_detect(file_list, "FA|Cer|GalCer|SM|AEG")]
# 
# # Sum and Ratio Phospholipids
# # Select all files with MG or MGDG in the name
# pattern <- "(?<=^Sum_)(PC|PE|PG|PA)(?=_)|(?<=_Sum_)(PC|PE|PG|PA)"
# file_list <- file_list[str_detect(file_list, pattern)]
# # Remove all files with FA or Cer or GalCer, SM, or AEG
# file_list <- file_list[!str_detect(file_list, "FA|Cer|GalCer|SM|AEG")]
# 
# 
# # Sum and Ratio Neutral
# # Select all files with MG or MGDG in the name
# pattern <- "(?<=^Sum_)(DG|MG|EG)(?=_)|(?<=_Sum_)(MG|DG|TG)"
# file_list <- file_list[str_detect(file_list, pattern)]
# # Remove all files with FA or Cer or GalCer, SM, or AEG
# file_list <- file_list[!str_detect(file_list, "FA|Cer|GalCer|SM|AEG")]
# 

# Empty vectors
results <- tibble()          
p_value_threshold <- 1e-7
-log10(p_value_threshold)

# Loop through the files
if (!exists("results")) results <- dplyr::tibble()

for (fp in paste0(dir, file_list)) {
  
  ## ---------- read file ---------------------------------------------------
  dat <- vroom(fp, show_col_types = FALSE)
  
  ## ---------- phenotype from file name ------------------------------------
  phen <- sub("\\.txt$", "", basename(fp))
  
  ## ---------- keep only significant rows & compute best SNP ---------------
  sig <- dat %>%
    dplyr::filter(MinPValue < p_value_threshold) %>%
    dplyr::mutate(
      Phenotype = phen,
      snp_vec  = stringr::str_split(SNPs, ","),                           # list<chr>
      pval_vec = stringr::str_split(PValues, ",") |>
        purrr::map(~ as.numeric(trimws(.))),                    # list<num>
      n_snps   = lengths(snp_vec),                                        # count SNPs
      idx_min  = purrr::map_int(pval_vec, ~ if (length(.x)) which.min(.x) else NA_integer_),
      best_snp = purrr::map2_chr(snp_vec, idx_min, ~
                                   ifelse(is.na(.y) || .y < 1 || .y > length(.x),
                                          NA_character_, stringr::str_trim(.x[.y]))),
      best_pval = purrr::map2_dbl(pval_vec, idx_min, ~
                                    ifelse(is.na(.y) || .y < 1 || .y > length(.x),
                                           NA_real_, .x[.y]))
    ) %>%
    dplyr::select(Chromosome, GeneID, SNPs, MinPValue,
                  Phenotype, n_snps, best_snp, best_pval)
  
  results <- dplyr::bind_rows(results, sig)
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
write.table(combined_results, "table/GWAS_results/Phospholipids_common_genes_individual_logpvalue_7_spat_fitted_BLUP.txt", row.names=F, quote=F, sep="\t")
getwd()


### Make a table for latex
combined_results <- combined_results[1:16,-c(4:7)]


# Escape underscores in the whole dataframe
combined_results_fixed <- combined_results %>%
  dplyr::mutate(across(everything(), ~ gsub("_", "\\\\_", .x)))

combined_results_fixed <- combine_results
# Now use kableExtra/xtable/gt to export safely
library(knitr)
library(kableExtra)

kable(combined_results_fixed, format = "latex", booktabs = TRUE,
      caption = "GWAS hits with associated lipid phenotypes") %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))


# write.table(combined_results, "table/GWAS_results/11_14_17_Eicosatrienoic_acid__Z_Z_Z___common_genes_lowinput_individual_logpvalue_5_spat_fitted_BLUP.txt", row.names=F, quote=F, sep="\t")
# write.table(combined_results, "table/GWAS_results/gibberellic_acid_common_genes_lowinput_individual_logpvalue_4_spat_fitted_BLUP.txt", row.names=F, quote=F, sep="\t")









###### MAKE FULL ANNOTATION

# Lowinput
# Gene pvalue file
ind <- vroom("/Users/nirwantandukar/Documents/Github/SoLD_paper/table/GWAS_results/All_lipids_Lowinput_common_genes_individual_logpvalue_7_spat_fitted_BLUP.txt")
sum_ratio <- vroom("/Users/nirwantandukar/Documents/Github/SoLD_paper/table/GWAS_results/All_lipids_Lowinput_common_genes_sum_ratio_logpvalue_7_spat_fitted_BLUP.txt")


# Control
ind <- vroom("/Users/nirwantandukar/Documents/Github/SoLD_paper/table/GWAS_results/All_lipids_Control_common_genes_individual_logpvalue_7_spat_fitted_BLUP.txt")
sum_ratio <- vroom("/Users/nirwantandukar/Documents/Github/SoLD_paper/table/GWAS_results/All_lipids_Control_common_genes_sum_ratio_logpvalue_7_spat_fitted_BLUP.txt")


# annotation file
# Lowinput
ind_annot <- vroom("/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/lowinput/gene_annotation/annotation_All_lipids_Lowinput_common_genes_individual_logpvalue_7_spat_fitted_BLUP.txt")
sum_ratio_annot <- vroom("/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/lowinput/gene_annotation/annotation_All_lipids_Lowinput_common_genes_sum_ratio_logpvalue_7_spat_fitted_BLUP.txt")
colnames(sum_ratio_annot)[1] <- "GeneID"
colnames(sum_ratio_annot)[3] <- "GeneName"


# Control
ind_annot <- vroom("/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/control/gene_annotation/annotation_All_lipids_Control_common_genes_individual_logpvalue_7_spat_fitted_BLUP.txt")
sum_ratio_annot <- vroom("/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/sum_ratio/control/gene_annotation/annotation_All_lipids_Control_common_genes_sum_ratio_logpvalue_7_spat_fitted_BLUP.txt")

# Combine the two files for ind and sum ratio separately with gene and GeneID
ind_full_keep <- ind %>%
  full_join(ind_annot, by = join_by(gene == GeneID), keep = TRUE) %>%
  dplyr::select(-GeneID)

sum_ratio_full_keep <- sum_ratio %>%
  full_join(sum_ratio_annot, by = join_by(gene == GeneID), keep = TRUE) %>%
  dplyr::select(-GeneID)



# Save the full annotation
write.table(ind_full_keep, "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/lowinput/full_annotation/All_lipids_Control_common_genes_individual_logpvalue_7_spat_fitted_BLUP_full_annotation.txt", row.names=F, quote=F, sep ="\t")

write.table(sum_ratio_full_keep, "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/final_results/lowinput/full_annotation/All_lipids_Control_common_genes_sum_ratio_logpvalue_7_spat_fitted_BLUP_full_annotation.txt", row.names=F, quote=F, sep = "\t")
