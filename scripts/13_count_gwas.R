# Load necessary library
library(dplyr)
library(stringr)



# List all .txt files in the directory

# Control
# Individual

# All
dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/annotation/lowinput/all/"
file_list <- list.files(path = dir,pattern = "*.txt")

# Traditional lipids
dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/annotation/lowinput/Traditional_lipids/"
file_list <- list.files(path = dir,pattern = "*.txt")

# Non-Traditional lipids
dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/annotation/lowinput/Nontraditional_lipids/"
file_list <- list.files(path = dir,pattern = "*.txt")

# Ratios

dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/sum_ratio_BLUP/annotation/"
file_list <- list.files(path = dir,pattern = "*.txt")

# file_list_sum_ratio <- list.files(path = dir_sum_ratio,pattern = "*.txt")
# file_list_sum_ratio <- file_list_sum_ratio[47]
# file_list <- file_list_sum_ratio
# 
# file_list <- c(file_list,file_list_sum_ratio)
# Select only strings with lowinput in their name
#file_list <- file_list[grep("lowinput", file_list)]

# Select only strings with DG, TG, DGDG, MG, MGDG, SQDG, PC, PE, LPC, LPE in their name
#file_list <- file_list[grep("DG|TG|DGDG|MG|MGDG|SQDG|PC|PE|LPC|LPE", file_list)]

#dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/control/"
#file_list <- list.files(path = "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/control",pattern = "*.txt")

#file_list <- file_list[223]
# Initialize an empty data frame to store the results
results <- data.frame()

# Define the p-value threshold for significance
# Total trait = 215
# boneferroni across traits
# = 0.05/215 = 0.00023
# log10 value is 

p_value_threshold <- 0.0000001
-log10(0.0000001)

#p_value_threshold <- 0.001

# Loop through each file and process the data
for (file in paste0(dir,file_list)) {
  # Read the file
  data <- read.table(file, header = TRUE, sep = " ")
  
  # Extract the phenotype from the file name (removing the .txt extension)
  phenotype <- sub("\\.txt$", "", file)
  
  # Subset the significant genes based on the p-value threshold
  significant_genes <- data %>% 
    dplyr::filter(MinPValue < p_value_threshold) %>% 
    dplyr::select(Chromosome, GeneID, SNPs, MinPValue) %>% 
    dplyr::mutate(Phenotype = phenotype)
  
  # Appexnd the significant genes to the results data frame
  results <- rbind(results, significant_genes)
}

head(results)
# Rename columns to match the desired output
colnames(results) <- c("Chromosome", "gene", "SNP", "pvalue", "phenotype")
head(results)

# Clean up the phenotype column
results$phenotype <- sub(paste0("^",dir), "", results$phenotype)
results$phenotype <- str_remove(results$phenotype, "_mod_sub_.*$")
results$gene <- sub("^gene:", "", results$gene)

# Display the results
head(results)

# Group by Chromosome, gene, and SNP, then summarize pvalue and phenotype
combined_results <- results %>%
  dplyr::group_by(Chromosome, gene) %>%
  dplyr::summarize(
    SNP = paste(unique(SNP), collapse = ","),
    pvalue = paste(pvalue, collapse = ","),
    phenotype = paste(unique(phenotype), collapse = ",")
  ) %>%
  dplyr::ungroup()

head(combined_results)

# Function to count unique phenotypes
count_unique_phenotypes <- function(phenotype_string) {
  phenotypes <- unlist(strsplit(phenotype_string, ","))
  return(length(unique(phenotypes)))
}

# Apply the function to each row to get the phenotype count
combined_results$phenotype_number <- sapply(combined_results$phenotype, count_unique_phenotypes)

# Display the final results with phenotype count
head(combined_results)


# Ensure the column is a character vector
#combined_results$phenotype_number <- as.character(combined_results$phenotype_number)

# Extract the last element after splitting by "/"
#combined_results$phenotype_number <- sapply(strsplit(combined_results$phenotype_number, "/"), tail, 1)

# Print the updated column to verify
#print(combined_results$phenotype_number)

# Print the updated column to verify
#print(combined_results$phenotype_number)

nrow(combined_results)


# sort by phenotype_number
combined_results <- combined_results %>%
  dplyr::arrange(desc(phenotype_number))


# Get 1% ot the top hits
#combined_results <- combined_results[c(1:276),]

# Get the max pvalue for each row
combined_results$max_pvalue <- sapply(strsplit(combined_results$pvalue, ","), function(x) max(as.numeric(x)))

combined_results <- combined_results[,c(1,2,7,5,6)]
combined_results
# # save
# write.table(combined_results, "table/GWAS_results/NonTradLipid_sum_ratio_common_genes_lowinput_individual_logpvalue_7_spat_fitted_BLUP.txt", row.names=F, quote=F, sep="\t")

write.table(combined_results, "table/GWAS_results/Sum_Ratio_lipids_indv_sum_ratio_common_genes_lowinput_individual_logpvalue_7_spat_fitted_BLUP.txt", row.names=F, quote=F, sep="\t")
getwd()
# write.table(combined_results, "table/GWAS_results/SQDG_32_0_common_genes_lowinput_individual_logpvalue_7_spat_fitted_BLUP.txt", row.names=F, quote=F, sep="\t")
