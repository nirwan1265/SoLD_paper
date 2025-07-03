# Load necessary library
library(dplyr)




# List all .txt files in the directory

# Control
# Individual
dir <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/lowinput/"
file_list <- list.files(path = dir,pattern = "*.txt")

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

p_value_threshold <- 0.00001
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
results$phenotype <- str_replace(results$phenotype, "mod_sub_lowinput_all_lipids_BLUPs_rename.part[0-9]+_SAP_bialleles_MAF_0.05.assoc_annotation", "")
#results$phenotype <- str_replace(results$phenotype, "mod_sub_lowinput_all_lipids_BLUPs_rename", "")


#results$phenotype <- str_replace(results$phenotype, "_mod_sub_summed_lipids_lowinput_with_ratios\\.part[0-9]*_SAP_bialleles_MAF_0.05\\.assoc_annotation", "")
results$phenotype <- str_replace(results$phenotype, "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/BLUP/lowinput/Traditional_lipids/", "")
#results$phenotype <- str_replace(results$phenotype, "_annotation", "")

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

# save
#write.table(combined_results, "LPC_lowinput_individual_logpvalue_7.txt", row.names=F, quote=F, sep="\t")
write.table(combined_results, "NonTradLipid_common_genes_lowinput_individual_logpvalue_5_BLUP.txt", row.names=F, quote=F, sep="\t")
#write.table(combined_results, "common_genes_control_individual_logpvalue_7.txt", row.names=F, quote=F, sep="\t")
getwd()
#write.table(combined_results, "VitaminK1_common_genes_control_individual_logpvalue_7.txt", row.names=F, quote=F, sep="\t")
#write.csv(combined_results, "common_genes_lowinput_individual_logpvalue_5.csv", row.names=F, quote=F)
getwd()
