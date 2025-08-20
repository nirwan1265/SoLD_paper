library(vroom)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(CMplot)
library(GenomicRanges)

################################################################################
####O-acyltransferase#####
################################################################################

# 1) point to your folder
dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas/O-acyltransferase"

# 2) grab all .txt files, name them by the basename sans “.txt”
files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
names(files) <- basename(files) %>% str_remove("\\.txt$")

# 3) read the first file as “master” for rs/chr/ps
master <- vroom(files[[1]]) %>% 
  dplyr::select(rs, chr, ps)

# 4) for each file, read just rs + p_wald, rename p_wald → file‑name
pval_list <- imap(files, ~ vroom(.x) %>%
                    dplyr::select(rs, p_wald) %>%
                    dplyr::rename(!!.y := p_wald)
)

# 5) find the SNPs common to *all* files
common_snps <- reduce(pval_list, inner_join, by = "rs") %>% pull(rs)

# 6) build your final table:
final_df <- master %>%
  filter(rs %in% common_snps) %>%      # keep only shared SNPs
  inner_join(reduce(pval_list, inner_join, by = "rs"), by = "rs") %>%
  select(rs, chr, ps, everything())

# 7) sanity check
dim(final_df)      # should be ~27k rows × (3 + number_of_files) cols
head(final_df)

final_df <- na.omit(final_df)

colnames(final_df) <- c("SNP","Chromosome","Position",
                        "TG_58_3","TG_58_4","TG_58_5","TG_60_4","TG_60_5")


final_df <- final_df %>%
  mutate(
    SNP        = factor(SNP),
    Chromosome = factor(Chromosome),
    Position   = as.integer(Position)
  )

# verify
str(final_df)


# 1) Define your threshold and the trait column names
thr    <- 1e-7
traits <- c("TG_58_3","TG_58_4","TG_58_5","TG_60_4","TG_60_5")


# 2) Build a list where each element is the vector of SNPs passing that trait’s cutoff
highlight_list <- lapply(traits, function(tr) {
  final_df[final_df[[tr]] < thr, "SNP"]
})

# 3) Name each element after its trait
names(highlight_list) <- traits


# 1) Pull each element’s SNP vector as character
snp_vecs <- lapply(highlight_list, function(df) as.character(df$SNP))

# 2) Intersect them all to get the common SNPs
common_SNPs <- Reduce(intersect, snp_vecs)

# 3) (Optional) If you want to keep only those in each sub‑list:
highlight_list_common <- lapply(highlight_list, function(df) {
  df %>% filter(as.character(SNP) %in% common_SNPs)
})


SNPs <- c("SNP_50088993")
SNPs <- as.factor(SNPs)

# SNPs <-  final_df[
#   final_df$Sum_TG < 1e-7 |
#     final_df$TG_54_6 < 1e-7 |
#     final_df$TG_54_7 < 1e-7 |
#     final_df$TG_56_4 < 1e-7 |
#     final_df$TG_56_6 < 1e-7, 1]
# 

CMplot(final_df,type="p",plot.type="m",LOG10=TRUE,highlight=SNPs,highlight.type="l",
       threshold=1e-7,threshold.col="black",threshold.lty=1,col=c("grey80","black"),
       signal.cex=0.3, signal.col="red", highlight.col="green",highlight.cex=0.5,
       file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

?CMplot()

################################################################################
#### SQDG(32:0) #####
################################################################################
# 
# # 1) point to your folder
# dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas/SQDG_32_0/"
# 
# # 2) grab all .txt files, name them by the basename sans “.txt”
# files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
# names(files) <- basename(files) %>% str_remove("\\.txt$")
# 
# # 3) read the first file as “master” for rs/chr/ps
# main <- vroom(files) %>%
#   dplyr::select(c(2, 1, 3, 12)) %>%
#   mutate(rs = paste0(rs,"_", chr))  # <- new unique ID
# colnames(main) <- c("SNP", "Chromosome", "Position","SQDG(32:0)")
# 
# 
# # Load the gff3 file
# # ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")
# # 
# # # Filter ref_GRanges for only genes
# # genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]
# 
# 
# # Define buffer size
# buffer <- 25000  # 25 kb
# 
# # Step 1: Extract gene info
# target_gene <- "SORBI_3002G000600"
# gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]
# 
# gene_chr <- as.character(seqnames(gene_row))
# gene_start <- start(gene_row)
# gene_end <- end(gene_row)
# 
# # Step 2: Subset main to same chromosome
# main_chr <- main[main$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]
# 
# # Step 3: Get SNPs within ±25 kb of gene boundaries
# snps_in_window <- main_chr %>%
#   filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
#   arrange(Position)
# 
# # Final output
# snps_in_window
# 
# # Update SNPs in `snps_in_window` with new ID format
# #snps_in_window$SNP <- paste0("SNP_", snps_in_window$Chromosome, "_", snps_in_window$Position)
# 
# # Match main again to character just in case
# main$SNP <- as.character(main$SNP)
# SNPs <- unique(snps_in_window$SNP)
# SNPs <- factor(SNPs)
# 
# CMplot(main,
#        plot.type = "m",
#        type = "p", 
#        LOG10 = TRUE,
#        col = c("grey80", "black"),          # Alternating chromosome colors
#        highlight = SNPs,                    # Vector of SNP IDs to highlight
#        highlight.col = "green",
#        highlight.cex = 1.5,
#        highlight.pch = 19,
#        signal.col = c("red"),    # Use same color for all points initially
#        signal.cex = 1.5,
#        threshold = 1e-7,
#        threshold.col = "red",               # Red line for threshold
#        threshold.lty = 1,
#        chr.border = FALSE,
#        dpi = 300,
#        file = "jpg",
#        file.output = TRUE,
#        verbose = TRUE,
#        width = 14,
#        height = 6)



################################################################################
#### Zeaxanthin #####
################################################################################

# 1) point to your folder
dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas/lowinput/Zeaxanthin/"

# 2) grab all .txt files, name them by the basename sans “.txt”
files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
names(files) <- basename(files) %>% str_remove("\\.txt$")

# 3) read the first file as “master” for rs/chr/ps
main <- vroom(files) %>%
  dplyr::select(c(2, 1, 3, 12)) %>%
  mutate(rs = paste0(rs,"_", chr))  # <- new unique ID
colnames(main) <- c("SNP", "Chromosome", "Position","Zeaxanthin")


# Load the gff3 file
# ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")
# 
# # Filter ref_GRanges for only genes
# genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


# Define buffer size
buffer <- 25000  # 25 kb

# Step 1: Extract gene info
target_gene <- c("SORBI_3002G033800")
gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]

gene_chr <- as.character(seqnames(gene_row))
gene_start <- start(gene_row)
gene_end <- end(gene_row)

# Step 2: Subset main to same chromosome
main_chr <- main[main$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]

# Step 3: Get SNPs within ±25 kb of gene boundaries
snps_in_window <- main_chr %>%
  filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
  arrange(Position)

# Final output
snps_in_window


# Step 1: Extract gene info
# target_gene <- c("SORBI_3002G033800")
# gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]
# 
# gene_chr <- as.character(seqnames(gene_row))
# gene_start <- start(gene_row)
# gene_end <- end(gene_row)
# 
# # Step 2: Subset main to same chromosome
# main_chr <- main[main$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]
# 
# # Step 3: Get SNPs within ±25 kb of gene boundaries
# snps_in_window2 <- main_chr %>%
#   filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
#   arrange(Position)
# 
# # Final output
# snps_in_window2
# 
# # Combine the 2 windoes
# snps_in_window <- bind_rows(snps_in_window, snps_in_window2)

# Update SNPs in `snps_in_window` with new ID format
#snps_in_window$SNP <- paste0("SNP_", snps_in_window$Chromosome, "_", snps_in_window$Position)

# Match main again to character just in case
main$SNP <- as.character(main$SNP)
SNPs <- unique(snps_in_window$SNP)
SNPs <- factor(SNPs)

CMplot(main,
       plot.type = "m",
       type = "p", 
       LOG10 = TRUE,
       col = c("grey80", "black"),          # Alternating chromosome colors
       highlight = SNPs,                    # Vector of SNP IDs to highlight
       highlight.col = "green",
       highlight.cex = 1.5,
       highlight.pch = 19,
       signal.col = c("red"),    # Use same color for all points initially
       signal.cex = 1.5,
       threshold = 1e-5,
       threshold.col = "red",               # Red line for threshold
       threshold.lty = 1,
       chr.border = FALSE,
       dpi = 300,
       file = "jpg",
       file.output = TRUE,
       verbose = TRUE,
       width = 14,
       height = 6)



################################################################################
#### beta_Sitosterol #####
################################################################################

# # 1) point to your folder
# dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas/beta__Sitosterol/"
# 
# # 2) grab all .txt files, name them by the basename sans “.txt”
# files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
# names(files) <- basename(files) %>% str_remove("\\.txt$")
# 
# # 3) read the first file as “master” for rs/chr/ps
# main <- vroom(files) %>%
#   dplyr::select(c(2, 1, 3, 12)) %>%
#   mutate(rs = paste0(rs,"_", chr))  # <- new unique ID
# colnames(main) <- c("SNP", "Chromosome", "Position","Beta Sitosterol")
# 
# 
# # Load the gff3 file
# # ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")
# # 
# # # Filter ref_GRanges for only genes
# # genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]
# 
# 
# # Define buffer size
# buffer <- 25000  # 25 kb
# 
# # Step 1: Extract gene info
# target_gene <- "SORBI_3003G049600"
# gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]
# 
# gene_chr <- as.character(seqnames(gene_row))
# gene_start <- start(gene_row)
# gene_end <- end(gene_row)
# 
# # Step 2: Subset main to same chromosome
# main_chr <- main[main$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]
# 
# # Step 3: Get SNPs within ±25 kb of gene boundaries
# snps_in_window <- main_chr %>%
#   filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
#   arrange(Position)
# 
# # Final output
# snps_in_window
# 
# # Update SNPs in `snps_in_window` with new ID format
# #snps_in_window$SNP <- paste0("SNP_", snps_in_window$Chromosome, "_", snps_in_window$Position)
# 
# # Match main again to character just in case
# main$SNP <- as.character(main$SNP)
# SNPs <- unique(snps_in_window$SNP)
# SNPs <- factor(SNPs)
# 
# CMplot(main,
#        plot.type = "m",
#        type = "p", 
#        LOG10 = TRUE,
#        col = c("grey80", "black"),          # Alternating chromosome colors
#        highlight = SNPs,                    # Vector of SNP IDs to highlight
#        highlight.col = "green",
#        highlight.cex = 1.5,
#        highlight.pch = 19,
#        signal.col = c("red"),    # Use same color for all points initially
#        signal.cex = 1.5,
#        threshold = 1e-5,
#        threshold.col = "red",               # Red line for threshold
#        threshold.lty = 1,
#        chr.border = FALSE,
#        dpi = 300,
#        file = "jpg",
#        file.output = TRUE,
#        verbose = TRUE,
#        width = 14,
#        height = 6)







################################################################################
#### Alpha_Carotene #####
################################################################################

# 1) point to your folder
dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas/lowinput/Alpha-Carotene/"

# 2) grab all .txt files, name them by the basename sans “.txt”
files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
names(files) <- basename(files) %>% str_remove("\\.txt$")

# 3) read the first file as “master” for rs/chr/ps
main <- vroom(files) %>%
  dplyr::select(c(2, 1, 3, 12)) %>%
  mutate(rs = paste0(rs,"_", chr))  # <- new unique ID
colnames(main) <- c("SNP", "Chromosome", "Position","Alpha Carotene")


# Load the gff3 file
# ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")
# 
# # Filter ref_GRanges for only genes
# genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


# Define buffer size
buffer <- 25000  # 25 kb

# Step 1: Extract gene info
target_gene <- "SORBI_3006G202500"
gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]

gene_chr <- as.character(seqnames(gene_row))
gene_start <- start(gene_row)
gene_end <- end(gene_row)

# Step 2: Subset main to same chromosome
main_chr <- main[main$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]

# Step 3: Get SNPs within ±25 kb of gene boundaries
snps_in_window <- main_chr %>%
  filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
  arrange(Position)

# Final output
snps_in_window

# Update SNPs in `snps_in_window` with new ID format
#snps_in_window$SNP <- paste0("SNP_", snps_in_window$Chromosome, "_", snps_in_window$Position)

# Match main again to character just in case
main$SNP <- as.character(main$SNP)
SNPs <- unique(snps_in_window$SNP)
SNPs <- factor(SNPs)

CMplot(main,
       plot.type = "m",
       type = "p", 
       LOG10 = TRUE,
       col = c("grey80", "black"),          # Alternating chromosome colors
       highlight = SNPs,                    # Vector of SNP IDs to highlight
       highlight.col = "green",
       highlight.cex = 1.5,
       highlight.pch = 19,
       signal.col = c("red"),    # Use same color for all points initially
       signal.cex = 1.5,
       threshold = 1e-5,
       threshold.col = "red",               # Red line for threshold
       threshold.lty = 1,
       chr.border = FALSE,
       dpi = 300,
       file = "jpg",
       file.output = TRUE,
       verbose = TRUE,
       width = 14,
       height = 6)





################################################################################
#### gibberellic_acid #####
################################################################################

# 1) point to your folder
dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas/lowinput/gibberellic_acid/"

# 2) grab all .txt files, name them by the basename sans “.txt”
files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
names(files) <- basename(files) %>% str_remove("\\.txt$")

# 3) read the first file as “master” for rs/chr/ps
main <- vroom(files) %>%
  dplyr::select(c(2, 1, 3, 12)) %>%
  mutate(rs = paste0(rs,"_", chr))  # <- new unique ID
colnames(main) <- c("SNP", "Chromosome", "Position","Gibberellic acid")


# Load the gff3 file
# ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")
# 
# # Filter ref_GRanges for only genes
# genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


# Define buffer size
buffer <- 25000  # 25 kb

# Step 1: Extract gene info
target_gene <- "SORBI_3007G090421"
gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]

gene_chr <- as.character(seqnames(gene_row))
gene_start <- start(gene_row)
gene_end <- end(gene_row)

# Step 2: Subset main to same chromosome
main_chr <- main[main$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]

# Step 3: Get SNPs within ±25 kb of gene boundaries
snps_in_window <- main_chr %>%
  filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
  arrange(Position)

# Final output
snps_in_window


# Step 1: Extract gene info
target_gene <- "SORBI_3010G030600"
gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]

gene_chr <- as.character(seqnames(gene_row))
gene_start <- start(gene_row)
gene_end <- end(gene_row)

# Step 2: Subset main to same chromosome
main_chr <- main[main$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]

# Step 3: Get SNPs within ±25 kb of gene boundaries
snps_in_window2 <- main_chr %>%
  filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
  arrange(Position)

# Final output
snps_in_window2

# Combine the 2 windows
snps_in_window <- bind_rows(snps_in_window, snps_in_window2)

# Update SNPs in `snps_in_window` with new ID format
#snps_in_window$SNP <- paste0("SNP_", snps_in_window$Chromosome, "_", snps_in_window$Position)

# Match main again to character just in case
main$SNP <- as.character(main$SNP)
SNPs <- unique(snps_in_window$SNP)
SNPs <- factor(SNPs)

CMplot(main,
       plot.type = "m",
       type = "p", 
       LOG10 = TRUE,
       col = c("grey80", "black"),          # Alternating chromosome colors
       highlight = SNPs,                    # Vector of SNP IDs to highlight
       highlight.col = "green",
       highlight.cex = 1.5,
       highlight.pch = 19,
       signal.col = c("red"),    # Use same color for all points initially
       signal.cex = 1.5,
       threshold = 1e-5,
       threshold.col = "red",               # Red line for threshold
       threshold.lty = 1,
       chr.border = FALSE,
       dpi = 300,
       file = "jpg",
       file.output = TRUE,
       verbose = TRUE,
       width = 14,
       height = 6)





################################################################################
#### DGDG/MGDG #####
################################################################################

# 1) point to your folder
dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas/lowinput/Sum_DGDG_Sum_MGDG/"

# 2) grab all .txt files, name them by the basename sans “.txt”
files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
names(files) <- basename(files) %>% str_remove("\\.txt$")

# 3) read the first file as “master” for rs/chr/ps
main <- vroom(files) %>%
  dplyr::select(c(2, 1, 3, 12)) %>%
  mutate(rs = paste0(rs,"_", chr))  # <- new unique ID
colnames(main) <- c("SNP", "Chromosome", "Position","DGDG_MGDG")


# Load the gff3 file
# ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")
# 
# # Filter ref_GRanges for only genes
# genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


# Define buffer size
buffer <- 25000  # 25 kb

# Step 1: Extract gene info
target_gene <- "SORBI_3005G045000"
gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]

gene_chr <- as.character(seqnames(gene_row))
gene_start <- start(gene_row)
gene_end <- end(gene_row)

# Step 2: Subset main to same chromosome
main_chr <- main[main$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]

# Step 3: Get SNPs within ±25 kb of gene boundaries
snps_in_window <- main_chr %>%
  filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
  arrange(Position)

# Final output
snps_in_window

# Update SNPs in `snps_in_window` with new ID format
#snps_in_window$SNP <- paste0("SNP_", snps_in_window$Chromosome, "_", snps_in_window$Position)

# Match main again to character just in case
main$SNP <- as.character(main$SNP)
SNPs <- unique(snps_in_window$SNP)
SNPs <- factor(SNPs)

CMplot(main,
       plot.type = "m",
       type = "p", 
       LOG10 = TRUE,
       col = c("grey80", "black"),          # Alternating chromosome colors
       highlight = SNPs,                    # Vector of SNP IDs to highlight
       highlight.col = "green",
       highlight.cex = 1.5,
       highlight.pch = 19,
       signal.col = c("red"),    # Use same color for all points initially
       signal.cex = 1.5,
       threshold = 1e-5,
       threshold.col = "red",               # Red line for threshold
       threshold.lty = 1,
       chr.border = FALSE,
       dpi = 300,
       file = "jpg",
       file.output = TRUE,
       verbose = TRUE,
       width = 14,
       height = 6)





################################################################################
#### FA(22:1) #####
################################################################################

# 1) point to your folder
dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas/lowinput/FA(22:1)/"

# 2) grab all .txt files, name them by the basename sans “.txt”
files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
names(files) <- basename(files) %>% str_remove("\\.txt$")

# 3) read the first file as “master” for rs/chr/ps
main <- vroom(files) %>%
  dplyr::select(c(2, 1, 3, 12)) %>%
  mutate(rs = paste0(rs,"_", chr))  # <- new unique ID
colnames(main) <- c("SNP", "Chromosome", "Position","FA(22:1)")


# Load the gff3 file
# ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")
# 
# # Filter ref_GRanges for only genes
# genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


# Define buffer size
buffer <- 25000  # 25 kb

# Step 1: Extract gene info
target_gene <- "SORBI_3001G097100"
gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]

gene_chr <- as.character(seqnames(gene_row))
gene_start <- start(gene_row)
gene_end <- end(gene_row)

# Step 2: Subset main to same chromosome
main_chr <- main[main$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]

# Step 3: Get SNPs within ±25 kb of gene boundaries
snps_in_window <- main_chr %>%
  filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
  arrange(Position)

# Final output
snps_in_window

# Update SNPs in `snps_in_window` with new ID format
#snps_in_window$SNP <- paste0("SNP_", snps_in_window$Chromosome, "_", snps_in_window$Position)

# Match main again to character just in case
main$SNP <- as.character(main$SNP)
SNPs <- unique(snps_in_window$SNP)
SNPs <- factor(SNPs)

CMplot(main,
       plot.type = "m",
       type = "p", 
       LOG10 = TRUE,
       col = c("grey80", "black"),          # Alternating chromosome colors
       highlight = SNPs,                    # Vector of SNP IDs to highlight
       highlight.col = "green",
       highlight.cex = 1.5,
       highlight.pch = 19,
       signal.col = c("red"),    # Use same color for all points initially
       signal.cex = 1.5,
       threshold = 1e-5,
       threshold.col = "red",               # Red line for threshold
       threshold.lty = 1,
       chr.border = FALSE,
       dpi = 300,
       file = "jpg",
       file.output = TRUE,
       verbose = TRUE,
       width = 14,
       height = 6)








SORBI_3002G282500

################################################################################
#### FA(22:1) #####
################################################################################

# # 1) point to your folder
# dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas/lowinput/"
# 
# # 2) grab all .txt files, name them by the basename sans “.txt”
# files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
# files <- files[2]
# names(files) <- basename(files) %>% str_remove("\\.txt$")
# 
# # 3) read the first file as “master” for rs/chr/ps
# main <- vroom(files) %>%
#   dplyr::select(c(2, 1, 3, 12)) %>%
#   mutate(rs = paste0(rs,"_", chr))  # <- new unique ID
# colnames(main) <- c("SNP", "Chromosome", "Position","Sum_PA")
# 
# # Remove NA
# main <- main %>%
#   filter(!is.na(SNP) & !is.na(Chromosome) & !is.na(Position) & !is.na(`DGMG`))
# main
# 
# # Load the gff3 file
# # ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")
# # 
# # # Filter ref_GRanges for only genes
# # genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]
# 
# 
# # Define buffer size
# buffer <- 25000  # 25 kb
# 
# # Step 1: Extract gene info
# target_gene <- "SORBI_3005G222500"
# gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]
# 
# gene_chr <- as.character(seqnames(gene_row))
# gene_start <- start(gene_row)
# gene_end <- end(gene_row)
# 
# # Step 2: Subset main to same chromosome
# main_chr <- main[main$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]
# 
# # Step 3: Get SNPs within ±25 kb of gene boundaries
# snps_in_window <- main_chr %>%
#   filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
#   arrange(Position)
# 
# # Final output
# snps_in_window
# 
# # Update SNPs in `snps_in_window` with new ID format
# #snps_in_window$SNP <- paste0("SNP_", snps_in_window$Chromosome, "_", snps_in_window$Position)
# 
# # Match main again to character just in case
# main$SNP <- as.character(main$SNP)
# SNPs <- unique(snps_in_window$SNP)
# SNPs <- factor(SNPs)
# 
# CMplot(main,
#        plot.type = "m",
#        type = "p", 
#        LOG10 = TRUE,
#        col = c("grey80", "black"),          # Alternating chromosome colors
#        highlight = SNPs,                    # Vector of SNP IDs to highlight
#        highlight.col = "green",
#        highlight.cex = 1.5,
#        highlight.pch = 19,
#        signal.col = c("red"),    # Use same color for all points initially
#        signal.cex = 1.5,
#        threshold = 1e-5,
#        threshold.col = "red",               # Red line for threshold
#        threshold.lty = 1,
#        chr.border = FALSE,
#        dpi = 300,
#        file = "jpg",
#        file.output = TRUE,
#        verbose = TRUE,
#        width = 14,
#        height = 6)







################################################################################
#### Sum SQDG
################################################################################

# 1) point to your folder
dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas/lowinput"

# 2) grab all .txt files, name them by the basename sans “.txt”
files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
files <- files[2]
names(files) <- basename(files) %>% str_remove("\\.txt$")

# 3) read the first file as “master” for rs/chr/ps
main <- vroom(files) %>%
  dplyr::select(c(2, 1, 3, 12)) %>%
  mutate(rs = paste0(rs,"_", chr))  # <- new unique ID
colnames(main) <- c("SNP", "Chromosome", "Position","Sum_SQDG")

# Remove NA
main <- main %>%
  filter(!is.na(SNP) & !is.na(Chromosome) & !is.na(Position) & !is.na(`Sum_SQDG`))
main

# Load the gff3 file
ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")

# Filter ref_GRanges for only genes
genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


# Define buffer size
buffer <- 25000  # 25 kb

# Step 1: Extract gene info
target_gene <- "SORBI_3005G195600"
gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]

gene_chr <- as.character(seqnames(gene_row))
gene_start <- start(gene_row)
gene_end <- end(gene_row)

# Step 2: Subset main to same chromosome
main_chr <- main[main$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]

# Step 3: Get SNPs within ±25 kb of gene boundaries
snps_in_window <- main_chr %>%
  filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
  arrange(Position)

# Final output
snps_in_window

# Update SNPs in `snps_in_window` with new ID format
#snps_in_window$SNP <- paste0("SNP_", snps_in_window$Chromosome, "_", snps_in_window$Position)

# Match main again to character just in case
main$SNP <- as.character(main$SNP)
SNPs <- unique(snps_in_window$SNP)
SNPs <- factor(SNPs)

CMplot(main,
       plot.type = "m",
       type = "p", 
       LOG10 = TRUE,
       col = c("grey80", "black"),          # Alternating chromosome colors
       highlight = SNPs,                    # Vector of SNP IDs to highlight
       highlight.col = "green",
       highlight.cex = 1.5,
       highlight.pch = 19,
       signal.col = c("red"),    # Use same color for all points initially
       signal.cex = 1.5,
       threshold = 1e-5,
       threshold.col = "red",               # Red line for threshold
       threshold.lty = 1,
       chr.border = FALSE,
       dpi = 300,
       file = "jpg",
       file.output = TRUE,
       verbose = TRUE,
       width = 14,
       height = 6)



################################################################################
#### Sum MG to Sum MGMG
################################################################################

# 1) point to your folder
dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas/lowinput"

# 2) grab all .txt files, name them by the basename sans “.txt”
files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
files <- files[1]
names(files) <- basename(files) %>% str_remove("\\.txt$")

# 3) read the first file as “master” for rs/chr/ps
main <- vroom(files) %>%
  dplyr::select(c(2, 1, 3, 12)) %>%
  mutate(rs = paste0(rs,"_", chr))  # <- new unique ID
colnames(main) <- c("SNP", "Chromosome", "Position","Sum_MG_Sum_MGDG")

# Remove NA
main <- main %>%
  filter(!is.na(SNP) & !is.na(Chromosome) & !is.na(Position) & !is.na(`Sum_MG_Sum_MGDG`))
main


# Remove SNPS = SNP_47069788_3 SNP_14421318_9 SNP_14421324_9 SNP_14421356_9 from main
main <- main %>%
  filter(!SNP %in% c("SNP_47069788_3", "SNP_14421318_9", "SNP_14421324_9", "SNP_14421356_9","SNP_14421346_9"))

# Load the gff3 file
ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")

# Filter ref_GRanges for only genes
genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


# Define buffer size
buffer <- 25000  # 25 kb

# Step 1: Extract gene info
target_gene <- "SORBI_3001G283700"
gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]

gene_chr <- as.character(seqnames(gene_row))
gene_start <- start(gene_row)
gene_end <- end(gene_row)

# Step 2: Subset main to same chromosome
main_chr <- main[main$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]

# Step 3: Get SNPs within ±25 kb of gene boundaries
snps_in_window <- main_chr %>%
  filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
  arrange(Position)

# Final output
snps_in_window

# Update SNPs in `snps_in_window` with new ID format
#snps_in_window$SNP <- paste0("SNP_", snps_in_window$Chromosome, "_", snps_in_window$Position)

# Match main again to character just in case
main$SNP <- as.character(main$SNP)
SNPs <- unique(snps_in_window$SNP)
SNPs <- factor(SNPs)

CMplot(main,
       plot.type = "m",
       type = "p", 
       LOG10 = TRUE,
       col = c("grey80", "black"),          # Alternating chromosome colors
       highlight = SNPs,                    # Vector of SNP IDs to highlight
       highlight.col = "green",
       highlight.cex = 1.5,
       highlight.pch = 19,
       signal.col = c("red"),    # Use same color for all points initially
       signal.cex = 1.5,
       threshold = 1e-5,
       threshold.col = "red",               # Red line for threshold
       threshold.lty = 1,
       chr.border = FALSE,
       dpi = 300,
       file = "jpg",
       file.output = TRUE,
       verbose = TRUE,
       width = 14,
       height = 6)




################################################################################
####Phospholipids#####
################################################################################

# 1) point to your folder
dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas/lowinput/Phospholipids/"


# 2) grab all .txt files, name them by the basename sans “.txt”
files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
names(files) <- basename(files) %>% str_remove("\\.txt$")
files <- files[c(1:4,8)]
files <- files[5:8]

# 3) read the first file as “master” for rs/chr/ps
master <- vroom(files[[1]]) %>% 
  dplyr::select(rs, chr, ps)

# 4) for each file, read just rs + p_wald, rename p_wald → file‑name
pval_list <- imap(files, ~ vroom(.x) %>%
                    dplyr::select(rs, p_wald) %>%
                    dplyr::rename(!!.y := p_wald)
)

# 5) find the SNPs common to *all* files
common_snps <-  purrr::reduce(pval_list, inner_join, by = "rs") %>% pull(rs)

# 6) build your final table:
final_df <- master %>%
  dplyr::filter(rs %in% common_snps) %>%      # keep only shared SNPs
  inner_join(purrr::reduce(pval_list, inner_join, by = "rs"), by = "rs") %>%
  dplyr::select(rs, chr, ps, everything())

# 7) sanity check
dim(final_df)      # should be ~27k rows × (3 + number_of_files) cols
head(final_df)

final_df <- na.omit(final_df)
colnames(final_df)



colnames(final_df) <- c("SNP","Chromosome","Position",
                        "PC(16:0_20:3)","PC(16:0_22:5)","PC(16:0_22:6)","PC(18:1_20:1)","PE(16:0_18:1)")

colnames(final_df) <- c("SNP","Chromosome","Position",
                        "PC(18:1_24:1)","PC(18:1_24:1)","PC(18:2_20:0)","PC(18:3_0:0)",)




final_df <- final_df %>%
  mutate(
    SNP        = factor(SNP),
    Chromosome = factor(Chromosome),
    Position   = as.integer(Position)
  )

# verify
str(final_df)


# 1) Define your threshold and the trait column names
thr    <- 1e-5
traits <- c("PC(16:0_20:3)","PC(16:0_22:5)","PC(16:0_22:6)","PC(18:1_20:1)","PE(16:0_18:1)")
traits <- c("PC(18:1_24:1)","PC(18:2_20:0)","PC(18:3_0:0)")


# Load the gff3 file
ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")

# Filter ref_GRanges for only genes
genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]


# Define buffer size
buffer <- 25000  # 25 kb

# Step 1: Extract gene info
target_gene <- "SORBI_3001G384300"
gene_row <- genes_only[mcols(genes_only)$gene_id == target_gene]

gene_chr <- as.character(seqnames(gene_row))
gene_start <- start(gene_row)
gene_end <- end(gene_row)

# Step 2: Subset main to same chromosome
main_chr <- final_df[final_df$Chromosome == as.numeric(gsub("chr", "", gene_chr)), ]

# Step 3: Get SNPs within ±25 kb of gene boundaries
snps_in_window <- main_chr %>%
  filter(Position >= (gene_start - buffer) & Position <= (gene_end + buffer)) %>%
  arrange(Position)

# Final output
snps_in_window


SNPs <- as.vector(snps_in_window$SNP)
# Remove duplictes
SNPs <- unique(SNPs)
SNPs <- as.factor(SNPs)

# SNPs <-  final_df[
#   final_df$Sum_TG < 1e-7 |
#     final_df$TG_54_6 < 1e-7 |
#     final_df$TG_54_7 < 1e-7 |
#     final_df$TG_56_4 < 1e-7 |
#     final_df$TG_56_6 < 1e-7, 1]
# 

CMplot(final_df,type="p",plot.type="m",LOG10=TRUE,highlight=SNPs,highlight.type="l",
       threshold=1e-7,threshold.col="black",threshold.lty=1,col=c("grey80","black"),
       signal.cex=0.3, signal.col="red", highlight.col="green",highlight.cex=0.5,
       file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

getwd()





