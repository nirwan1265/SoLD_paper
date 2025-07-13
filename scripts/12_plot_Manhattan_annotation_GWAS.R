# Load necessary libraries
library(vroom)
library(dplyr)
library(CMplot)
library(GenomicRanges)
library(rtracklayer)


# Conda env
#conda activate /usr/local/usrapps/maize/ntanduk/seqanal

# Set the base directory where the phenotype folders are located
base_dir_LMM <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/sum_ratio_BLUP/raw_gwas"

# List all .txt files in the directory
files <- list.files(path=base_dir_LMM, pattern = "\\.txt$", full.names = TRUE)

# Combine the two lists
file_all <- files


# Function to create Manhattan plots and annotate them
create_manhattan_and_annotate <- function(file, ref_GRanges, output_prefix) {
  # Extract the file name without the directory and .txt extension
  file_name <- gsub(".txt$", "", basename(file))
  
  # Read the file into a data frame
  data <- vroom(file)
  
  # Select and rename columns
  data <- data %>% dplyr::select("rs", "chr", "ps", "p_wald")
  
  colnames(data) <- c("SNP", "Chromosome", "Position", "PValue")
  data_annotate <- data[which(data$PValue <= 0.05), ]
  data$SNP <- paste0("SNP_", data$Position)
  max_pval <- -log10(min(data$PValue)) + 1
  #Create Manhattan plot
  # CMplot(data, plot.type="m", col=c("grey30", "grey60"), LOG10=TRUE, ylim=c(2, max_pval), threshold=c(1e-7, 1e-5),
  #        threshold.lty=c(1, 2), threshold.lwd=c(1, 1), threshold.col=c("black", "grey"), amplify=TRUE,
  #        chr.den.col=NULL, signal.col=c("red", "green"), signal.cex=c(0.5, 0.5), signal.pch=c(19, 19),
  #        file="jpg", file.name=paste0(output_prefix, "_manhattan"), dpi=300, file.output=TRUE, verbose=TRUE, width=14, height=6, band=0, cex=0.5)
  # 
  # Create GRanges object
  gr <- GRanges(
    seqnames = Rle(data_annotate$Chromosome),
    ranges = IRanges(start = data_annotate$Position, end = data_annotate$Position),
    marker = data_annotate$SNP,
    pvalue = data_annotate$PValue
  )
  
  # Filter and resize SNPs for annotation
  snp_GRanges_filtered <- gr[mcols(gr)$pvalue < 0.05]
  snp_GRanges_filtered <- gr
  expanded_snp_GRanges <- resize(snp_GRanges_filtered, width = 25001, fix = "center")
  
  # Find overlaps with reference genes
  overlaps <- findOverlaps(genes_only, expanded_snp_GRanges)
  overlap_data <- data.frame(
    Chromosome = seqnames(expanded_snp_GRanges[subjectHits(overlaps)]),
    GeneID = mcols(ref_GRanges)[queryHits(overlaps), "ID"],
    SNP = mcols(expanded_snp_GRanges)[subjectHits(overlaps), "marker"],
    SNP_Position = start(snp_GRanges_filtered[subjectHits(overlaps)]),
    Gene_Start = start(ref_GRanges[queryHits(overlaps)]),
    Gene_End = end(ref_GRanges[queryHits(overlaps)]),
    PValue = mcols(snp_GRanges_filtered)[subjectHits(overlaps), "pvalue"]
  )
  
  # Categorize SNP position relative to the gene
  overlap_data$Relation <- ifelse(
    overlap_data$Gene_Start <= overlap_data$SNP_Position & overlap_data$Gene_End >= overlap_data$SNP_Position,
    "within",
    ifelse(
      overlap_data$SNP_Position < overlap_data$Gene_Start,
      "upstream",
      "downstream"
    )
  )
  
  # Calculate distance from gene
  overlap_data$DistanceToGene <- ifelse(
    overlap_data$Relation == "within",
    "within",
    ifelse(
      overlap_data$Relation == "upstream",
      overlap_data$Gene_Start - overlap_data$SNP_Position,
      overlap_data$SNP_Position - overlap_data$Gene_End
    )
  )
  # remove "gene:" string from the second column of overlap_data
  overlap_data$GeneID <- gsub("gene:", "", overlap_data$GeneID)
  
  # Only select rows starting with SORBI_ in GeneID
  overlap_data <- overlap_data[grep("^SORBI_", overlap_data$GeneID), ]
  
  
  collapsed_data <- overlap_data %>%
    group_by(GeneID) %>%
    summarise(
      Chromosome = unique(Chromosome),
      SNPs = paste(unique(SNP), collapse = ","),
      SNP_Positions = paste(unique(SNP_Position), collapse = ","),
      Gene_Start = unique(Gene_Start),
      Gene_End = unique(Gene_End),
      PValues = paste(unique(as.character(PValue)), collapse = ","),
      Relation = paste(unique(Relation), collapse = ","),
      DistanceToGene = paste(unique(DistanceToGene), collapse = ","),
      MinPValue = min(as.numeric(PValue), na.rm = TRUE),  # new column for sorting
      .groups = "drop"
    ) %>%
    arrange(MinPValue)  # sort by minimum P-value (smallest = strongest association)
  
  # Write annotated data to CSV
  write.table(collapsed_data, paste0(output_prefix, "_annotation.txt"), row.names = FALSE, quote = FALSE)
}

# Load your reference GFF file - replace with the actual path
ref_GRanges <- rtracklayer::import("/rsstu/users/r/rrellan/DOE_CAREER/SorghumGEA/data/SAP/gene_annotation/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")

ref_GRanges <- rtracklayer::import("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Sorghum.annotation/ensemblgenomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3")

# Filter ref_GRanges for only genes
genes_only <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]

# Loop through each file, create Manhattan plots, and annotate
for (lipids in file_all) {
  # Extract the file name without the directory and .txt extension for output prefix
  file_name <- gsub(".txt$", "", basename(lipids))
  
  # Create Manhattan plot and annotate
  create_manhattan_and_annotate(lipids, genes_only, file_name)
}

