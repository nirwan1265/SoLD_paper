################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PCA ANALYSIS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################

# Libraries
library(patchwork) # For combining plots
library(SNPRelate)

# Set working directory where you want GDS created
setwd("/Users/nirwantandukar/Documents")

# Define input VCF file (can be .vcf or .vcf.gz)
vcf_file <- "/Users/nirwantandukar/Documents/SAP_Chr10.vcf.gz"

# Convert VCF to GDS
snpgdsVCF2GDS(vcf.fn = vcf_file, 
              out.fn = "SAP_chr10.gds", 
              method = "copy.num.of.ref", 
              ignore.chr.prefix = "chr")

# Open GDS file
gds_file <- snpgdsOpen("SAP_chr10.gds")

# LD pruning
set.seed(1000)
#snpset <- snpgdsLDpruning(gds_file, ld.threshold = 0.2)
snpset <- snpgdsLDpruning(
  gds_file,
  ld.threshold = 0.2,
  slide.max.bp = 50000,  # 50 kb window
  slide.max.n = 1000     # (optional) max SNPs to compare in the window
)
snpset.id <- unlist(snpset)

# Run PCA
pca_result <- snpgdsPCA(gds_file, snp.id = snpset.id, num.thread = 10)

# Save eigenvectors and percentages
pca_tab <- data.frame(
  sample.id = pca_result$sample.id,
  EV1 = pca_result$eigenvect[, 1],
  EV2 = pca_result$eigenvect[, 2],
  EV3 = pca_result$eigenvect[, 3],
  EV4 = pca_result$eigenvect[, 4],
  EV5 = pca_result$eigenvect[, 5],
  EV6 = pca_result$eigenvect[, 6],
  EV7 = pca_result$eigenvect[, 7],
  EV8 = pca_result$eigenvect[, 8],
  EV9 = pca_result$eigenvect[, 9],
  EV10 = pca_result$eigenvect[, 10],
  stringsAsFactors = FALSE
)

# Write to CSV
write.csv(pca_tab, "PCA_SAP_chr10.csv", row.names = FALSE)

# Scree plot data (percentage variance explained)
scree <- data.frame(
  PC = paste0("PC", 1:length(pca_result$varprop)),
  Variance = pca_result$varprop * 100
)

#write.csv(scree, "results/PCA_variance_explained_maize_romerro.csv", row.names = FALSE)



library(ggplot2)
library(ggrepel) # For clean labels

# non NA scree
scree_clean <- scree$Variance[!is.na(scree$Variance)]


# Data prep
scree_data <- data.frame(
  PC = 1:length(scree_clean),
  Variance = scree_clean,
  Cumulative = cumsum(scree_clean)
)

# Scree Plot (Variance per PC)
scree_plot_pc3 <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  # Line and points
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(
    color = ifelse(scree_data$PC == 8, "red", "steelblue"), # Red for PC3
    size = ifelse(scree_data$PC == 8, 5, 3) # Larger for PC3
  ) +
  # Reference line at PC3's variance level
  geom_hline(
    yintercept = scree_data$Variance[8], 
    color = "red", 
    linetype = "dashed", 
    linewidth = 0.8
  ) +
  # Label only PC3 explicitly
  geom_text_repel(
    aes(label = ifelse(PC == 8, 
                       paste0("PC8 (", round(Cumulative[8], 1), "%)"), 
                       "")),
    nudge_y = 1.5,
    nudge_x = 0.3,
    size = 5,
    color = "red",
    fontface = "bold"
  ) +
  # Axis labels and title
  labs(
    x = "Principal Component (PC)",
    y = "Variance Explained (%)",
    title = "PCA Scree Plot: Maize Landrace Diversity",
    subtitle = "Red line indicates optimal cutoff at PC8"
  ) +
  # Theme adjustments
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "red", size = 11),
    axis.text = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
  )


quartz()
scree_plot_pc3

# Save the plot
ggsave("results/figures/scree_plot_maize_romerro.png", scree_plot_pc3, 
       width = 6, height = 5, dpi = 300)



# Cumulative
cumulative_plot <- ggplot(scree_data, aes(x = PC, y = Cumulative)) +
  geom_line(color = "darkorange", linewidth = 1) +
  geom_point(color = "darkorange", size = 3) +
  geom_text_repel(
    aes(label = ifelse(PC %% 5 == 0, paste0(round(Cumulative, 1), "%"), "")),
    nudge_y = 2, size = 4, color = "black"
  ) +
  labs(
    x = "Principal Component (PC)",
    y = "Cumulative Variance Explained (%)",
    title = "Cumulative Variance: Maize Landrace PCA"
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
  )

quartz()
cumulative_plot
# Save as PDF
ggsave("results/figures/cumulative_variance_nature_maize_romerro.png", cumulative_plot, width = 6, height = 4, dpi = 300)




combined_plot <- (scree_plot_pc3 / cumulative_plot) + 
  plot_annotation(tag_levels = 'A') # Add panel labels (A, B)

plot(combined_plot)
ggsave("results/figures/combined_pca_plots_nature_maize_romerror.png", combined_plot, width = 6, height = 8, dpi = 300)

