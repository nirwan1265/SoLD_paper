#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
# PCA analysis
################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(SNPRelate)

###PCA analysis
#Need a vcf file format of the hapmap which is converted to GDS format
#TASSEL or PLINK is used for converting hapmap to VCF file format
#Need a directory to  create the gds file. If working on the server, we might need to define this before starting
setwd("/Users/nirwantandukar/Documents/Research/data/sorghum_genotype/SAP")

# Converting vcf to gds
vcf_ij <- snpgdsVCF2GDS("SAP_Chr10.vcf.gz", "SAP_chr10.gds", method = "copy.num.of.ref")

# get the gds fil
gdsfile <- snpgdsOpen("SAP_chr10.gds")


##LD-based SNP pruning
set.seed(1000)

# Try different LD thresholds for sensitivity analysis but read in a paper somewhere that 0.2 was used for GBJ
snpset <- snpgdsLDpruning(gdsfile, ld.threshold = 0.2)

# Get all the SNPs
snpset.id <- unlist(unname(snpset))

# Run PCA
pca <- snpgdsPCA(gdsfile, snp.id = snpset.id, num.thread = 10)

# Get the eigenvalues
pca$eigenval


# Save eigenvectors and percentages
pca_tab <- data.frame(
  sample.id = pca$sample.id,
  EV1 = pca$eigenvect[, 1],
  EV2 = pca$eigenvect[, 2],
  EV3 = pca$eigenvect[, 3],
  EV4 = pca$eigenvect[, 4],
  EV5 = pca$eigenvect[, 5],
  EV6 = pca$eigenvect[, 6],
  EV7 = pca$eigenvect[, 7],
  EV8 = pca$eigenvect[, 8],
  EV9 = pca$eigenvect[, 9],
  EV10 = pca$eigenvect[, 10],
  stringsAsFactors = FALSE
)

# NOTE: look at the scree plot and save the correct no. of PCA. which is 5 here
pca_tab
# Write to CSV
write.csv(pca_tab[,1:6], "PCA_SAP.csv", row.names = FALSE)
getwd()
# Scree plot data (percentage variance explained)
scree <- data.frame(
  PC = paste0("PC", 1:length(pca$varprop)),
  Variance = pca$varprop * 100
)
scree_clean <- scree$Variance[!is.na(scree$Variance)]
#write.csv(scree, "results/PCA_variance_explained_maize_romerro.csv", row.names = FALSE)



library(ggplot2)
library(ggrepel) # For clean labels

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
    color = ifelse(scree_data$PC == 3, "red", "steelblue"), # Red for PC3
    size = ifelse(scree_data$PC == 3, 5, 3) # Larger for PC3
  ) +
  # Reference line at PC3's variance level
  geom_hline(
    yintercept = scree_data$Variance[3], 
    color = "red", 
    linetype = "dashed", 
    linewidth = 0.8
  ) +
  # Label only PC3 explicitly
  geom_text_repel(
    aes(label = ifelse(PC == 5, 
                       paste0("PC3 (", round(Variance[3], 1), "%)"), 
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
    subtitle = "Red line indicates optimal cutoff at PC5"
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


library(ggplot2)

# Make scree plot
scree_plot <- ggplot(data.frame(PC = 1:length(pca$eigenval), Variance = pca$eigenval), aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Scree Plot of PCA", x = "Principal Component", y = "Eigenvalue") +
  theme_minimal(base_size = 16)


# Convert PCA results to a data frame
pc_df <- data.frame(
  sample.id = pca$sample.id,  # Sample IDs
  PC1 = pca$eigenvect[,1],    # First Principal Component
  PC2 = pca$eigenvect[,2]     # Second Principal Component
)

# View the first few rows
head(pc_df)

# Add a column to categorize samples by prefix
pc_df$group <- case_when(
  grepl("^I01", pc_df$sample.id) ~ "I0",
  grepl("^I14", pc_df$sample.id) ~ "I14",
  grepl("^J01", pc_df$sample.id) ~ "J0",
  grepl("^J14", pc_df$sample.id) ~ "J14",
  TRUE ~ "Other"
)


# PCA Scatter Plot
# PCA Plot with Labels by Group
# PCA Plot without dot labels and fixed legend
pca_plot <- ggplot(pc_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(alpha = 0.8, size = 3, stroke = 0.5) +  # Slightly larger points with subtle stroke
  theme_minimal(base_size = 18) +
  labs(
    x = paste0("PC 1 (", round(pca$varprop[1] * 100, 2), "% variance)"),
    y = paste0("PC 2 (", round(pca$varprop[2] * 100, 2), "% variance)"),
    title = "PCA of Indian Chief and Jarvis Lines",
    subtitle = "Comparison between generations 0 and 14"
  ) +
  scale_color_manual(
    name = "Lines and Generation",
    values = c("I0" = "#E41A1C", "I14" = "#377EB8", "J0" = "#4DAF4A", "J14" = "#984EA3"),
    labels = c("I0" = "Indian Chief G0", "I14" = "Indian Chief G14", 
               "J0" = "Jarvis G0", "J14" = "Jarvis G14")
  ) +
  theme(
    text = element_text(family = "Arial"),  # Use a clean font
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold", color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 16, face = "bold", color = "black"),
    legend.position = "right",
    legend.key = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 20)),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +  # Make legend symbols more visible
  coord_fixed(ratio = 1)  # Keep aspect ratio square for proper PCA interpretation

quartz()
pca_plot



pca_plot <- ggplot(pc_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(alpha = 0.8, size = 3, stroke = 0.5) +
  theme_minimal(base_size = 18) +
  labs(
    x = paste0("PC 1 (", round(pca$varprop[1] * 100, 2), "% variance)"),
    y = paste0("PC 2 (", round(pca$varprop[2] * 100, 2), "% variance)"),
    title = "PCA of Indian Chief and Jarvis Lines",
    subtitle = "Comparison between generations 0 and 14"
  ) +
  scale_color_manual(
    name = "Line and Generation",
    values = c("I0" = "#E41A1C", "I14" = "#377EB8", "J0" = "#4DAF4A", "J14" = "#984EA3"),
    labels = c("I0" = "Indian Chief G0", "I14" = "Indian Chief G14", 
               "J0" = "Jarvis G0", "J14" = "Jarvis G14")
  ) +
  theme(
    text = element_text(family = "Arial"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black"),  # Slightly smaller for space
    legend.title = element_text(size = 14, face = "bold", color = "black"),  # Slightly smaller
    legend.position = c(0.95, 0.05),  # Bottom right position
    legend.justification = c(1, 0),  # Anchored to bottom right
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.key = element_blank(),
    legend.spacing = unit(0.2, "cm"),
    legend.box.margin = margin(5, 5, 5, 5),  # Small margin around legend
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 16, hjust = 0.5, margin = margin(b = 20)),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +  # Slightly smaller legend symbols
  coord_fixed(ratio = 1)



# Save the plot
ggsave("PCA_Indian_Jarvis.png", plot = pca_plot, width = 10, height = 8, dpi = 300)

