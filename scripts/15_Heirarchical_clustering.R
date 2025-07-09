# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)  # For heatmap visualization
library(dendextend)       # For hierarchical clustering visualization

### GET THE DATA:
sap_li <- vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/GWAS_phenotype/individual_lipids_lowinput.csv")
getwd()
ncol(sap_li)

sap_li <- vroom("data/SPATS_fitted/non_normalized_intensities/lowinput_all_lipids_fitted_phenotype_non_normalized.csv")

# Remove the first column (genotype IDs)
lipid_data <- sap_li %>% dplyr::select(-c(2,3,4))  
lipid_data <- as.data.frame(t(lipid_data))
colnames(lipid_data) <- lipid_data[1, ]  # Set the first row as column names
lipid_data <- lipid_data[-1, ]  # Remove the first row (now redundant)
str(lipid_data)

# Convert all columns to numeric
lipid_data <- lipid_data %>% mutate(across(everything(), as.numeric))

# Valid classes
# valid_classes <- c("TG","DG","MG",
#                    "PC","PE","PG","PI",
#                    "LPC","LPE",
#                    "DGDG","MGDG","GalCer",
#                    "Cer","SM","FA","DGDG","SQDG","AEG","GalCer","FA")


valid_classes <- c("TG","DG","MG",
                   "PC","PE","PG","PI",
                   "LPC","LPE",
                   "DGDG","MGDG",
                   "SQDG")


class_pattern <- paste0("^(", paste(valid_classes, collapse = "|"), ")")

# 3) Find matching column names
lipid_cols <- grep(class_pattern, names(lipid_data), value = TRUE)
lipid_data <- lipid_data %>% select(all_of(lipid_cols))

# Extract lipid classes from column names
# Compute the Pearson correlation matrix
cor_matrix <- cor(lipid_data, use = "pairwise.complete.obs", method = "pearson")

#### PERFORM HIERARCHICAL CLUSTERING
# Convert correlation into distance (1 - correlation)
dist_matrix <- stats::as.dist(1 - cor_matrix)

# Hierarchical clustering using Ward's method
hc <- hclust(dist_matrix, method = "ward.D2")

#### DEFINE CLUSTERS
# Cut the dendrogram at a threshold (adjust as needed)
cluster_assignments <- cutree(hc, h = 0.9)  # Adjust height to change grouping

# Create a dataframe of lipids and their cluster assignments
cluster_df <- data.frame(
  Lipid = colnames(lipid_data),
  Cluster = cluster_assignments
)

# Step 4: Save or Display the Clusters
print(cluster_df)  # View the grouped lipids
table(cluster_df$Cluster)
#write.csv(cluster_df, "lipid_clusters_lowinput.csv", row.names = FALSE)  # Save to CSV

#### VISUALIZE THE CLUSTERING
# Plot dendrogram with cluster assignments
dend <- as.dendrogram(hc) %>% color_branches(k = length(unique(cluster_assignments)))
quartz()
plot(dend)

# Open a PNG device to save the plot
png("Genotype_Correlation_Dendrogram_lowinput.png", width = 44, height = 12, res = 150, units = "in")

# Plot the dendrogram
plot(dend, main = "Hierarchical Clustering of Lipid Correlations")

# Close the device to save the file
dev.off()

# Open a PNG device to save the plot
png("Genotype_Correlation_Heatmap_lowinput.png", width = 44, height = 44, res = 150, units = "in")

# Heatmap of the correlation matrix
#quartz()
plot(Heatmap(cor_matrix, name = "Correlation", show_row_names = TRUE, show_column_names = TRUE,
        cluster_rows = hc, cluster_columns = hc))

dev.off()

# save
png("Lipid_Correlation_Heatmap_lowinput.png", width = 24, height = 24, res = 150, units = "in")
draw(heatmap_li)
dev.off()



#### SUBSETING THE DATA
# Ensure cluster_df has the right lipid names
colnames(cluster_df) <- c("Lipid", "Cluster")

# Reshape sap_li: Extract genotype column and pivot lipids into long format
sap_li_long <- sap_li %>%
  pivot_longer(cols = -1, names_to = "Lipid", values_to = "Value")

# Merge the data with cluster assignments
sap_li_clustered <- left_join(sap_li_long, cluster_df, by = "Lipid")

# Split dataset by cluster and save each one separately
clusters <- unique(sap_li_clustered$Cluster)

for (cluster in clusters) {
  # Extract data for the current cluster
  cluster_data <- sap_li_clustered %>% filter(Cluster == cluster)
  
  # Convert back to wide format (genotypes as rows, lipids as columns)
  cluster_wide <- cluster_data %>%
    dplyr::select(-Cluster) %>%
    pivot_wider(names_from = "Lipid", values_from = "Value")
  
  # Save the cluster as a CSV file
  write.csv(cluster_wide, paste0("Cluster_", cluster, "_lowinput.csv"), row.names = FALSE)
}


## PCA for cluster 6
cluster6 <- vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/results/Clustering_analysis/lowinput/cluster_gwas/Cluster_6_lowinput.csv")

# Remove genotype column (assuming it's the first column)
cluster6_numeric <- cluster6 %>% dplyr::select(-1)

# Perform PCA
pca_result <- prcomp(cluster6_numeric, center = TRUE, scale. = TRUE)

# Extract variance explained by each PC
explained_var <- summary(pca_result)$importance[2, ] * 100  # Convert to percentage

# Save PC1 as a new dataset
cluster6_PCA <- data.frame(
  Genotype = cluster6$...1,  # Genotype column
  PC1 = pca_result$x[, 1]    # First principal component
)

# Save the PC1 values
write.csv(cluster6_PCA, "Cluster6_PC1.csv", row.names = FALSE)

# Print variance explained by PC1
cat("Variance explained by PC1:", round(explained_var[1], 2), "%\n")

# Check if we should use multiple PCs
if (sum(explained_var[1:3]) > 80) {
  cat("PC1–PC3 explain >80% variance, consider multi-trait GWAS.\n")
  
  # Save PC1–PC3
  cluster6_PCs <- data.frame(
    Genotype = cluster6$...1,
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    PC3 = pca_result$x[, 3]
  )
  
  # Save the multiple PCs
  write.csv(cluster6_PCs, "Cluster6_PCs.csv", row.names = FALSE)
}

# Plot variance explained (scree plot)
png("PCA_Scree_Cluster6.png", width = 1000, height = 700)
barplot(explained_var[1:10], names.arg = paste0("PC", 1:10), col = "steelblue",
        main = "PCA Scree Plot - Cluster 6", ylab = "Variance Explained (%)")
dev.off()

