#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
###################### SOrghum Lipidomics Database (SOLD) ######################
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
################################################################################
#### DATA PREPROCESSING 
################################################################################
#-------------------------------------------------------------------------------



################################################################################
############################# LIBRARIES  #######################################
################################################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

################################################################################
########################## READ THE FILES  #####################################
################################################################################


# Loading datasets
A <- read.csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/SERRF_normalization/SERRF_Result_A_new/normalized by - SERRF.csv")
B <- read.csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/SERRF_normalization/SERRF_Result_B_new/normalized by - SERRF.csv")
glimpse(A)
###### Filtering just the Peak intensities 
# Find column names that match the pattern "S1_Run" at the end
#selected_columns_A <- grep("S1_Run\\d+", colnames(A), value = TRUE)
#selected_columns_B <- grep("S1_Run\\d+", colnames(B), value = TRUE)

# Subset the CHECKS; for the checks, use the SERRF_Result_B_new else use the old one
selected_columns_A_CHECKS <- grep("CHECK", colnames(A), value = TRUE)
selected_columns_B_CHECKS <- grep("Check", colnames(B), value = TRUE)

# Subset the PI columns
selected_columns_A_PI <- grep("PI", colnames(A), value = TRUE)
selected_columns_B_PI <- grep("PI", colnames(B), value = TRUE)


# Combine the two
selected_columns_A <- c(selected_columns_A_CHECKS, selected_columns_A_PI)
selected_columns_B <- c(selected_columns_B_CHECKS, selected_columns_B_PI)

# Select the corresponding columns in your data frame A
A_filtered <- A[, colnames(A) %in% selected_columns_A]
B_filtered <- B[, colnames(B) %in% selected_columns_B]


# Adding colnames
rownames(A_filtered) <- A[,1]
rownames(B_filtered) <- B[,1]
rm(A,B)

# Saving CHECKs (removed in 2)
#write.csv(A_filtered, "Control_SERRF_checks.csv", row.names = TRUE)
#write.csv(B_filtered, "LowInput_SERRF_checks.csv", row.names = TRUE)


# Set the threshold for proportion of zeroes you want to remove
threshold <- 0.5

# Calculate the proportion of zeroes in each column
zero_proportions_A <- colMeans(A_filtered == 0, na.rm = TRUE)
zero_proportions_B <- colMeans(B_filtered == 0, na.rm = TRUE)

# Get the column indices to keep (where less than 50% are zeroes)
columns_to_keep_A <- which(zero_proportions_A < threshold)
columns_to_keep_B <- which(zero_proportions_B < threshold)

# Subset A_filtered to include only the selected columns
A_filtered <- A_filtered[, columns_to_keep_A]
A_filtered$X.Scan. <- rownames(A_filtered)
B_filtered <- B_filtered[, columns_to_keep_B]
B_filtered$X.Scan. <- rownames(B_filtered)


################################################################################
####################### OBTAINING LIPID NAMES  #################################
################################################################################

##### Filtering based on lipids
lipid_A <- read.csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/lipid_class_A.csv")
lipid_A$X.Scan. <- as.character(lipid_A$X.Scan.)
lipid_B <- read.csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/lipid_class_B.csv")
lipid_B$X.Scan. <- as.character(lipid_B$X.Scan.)

colnames(lipid_A)

# Unique ID's  
unique_scans_all_lipid_A <- lipid_A %>%
  dplyr::select(Compound_Name,X.Scan.)
unique_scans_all_lipid_B <- lipid_B %>%
  dplyr::select(Compound_Name,X.Scan.)

# Process duplicates
unique_scans_all_lipid_A$Compound_Name <- with(unique_scans_all_lipid_A, ave(Compound_Name, Compound_Name, FUN = function(x) {
  if (length(x) > 1) {
    return(paste0(x, "_", 1:length(x)))
  } else {
    return(x)
  }
}))

unique_scans_all_lipid_B$Compound_Name <- with(unique_scans_all_lipid_B, ave(Compound_Name, Compound_Name, FUN = function(x) {
  if (length(x) > 1) {
    return(paste0(x, "_", 1:length(x)))
  } else {
    return(x)
  }
}))


# Convert X.Scan. column to character type
unique_scans_all_lipid_A$X.Scan. <- as.character(unique_scans_all_lipid_A$X.Scan.)
unique_scans_all_lipid_B$X.Scan. <- as.character(unique_scans_all_lipid_B$X.Scan.)

# Subset rows in A_filtered using row names from unique_scans
subset_A <- inner_join(A_filtered, unique_scans_all_lipid_A)

subset_B <- inner_join(B_filtered, unique_scans_all_lipid_B)

# Collapse by X.Scan. and combine Compound_Name values
subset_A_collapsed <- subset_A %>% 
  dplyr::group_by(X.Scan.) %>% 
  dplyr::summarise(
    ## keep (or sum) numeric sample columns
    across(where(is.numeric), dplyr::first),
    
    ## main name column
    Compound_Name = paste(unique(Compound_Name), collapse = "%"),
    
    ## NP‑classifier taxonomy
    npclassifier_superclass = paste(unique(npclassifier_superclass), collapse = ";"),
    npclassifier_class      = paste(unique(npclassifier_class),      collapse = ";"),
    npclassifier_pathway    = paste(unique(npclassifier_pathway),    collapse = ";"),
    
    .groups = "drop"
  )


subset_B_collapsed <- subset_B %>% 
  dplyr::group_by(X.Scan.) %>% 
  dplyr::summarise(
    across(where(is.numeric), dplyr::first),
    Compound_Name = paste(unique(Compound_Name), collapse = "%"),
    npclassifier_superclass = paste(unique(npclassifier_superclass), collapse = ";"),
    npclassifier_class      = paste(unique(npclassifier_class),      collapse = ";"),
    npclassifier_pathway    = paste(unique(npclassifier_pathway),    collapse = ";"),
    .groups = "drop"
  )


# Remove the X.Scan. and order
subset_A <- subset_A_collapsed %>%            # reorder: name first, then everything
  dplyr::select(-X.Scan.) %>% 
  dplyr::select(
    Compound_Name,
    npclassifier_superclass,
    npclassifier_class,
    npclassifier_pathway,
    everything()
  )

subset_B <- subset_B_collapsed %>% 
  dplyr::select(-X.Scan.) %>% 
  dplyr::select(
    Compound_Name,
    npclassifier_superclass,
    npclassifier_class,
    npclassifier_pathway,
    everything()
  )

# Save the subset data
write.csv(subset_A, "Control_All_Lipids.csv", row.names = FALSE)
write.csv(subset_B, "LowInput_All_Lipids.csv", row.names = FALSE)



# Data wrangle
# Go to excel and remove Spectra Match to and NIST14
# Combine the Lipids C and Double bonds
# Only put 1 name if they have multiple names







################################################################################
############################# NORMALIZING  #####################################
################################################################################

# Read the data again
# subset_A <- read.csv("/Users/nirwantandukar/Documents/Research/results/SAP/Table/Control_All_Lipids.csv")
# subset_B <- read.csv("/Users/nirwantandukar/Documents/Research/results/SAP/Table/LowInput_All_Lipids.csv")
# str(subset_A)

subset_A <- read.csv("/Users/nirwantandukar/Documents/Research/results/SAP/Table/Control_all_lipids_raw_final.csv")
subset_B <- read.csv("/Users/nirwantandukar/Documents/Research/results/SAP/Table/Control_all_lipids_raw_final.csv")

# Averaging the repeats compounds
subset_A <- subset_A %>%
  dplyr::group_by(Compound_Name) %>%
  dplyr::summarize(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")


subset_B <- subset_B %>%
  dplyr::group_by(Compound_Name) %>%
  dplyr::summarize(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

# Save the table
write.csv(subset_A, "Control_All_Lipids_non_normailzed.csv", row.names = FALSE)
write.csv(subset_B, "LowInput_All_Lipids_non_normailzed.csv", row.names = FALSE)
getwd()
# NOTE GO TO EXCEL AND SEPARATE THE VALUES 
# Remove Mass, Collission energy and then add the Carbon and bonds for lipids











# Z-score normalization
subset_A_zscore <- subset_A %>%
  mutate(across(where(is.numeric), ~scale(.)[,1]))

subset_B_zscore <- subset_B %>%
  mutate(across(where(is.numeric), ~scale(.)[,1]))


# Log10 transform and median centering
subset_A_log10_median <- subset_A %>%
  mutate(across(where(is.numeric), ~log10(.) - median(log10(.), na.rm = TRUE)))

subset_B_log10_median <- subset_B %>%
  mutate(across(where(is.numeric), ~log10(.) - median(log10(.), na.rm = TRUE)))




# Define the compound of interest
compound <- ".beta.-Amyrin acetate"

# --------- Extract values ---------
# Raw
A_raw <- subset_A %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name) %>% unlist()
B_raw <- subset_B %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name) %>% unlist()

# Z-score
A_z <- subset_A_zscore %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name) %>% unlist()
B_z <- subset_B_zscore %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name) %>% unlist()

# Log10-median
A_log <- subset_A_log10_median %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name) %>% unlist()
B_log <- subset_B_log10_median %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name) %>% unlist()

# --------- Combine for plotting ---------
combine_for_plot <- function(A, B, label) {
  data.frame(
    Value = c(A, B),
    Group = rep(c("A", "B"), times = c(length(A), length(B))),
    Type = label
  )
}

df_all <- bind_rows(
  combine_for_plot(A_raw, B_raw, "Raw"),
  combine_for_plot(A_z, B_z, "Z-score"),
  combine_for_plot(A_log, B_log, "Log10-Median")
)

# --------- Boxplot ---------
quartz()
ggplot(df_all, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~Type, scales = "free_y") +
  theme_minimal() +
  ggtitle(paste("Comparison of", compound)) +
  ylab("Normalized Intensity") +
  theme(plot.title = element_text(hjust = 0.5))

# --------- T-tests ---------
t_raw <- t.test(A_raw, B_raw)
t_z <- t.test(A_z, B_z)
t_log <- t.test(A_log, B_log)

# Print results
cat("T-test Results for", compound, "\n")
cat("Raw: p =", t_raw$p.value, "\n")
cat("Z-score: p =", t_z$p.value, "\n")
cat("Log10-Median: p =", t_log$p.value, "\n")


compound <- ".beta.-Amyrin acetate"

# Extract and label each set
dist_df <- bind_rows(
  data.frame(Value = unlist(subset_A %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name)),
             Type = "Raw", Group = "A"),
  data.frame(Value = unlist(subset_B %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name)),
             Type = "Raw", Group = "B"),
  
  data.frame(Value = unlist(subset_A_zscore %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name)),
             Type = "Z-score", Group = "A"),
  data.frame(Value = unlist(subset_B_zscore %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name)),
             Type = "Z-score", Group = "B"),
  
  data.frame(Value = unlist(subset_A_log10_median %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name)),
             Type = "Log10-Median", Group = "A"),
  data.frame(Value = unlist(subset_B_log10_median %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name)),
             Type = "Log10-Median", Group = "B")
)

# Plot
quartz()
ggplot(dist_df, aes(x = Value, fill = Group)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Type, scales = "free") +
  theme_minimal() +
  ggtitle(paste("Density Plot of", compound)) +
  xlab("Value") +
  ylab("Density") +
  theme(plot.title = element_text(hjust = 0.5))



### QQ plot
compound <- ".beta.-Amyrin acetate"  # or any compound you want to analyze

# Extract all 6 versions
raw_A <- unlist(subset_A %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name))
z_A   <- unlist(subset_A_zscore %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name))
log_A <- unlist(subset_A_log10_median %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name))

raw_B <- unlist(subset_B %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name))
z_B   <- unlist(subset_B_zscore %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name))
log_B <- unlist(subset_B_log10_median %>% dplyr::filter(Compound_Name == compound) %>% dplyr::select(-Compound_Name))

# Open plotting window (macOS/quartz, otherwise use windows() on Windows)
quartz()

# Set up 3x2 layout
par(mfrow = c(2, 3))

# Generate QQ plots
qqnorm(raw_A, main = "QQ Plot - Raw (Group A)"); qqline(raw_A, col = "blue", lwd = 2)
qqnorm(z_A, main = "QQ Plot - Z-score (Group A)"); qqline(z_A, col = "blue", lwd = 2)
qqnorm(log_A, main = "QQ Plot - Log10-Median (Group A)"); qqline(log_A, col = "blue", lwd = 2)

qqnorm(raw_B, main = "QQ Plot - Raw (Group B)"); qqline(raw_B, col = "blue", lwd = 2)
qqnorm(z_B, main = "QQ Plot - Z-score (Group B)"); qqline(z_B, col = "blue", lwd = 2)
qqnorm(log_B, main = "QQ Plot - Log10-Median (Group B)"); qqline(log_B, col = "blue", lwd = 2)




### Seems Z-scores looks the best from the qqplots

# Save the Z-score data
# Transform the data
subset_A_zscore <- as.data.frame(t(subset_A_zscore))
subset_B_zscore <- as.data.frame(t(subset_B_zscore))

# Add the Compound_Name column
#subset_A_zscore$Compound_Name <- rownames(subset_A_zscore)
#subset_B_zscore$Compound_Name <- rownames(subset_B_zscore)

# Change Compound_Name row to Colnames and remove that row
colnames(subset_A_zscore) <- subset_A_zscore[1,]
subset_A_zscore <- subset_A_zscore[-1,]

colnames(subset_B_zscore) <- subset_B_zscore[1,]
subset_B_zscore <- subset_B_zscore[-1,]


# Save the data
# write.csv(subset_A_zscore, "Control_All_Lipids_Zscore.csv", row.names = T)
# write.csv(subset_B_zscore, "LowInput_All_Lipids_Zscore.csv", row.names = T)
# REMOVE THE LAST COLUMN


################################################################################
############################# NAME CHANGE  #####################################
################################################################################


### NOTE
# For raw samples df, use subset_A and B without Z scores
# this is used for Linex2 


# Load the data for the name change
# done using PubChem
name_change <- read.csv("/Users/nirwantandukar/Documents/Research/results/SAP/Lipid Class/lipid_class.csv",check.names = FALSE)
str(name_change)

# Common name binary
name_map <- name_change %>%
  dplyr::filter(CommonName != "") %>%
  dplyr::select(Compound_Name, CommonName) %>%
  deframe()  # creates named vector: names = old, values = new

# Load subset_A_zcore and subset_B_zscore
subset_A_zscore <- read.csv("/Users/nirwantandukar/Documents/Research/results/SAP/Table/Control_All_Lipids_Zscore.csv",
                            check.names = FALSE, header = F)
subset_A_zscore <- as.data.frame(t(subset_A_zscore))
colnames(subset_A_zscore) <- subset_A_zscore[1,]
subset_A_zscore <- subset_A_zscore[-1,]
colnames(subset_A_zscore)[1] <- "Compound_Name"


subset_B_zscore <- read.csv("/Users/nirwantandukar/Documents/Research/results/SAP/Table/LowInput_All_Lipids_Zscore.csv",check.names = FALSE, header = F)
subset_B_zscore <- as.data.frame(t(subset_B_zscore))
colnames(subset_B_zscore) <- subset_B_zscore[1,]
subset_B_zscore <- subset_B_zscore[-1,]
colnames(subset_B_zscore)[1] <- "Compound_Name"



# Load the un-normalized raw data
# For the sake of not changing the names i am just gonna use the _zscore naming system
subset_A_zscore <- read.csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/non_normalized/Control_All_Lipids_non_normailzed.csv",
                            check.names = FALSE, header = T)
# subset_A_zscore <- as.data.frame(t(subset_A_zscore))
# colnames(subset_A_zscore) <- subset_A_zscore[1,]
# subset_A_zscore <- subset_A_zscore[-1,]
# colnames(subset_A_zscore)[1] <- "Compound_Name"

subset_B_zscore <- read.csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/non_normalized/LowInput_All_Lipids_non_normailzed.csv",check.names = FALSE, header = T)
# subset_B_zscore <- as.data.frame(t(subset_B_zscore))
# colnames(subset_B_zscore) <- subset_B_zscore[1,]
# subset_B_zscore <- subset_B_zscore[-1,]
# colnames(subset_B_zscore)[1] <- "Compound_Name"

# Rownames NULL
rownames(subset_A_zscore) <- NULL
rownames(subset_B_zscore) <- NULL

str(subset_A_zscore)


### NOTE
# For raw samples
# subset_A_zscore <- subset_A
# subset_B_zscore <- subset_B


# Replace only those that exist in the map
subset_A_zscore$Compound_Name <- ifelse(
  subset_A_zscore$Compound_Name %in% names(name_map),
  name_map[subset_A_zscore$Compound_Name],
  subset_A_zscore$Compound_Name
)
str(subset_A_zscore)
subset_B_zscore$Compound_Name <- ifelse(
  subset_B_zscore$Compound_Name %in% names(name_map),
  name_map[subset_B_zscore$Compound_Name],
  subset_B_zscore$Compound_Name
)

# Convert all rows to numeric except the first column
subset_A_zscore[-1] <- lapply(subset_A_zscore[-1], function(x) as.numeric(as.character(x)))
subset_B_zscore[-1] <- lapply(subset_B_zscore[-1], function(x) as.numeric(as.character(x)))
str(subset_A_zscore)


# Average the rows with same name in Compound_Name
subset_A_zscore <- subset_A_zscore %>%
  dplyr::group_by(Compound_Name) %>%
  dplyr::summarize(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

subset_B_zscore <- subset_B_zscore %>%
  dplyr::group_by(Compound_Name) %>%
  dplyr::summarize(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")


# Save the data
write.csv(as.data.frame(t(subset_A_zscore)), "Control_all_lipids_final_non_normalized.csv", row.names = T)
write.csv(as.data.frame(t(subset_B_zscore)), "Lowinput_all_lipids_final_non_normalized.csv", row.names = T)

getwd()
