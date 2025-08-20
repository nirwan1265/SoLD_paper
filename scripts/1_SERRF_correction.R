library(stringr)
################################################################################
############ Loading the Quantitative Peaks Intensities.  ######################
################################################################################

# Control
A <- read.csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/SetA_lipid_FLO2019Control.csv")

# Retention time > 1 min
A <- A[which(A$row.retention.time >= 1), ]

# Separating the tables
A_peaks <- A[14:882]
A_id <- A[1]

# 1. grab your original names
orig <- colnames(A_peaks)

# 2. pull out the feature ID: either PI<number> or CHECK<1|2>_<repeat>
feat <- str_extract(orig, "PI\\d+|CHECK[12]_\\d+")

# 3. pull out the run number
run  <- str_extract(orig, "(?<=Run)\\d+")

# 4. build the new names by type
new_names <- ifelse(
  str_detect(orig, "\\.S1_Run"), 
  # for real samples:
  paste0("S1_", feat, "_Run", run),
  
  ifelse(
    str_detect(orig, "InjBL"), 
    # for blanks:
    paste0("InjBL_Run", run),
    
    ifelse(
      str_detect(orig, "ISTD"), 
      # for internal standards:
      paste0("ISTD_Run", run),
      
      ifelse(
        str_detect(orig, "QC"), 
        # for QC samples:
        paste0("QC_Run", run),
        NA_character_
      )
    )
  )
)

# 5. assign back
colnames(A_peaks) <- new_names

# 6. (optional) reorder left-to-right by run number
run2 <- as.integer(str_extract(new_names, "(?<=Run)\\d+"))
ord  <- order(run2, na.last = TRUE)
A_peaks <- A_peaks[, ord]

# Final Table
A <- cbind(A_id,A_peaks)
names(A)
str(A)



################################################################################
########################## Quality Control  ####################################
################################################################################

# Step 1: Calculate the average of all the InjBL values in each row
A$InjBL_average <- rowMeans(A[, grep("^InjBL", names(A))])

# Step 2 & 3: Count the number of S1 values greater than 10 times the average of InjBL
A$S1_gt_10x_InjBL <- rowSums(A[, grep("^S1_", names(A))] > 10 * A$InjBL_average)

# Step 4: Check if 70% or more of the S1 values meet the criterion
threshold <- 0.7 * ncol(A[, grep("^S1_", names(A))])  # 70% threshold

# Filter rows based on the criterion
cleaned_data <- A[A$S1_gt_10x_InjBL >= threshold, ]

# Remove auxiliary columns
cleaned_data <- cleaned_data[, !(names(cleaned_data) %in% c("InjBL_average", "S1_gt_10x_InjBL"))]

A <- cleaned_data


#Removing columns with InjBL cause SERRF does not want it
A <- A %>%
  dplyr::select(-contains("InjBL"))

### Removing the CHECKS. all the columns with S1_Run*_NA are CHECKS
# Use grep to find column names that match the pattern "S1_Run" followed by any number and "_NA"
# columns_to_remove <- grep("S1_Run\\d+_NA", names(A), value = TRUE)
# 
# # Remove these columns from the data frame
# A <- A[, !(names(A) %in% columns_to_remove)]


names(A)


### Replace zero with 2/3 of the lowest value in a row. 
# Find the smallest non-zero value for each column
smallest_nonzero <- apply(A[, -1], 2, function(x) min(x[x > 0]))

# Replace zeros with 2⁄3 of smallest non-zero value for each column
for (col in names(A)[-1]) {
  smallest_nonzero <- min(A[A[, col] > 0, col])
  A[A[, col] == 0, col] <- 2/3 * smallest_nonzero
}
# View the first few rows of the cleaned data
head(A)



################################################################################
############################# For SERRF  #######################################
################################################################################

# Create a new row with labels based on column names
new_row <- data.frame("A")

# Add QC, Sample, or Validate labels based on column names
col_names <- colnames(A)
for (col_name in col_names) {
  if (grepl("QC_", col_name)) {
    new_label <- "qc"
  } else if (grepl("S1_", col_name)) {
    new_label <- "sample"
  } else if (grepl("ISTD_", col_name)) {
    new_label <- "validate"
  } else if (grepl("row.ID", col_name)) {
    new_label <- "label"
  } else {
    new_label <- ""
  }
  
  new_row[[col_name]] <- new_label
  
}
new_row <- as.data.frame(new_row[,-1])

# Add the new row to the data frame
A <- rbind(new_row, A)

# Adding sample numbers
# Get the number of columns in A_peaks
num_cols <- ncol(A) - 1

# Create a new row with increasing numbers from 1 to num_cols
new_row <- c("time", 1:num_cols)

# Add the new row as the first row in A_peaks
A <- rbind(new_row, A)

# Count the column names that start with "QC"
qc_count <- sum(grepl("^QC", colnames(A)))
qc_count

# Count columns with names starting with "S1_Run"
sum(grepl("^S1_PI", names(A)))
sum(grepl("^S1_CHECK", names(A)))

# Adding batch
batch <- read.csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/batch_run.csv")
batch_A <- batch[which(batch$Set == "A"), 1:4]

head(batch_A)

# Subsetting the numbers from Run
# Extract the numbers from the column names and create a new row
extracted_numbers <- str_extract(names(A), "Run(\\d+)")

# Remove "Run" and convert to numeric
extracted_numbers <- as.numeric(gsub("Run", "", extracted_numbers))

# Initialize a new vector to store the batch values
batch_values <- character(length(extracted_numbers))

# Iterate through each value in extracted_numbers
for (i in 2:length(extracted_numbers)) {
  run_value <- as.numeric(extracted_numbers[i])
  
  # Find the corresponding batch value based on the range
  batch_value <- batch_A$Batch[which(run_value >= batch_A$Run_start & run_value <= batch_A$Run_stop)]
  
  # If a batch value is found, assign it; otherwise, use NA
  if (length(batch_value) > 0) {
    batch_values[i] <- batch_value
  } else {
    batch_values[i] <- NA
  }
}

# Replace the first element with "batch"
batch_values[1] <- "batch"

# SERRF needs atleast 6 batches for 1 batch. so adding up the batches:
# Calculate the counts of each letter except "batch"
letter_counts <- table(batch_values)
letter_counts <- letter_counts[!names(letter_counts) %in% "batch"]
letter_counts

# A and B are 1 batch and G H I J are another 
# Replace A and B with A
batch_values[batch_values %in% c("A", "B","C","D","E","F","G","H","I","J","K")] <- "A"

# Replace G, H, I, J with G
batch_values[batch_values %in% c("L","M","N","O","P","Q","R")] <- "B"

# Replace G, H, I, J with G
batch_values[batch_values %in% c("S","T","U","V")] <- "C"


# Find unique non-"batch" letters in the batch_values vector
unique_letters <- unique(batch_values[batch_values != "batch"])

# Create a mapping from unique letters to serial letters (e.g., A, B, C, ...)
serial_mapping <- setNames(LETTERS[seq_along(unique_letters)], unique_letters)

# Update the batch_values vector using the mapping (excluding "batch")
batch_values[batch_values != "batch"] <- serial_mapping[batch_values[batch_values != "batch"]]


# Add the extracted numbers as a new row to your data frame
A <- rbind(batch_values, A)

# Saving
write.csv(A,"/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/results/SERRF/A_SERRF_new.csv")

# Need some post processing for the SERRF
#batch	A	A	A	A	A	A	A	A	A	A
#sampleType	qc	validate	sample	sample	sample	sample	sample	sample	sample	sample
#time	1	2	3	4	5	6	7	8	9	10
#No	label	QC000	sample01	GB001617	GB001333	GB001191	GB001827	GB001722	GB001468	GB001543	GB001347
#1	1_ISTD Ceramide (d18:1/17:0) [M+HCOO]- 	167879	185671	158256	164492	155000	150957	134195	184272	165878	157758

# Website for running SERRF:
#https://slfan.shinyapps.io/ShinySERRF/







################################################################################
############ Loading the Quantitative Peaks Intensities.  ######################
################################################################################
## FOR LOWINPUT



library(stringr)
################################################################################
############ Loading the Quantitative Peaks Intensities.  ######################
################################################################################

B <- read.csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/SetB_lipid_FLO2022_lowP.csv")
B <- B[,-791]

# Retention time > 1 min
B <- B[which(B$row.retention.time >= 1), ]

# Separating the tables
B_peaks <- B[14:790]
B_id <- B[1]
names(B_peaks)

x <- as.data.frame(colnames(B_peaks))


# 1. grab the original column names
orig <- colnames(B_peaks)

# 2. extract the “feature” part: either PI<number>_<rep> or Check_<rep>
#feat <- str_extract(orig, "PI\\d+_\\d+|[Cc]heck_\\d+")
feat <- str_extract(orig, "PI\\d+|[Cc]heck_\\d+")
# 3. extract the Run number
run  <- str_extract(orig, "(?<=Run)\\d+")
run_int <- as.integer(run)

# 4. build the new names by injection type
new_names <- ifelse(
  str_detect(orig, "S1_Run"), 
  paste0("S1_",    feat, "_Run", run),      # real samples
  ifelse(
    str_detect(orig, "InjBL"),  
    paste0("InjBL_Run",    run),            # injection blanks
    ifelse(
      str_detect(orig, "ISTD"),  
      paste0("ISTD_Run",     run),          # internal standards (if any)
      ifelse(
        str_detect(orig, "QC_Run"),  
        paste0("QC_Run",       run),        # QC samples
        NA_character_
      )
    )
  )
)

# 5. reassign and reorder by Run
colnames(B_peaks) <- new_names
ord <- order(run_int, na.last = TRUE)

B_peaks <- B_peaks[, ord]
colnames(B_peaks)

# Final Table
B <- cbind(B_id,B_peaks)
names(B)
str(B)


################################################################################
########################## Quality Control  ####################################
################################################################################

# Step 1: Calculate the average of all the InjBL values in each row
B$InjBL_average <- rowMeans(B[, grep("^InjBL", names(B))])

# Step 2 & 3: Count the number of S1 values greater than 10 times the average of InjBL
B$S1_gt_10x_InjBL <- rowSums(B[, grep("^S1_", names(B))] > 10 * B$InjBL_average)

# Step 4: Check if 70% or more of the S1 values meet the criterion
threshold <- 0.7 * ncol(B[, grep("^S1_", names(B))])  # 70% threshold

# Filter rows based on the criterion
cleaned_data <- B[B$S1_gt_10x_InjBL >= threshold, ]

# Remove auxiliary columns
cleaned_data <- cleaned_data[, !(names(cleaned_data) %in% c("InjBL_average", "S1_gt_10x_InjBL"))]

B <- cleaned_data

#Removing columns with InjBL cause SERRF does not want it
B <- B %>%
  dplyr::select(-contains("InjBL"))

### Removing the CHECKS. all the columns with S1_Run*_NA are CHECKS
# Use grep to find column names that match the pattern "S1_Run" followed by any number and "_NA"
columns_to_remove <- grep("S1_Run\\d+_NA", names(B), value = TRUE)

# Remove these columns from the data frame
B <- B[, !(names(B) %in% columns_to_remove)]

names(B)

### Replace zero with 2/3 of the lowest value in a row. 
# Find the smallest non-zero value for each column
smallest_nonzero <- apply(B[, -1], 2, function(x) min(x[x > 0]))

# Replace zeros with 2⁄3 of smallest non-zero value for each column
for (col in names(B)[-1]) {
  smallest_nonzero <- min(B[B[, col] > 0, col])
  B[B[, col] == 0, col] <- 2/3 * smallest_nonzero
}
# View the first few rows of the cleaned data
head(B)



################################################################################
############################# For SERRF  #######################################
################################################################################

# Create B new row with labels based on column names
new_row <- data.frame("B")


# Add QC, Sample, or Validate labels based on column names
col_names <- colnames(B)
for (col_name in col_names) {
  if (grepl("QC_", col_name)) {
    new_label <- "qc"
  } else if (grepl("S1_", col_name)) {
    new_label <- "sample"
  } else if (grepl("ISTD_", col_name)) {
    new_label <- "validate"
  } else if (grepl("row.ID", col_name)) {
    new_label <- "label"
  } else {
    new_label <- ""
  }
  
  new_row[[col_name]] <- new_label
  
}
new_row <- as.data.frame(new_row[,-1])

# Add the new row to the data frame
B <- rbind(new_row, B)

# Adding sample numbers
# Get the number of columns in B_peaks
num_cols <- ncol(B) - 1

# Create B new row with increasing numbers from 1 to num_cols
new_row <- c("time", 1:num_cols)

# Add the new row as the first row in B_peaks
B <- rbind(new_row, B)

# Count the column names that start with "QC"
qc_count <- sum(grepl("^QC", colnames(B)))

# Count columns with names starting with "S1_Run"
sum(grepl("^S1_Run", names(B)))


# Adding batch
batch <- read.csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/batch_run.csv")
batch_B <- batch[which(batch$Set == "B"), 1:4]

head(batch_B)


# Subsetting the numbers from Run
# Extract the numbers from the column names and create a new row
extracted_numbers <- str_extract(names(B), "Run(\\d+)")

# Remove "Run" and convert to numeric
extracted_numbers <- as.numeric(gsub("Run", "", extracted_numbers))

# Initialize a new vector to store the batch values
batch_values <- character(length(extracted_numbers))

# Iterate through each value in extracted_numbers
for (i in 2:length(extracted_numbers)) {
  run_value <- as.numeric(extracted_numbers[i])
  
  # Find the corresponding batch value based on the range
  batch_value <- batch_B$Batch[which(run_value >= batch_B$Run_start & run_value <= batch_B$Run_stop)]
  
  # If a batch value is found, assign it; otherwise, use NA
  if (length(batch_value) > 0) {
    batch_values[i] <- batch_value
  } else {
    batch_values[i] <- NA
  }
}


# Replace the first element with "batch"
batch_values[1] <- "batch"


# SERRF needs atleast 6 batches for 1 batch. so adding up the batches:
# Calculate the counts of each letter except "batch"
letter_counts <- table(batch_values)
letter_counts <- letter_counts[!names(letter_counts) %in% "batch"]
letter_counts

# A and B are 1 batch and G H I J are another 
# Replace A and B with A
batch_values[batch_values %in% c("A", "B","C","D","E","F")] <- "A"

# Replace A and B with A
batch_values[batch_values %in% c( "G","H","I","J","L","K", "M")] <- "B"

# Replace A and B with A
#batch_values[batch_values %in% c("G", "H","I")] <- "C"

# Replace L and M with L
#batch_values[batch_values %in% c("I","J","L","K", "M")] <- "C"

# Find unique non-"batch" letters in the batch_values vector
unique_letters <- unique(batch_values[batch_values != "batch"])

# Create a mapping from unique letters to serial letters (e.g., A, B, C, ...)
serial_mapping <- setNames(LETTERS[seq_along(unique_letters)], unique_letters)

# Update the batch_values vector using the mapping (excluding "batch")
batch_values[batch_values != "batch"] <- serial_mapping[batch_values[batch_values != "batch"]]



# Add the extracted numbers as a new row to your data frame
B <- rbind(batch_values, B)


# Saving
write.csv(B,"/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/results/SERRF/B_SERRF_new.csv")

# Need some post processing for the SERRF
#batch	B	B	B	B	B	B	B	B	B	B
#sampleType	qc	validate	sample	sample	sample	sample	sample	sample	sample	sample
#time	1	2	3	4	5	6	7	8	9	10
#No	label	QC000	sample01	GB001617	GB001333	GB001191	GB001827	GB001722	GB001468	GB001543	GB001347
#1	1_ISTD Ceramide (d18:1/17:0) [M+HCOO]- 	167879	185671	158256	164492	155000	150957	134195	184272	165878	157758

# Website for running SERRF:
#https://slfan.shinyapps.io/ShinySERRF/

x <- as.vector(unlist(B[3,]))
table(x)






###### GWAS phenotype
# Log10 and cenetered