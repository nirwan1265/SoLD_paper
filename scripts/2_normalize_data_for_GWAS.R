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
#### DATA WRANGLE
################################################################################
#-------------------------------------------------------------------------------


###############################################################################
# REQUIRED LIBRARIES
###############################################################################
suppressPackageStartupMessages({
  library(vroom);   library(dplyr);   library(tidyr);  library(stringr)
  library(ggplot2); library(viridis); library(gridExtra); library(grid)
  library(tibble); library(patchwork); library(scales)
})


###############################################################################
# LOAD THE RAW READS 
###############################################################################

control <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Control_all_lipids_final_non_normalized.csv")
lowinput <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv")


###############################################################################
# LOAD THE LIPID CLASSES
###############################################################################

lipid_class_info <- vroom::vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/data/lipid_class.csv", show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Class))
lipid_class_info

unique(lipid_class_info$Class)


###############################################################################
# MAP THE LIPIDS TO THE CLASSES
###############################################################################

# Get lipid names from your data (i.e. column names except for the first column 'Compound_Name')
lipid_names_lowinput <- colnames(lowinput)[-1]
lipid_names_control <- colnames(control)[-1]

# Match to the lipid_class_info
matched_class_info_lowinput <- lipid_class_info %>%
  filter(Lipids %in% lipid_names_lowinput)
matched_class_info_control <- lipid_class_info %>%
  filter(Lipids %in% lipid_names_control)

# Group by Class and concatenate lipid names
lipids_by_class_lowinput <- matched_class_info_lowinput %>%
  group_by(Class) %>%
  summarise(
    Lipids = paste(Lipids, collapse = ", "),
    Count = n()
  ) %>%
  ungroup()

lipids_by_class_control <- matched_class_info_control %>%
  group_by(Class) %>%
  summarise(
    Lipids = paste(Lipids, collapse = ", "),
    Count = n()
  ) %>%
  ungroup()

# Change the Lipid names to Common Name
lipids_by_class_lowinput <- lipids_by_class_lowinput %>%
  # 1. Turn each comma‐separated “Lipids” string into one row per lipid
  separate_rows(Lipids, sep = ",\\s*") %>%
  # 2. Join in “CommonName” from lipid_class_info
  left_join(
    lipid_class_info %>% select(Lipids, CommonName),
    by = c("Lipids" = "Lipids")
  ) %>%
  # 3. If CommonName is non‐NA/non‐blank, use it; otherwise keep the original Lipids
  mutate(
    Lipid_display = if_else(
      !is.na(CommonName) & CommonName != "",
      CommonName,
      Lipids
    )
  ) %>%
  # 4. Re‐group by Class and collapse back into comma‐separated string
  group_by(Class) %>%
  summarise(
    Lipids         = paste(Lipid_display, collapse = ", "),
    Count          = n(),              # recomputed count of lipids
    .groups        = "drop"
  )

# Do the same for Control:
lipids_by_class_control <- lipids_by_class_control %>%
  separate_rows(Lipids, sep = ",\\s*") %>%
  left_join(
    lipid_class_info %>% select(Lipids, CommonName),
    by = c("Lipids" = "Lipids")
  ) %>%
  mutate(
    Lipid_display = if_else(
      !is.na(CommonName) & CommonName != "",
      CommonName,
      Lipids
    )
  ) %>%
  group_by(Class) %>%
  summarise(
    Lipids         = paste(Lipid_display, collapse = ", "),
    Count          = n(),
    .groups        = "drop"
  )


# View result
print(lipids_by_class_lowinput, n = 30)
print(lipids_by_class_control, n = 30)


### Remove Herbicide, Plasticizer, and Surfactant from control and lowinput
# ── 1.  Classes you want gone ─────────────────────────────────────
rm_classes <- c("Herbicide", "Plasticizer", "Surfactant")

# helper: turn the comma-separated column into a clean character vector
split_names <- function(tbl) {
  tbl %>%                                           # lipids_by_class_…
    filter(Class %in% rm_classes) %>%               # keep target classes
    pull(Lipids) %>%                                # character strings
    str_split(",\\s*") %>%                          # split on comma + spaces
    unlist() %>%                                    # flatten
    str_trim() %>%                                  # remove leading/trailing spaces
    unique()
}

lipids_to_remove_control  <- split_names(lipids_by_class_control)
lipids_to_remove_lowinput <- split_names(lipids_by_class_lowinput)

# finally drop them (any_of = safe)
control  <- control  %>% select(-any_of(lipids_to_remove_control))
lowinput <- lowinput %>% select(-any_of(lipids_to_remove_lowinput))



## Z SCORE CALCULATION

# 1. Extract the numeric matrix (exclude Compound_Name column)
lipid_matrix_control <- control %>%
  column_to_rownames("Compound_Name") %>%  # optional: set sample ID as rowname
  as.matrix()

lipid_matrix_lowinput <- lowinput %>%
  column_to_rownames("Compound_Name") %>%  # optional: set sample ID as rowname
  as.matrix()

# 2. Calculate Z-scores per column (lipid-wise)
zscore_matrix_control <- scale(lipid_matrix_control, center = TRUE, scale = TRUE)
zscore_matrix_lowinput <- scale(lipid_matrix_lowinput, center = TRUE, scale = TRUE)

# 3. Turn it back into a tibble with sample names
zscore_df_control <- as_tibble(zscore_matrix_control, rownames = "Compound_Name")
zscore_df_lowinput <- as_tibble(zscore_matrix_lowinput, rownames = "Compound_Name")

# Save 
#write_csv(zscore_df_control, "SAP_zscore_control.csv")
#write_csv(zscore_df_lowinput, "SAP_zscore_lowinput.csv")


## LOG₁₀ TRANSFORMATION  ─────────────────────────────────────────────

# Choose a small pseudocount to prevent log10(0); adjust if needed
pseudo <- 0.5    # or e.g. 0.5 * min(non-zero intensity)

# 1. Extract numeric matrices (as you did for Z-scores)
lipid_matrix_control  <- control  %>% column_to_rownames("Compound_Name") %>% as.matrix()
lipid_matrix_lowinput <- lowinput %>% column_to_rownames("Compound_Name") %>% as.matrix()

# 2. Apply log10(x + pseudo) element-wise
log10_matrix_control  <- log10(lipid_matrix_control  + pseudo)
log10_matrix_lowinput <- log10(lipid_matrix_lowinput + pseudo)

# 3. Convert back to tibbles with sample IDs
log10_df_control  <- as_tibble(log10_matrix_control,  rownames = "Compound_Name")
log10_df_lowinput <- as_tibble(log10_matrix_lowinput, rownames = "Compound_Name")

# 4. Inspect
print(log10_df_control)
print(log10_df_lowinput)

# save
# write_csv(log10_df_control, "SAP_log10_control.csv")
# write_csv(log10_df_lowinput, "SAP_log10_lowinput.csv")


# # ── Distribution of “gibberellic acid” ────────────────────────────
# 
# # This script builds one tibble that holds the **raw**, **log10(+c)**, 
# # and **Z-score** values for the compound across both conditions, then 
# # draws faceted histograms (one panel per transformation).
# 
# 
# target  <- "gibberellic acid"   # exact column label
# pseudo  <- 0.5                    # pseudocount for log10
# 
# # 1. Fetch raw intensities from each data frame --------------------
# if (!target %in% colnames(control)  ||
#     !target %in% colnames(lowinput)) {
#   stop("Compound not found in both tables — check the column name.")
# }
# 
# raw_df <- bind_rows(
#   control  %>% select(Compound_Name, all_of(target)) %>%
#     mutate(Condition = "Control"),
#   lowinput %>% select(Compound_Name, all_of(target)) %>%
#     mutate(Condition = "LowInput")
# ) %>% 
#   rename(Raw = !!target)                         # “Raw” column
# 
# # 2. Add log10(+pseudo) --------------------------
# raw_df <- raw_df %>% 
#   mutate(Log10 = log10(Raw + pseudo))
# 
# # 3. Add per-condition Z-scores ------------------
# z_control  <- scale(control[[target]],  center = TRUE, scale = TRUE)
# z_lowinput <- scale(lowinput[[target]], center = TRUE, scale = TRUE)
# 
# raw_df <- raw_df %>% 
#   mutate(Z = c(as.vector(z_control), as.vector(z_lowinput)))
# 
# # 4. Long format for ggplot ----------------------
# plot_df <- raw_df %>% 
#   pivot_longer(c(Raw, Log10, Z), names_to = "Metric", values_to = "Value")
# 
# # 5. Draw faceted histograms ---------------------
# quartz()
# ggplot(plot_df, aes(Value, fill = Condition)) +
#   geom_histogram(bins = 30, alpha = 0.55, position = "identity") +
#   facet_wrap(~ Metric, scales = "free") +
#   scale_fill_manual(values = c(Control = "#4C72B0",
#                                LowInput = "#DD8452")) +
#   labs(title = target,
#        y      = "Count",
#        x      = NULL,
#        fill   = NULL) +
#   theme_minimal(base_size = 12)
# 
# 
# 
