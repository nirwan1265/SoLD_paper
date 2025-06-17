# ──────────────────────────────────────────────────────────────
#  0)  PACKAGES  (only ggplot2 & dplyr beyond what you loaded)
# ──────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)
library(scales)

# ──────────────────────────────────────────────────────────────
#  1)  READ RAW INTENSITIES  (same paths as before)
# ──────────────────────────────────────────────────────────────

ctrl_file  <- "/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Control_all_lipids_final_non_normalized.csv"

lowp_file <- "/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv"


control   <- vroom::vroom(ctrl_file,  show_col_types = FALSE) %>%
  mutate(Condition = "Control")
lowinput  <- vroom::vroom(lowp_file,  show_col_types = FALSE) %>%
  mutate(Condition = "LowInput")

# Remove PC(17:0) from control and lowinput (internal standard) from the columns
control  <- control %>% dplyr::select(-`PC(17:0)`)
lowinput <- lowinput %>% dplyr::select(-`PC(17:0)`)


raw_all   <- bind_rows(control, lowinput)


valid_classes <- c("TG","DG","MG",
                   "PC","PE","PI",
                   "DGDG","MGDG",
                   "SQDG","SM","AEG",
                   "LPC","LPE","PG","PA")
class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

library(vroom)
library(dplyr)
library(tidyr)
library(stringr)

# 1. Read data (you have these already)
# control  <- vroom(ctrl_file,  show_col_types = FALSE)
# lowinput <- vroom(lowp_file, show_col_types = FALSE)

control2  <- control  %>% mutate(Condition = "Control")
lowinput2 <- lowinput %>% mutate(Condition = "LowInput")
combined <- bind_rows(control2, lowinput2)

valid_classes <- c("TG","DG","MG",
                   "PC","PE","PI",
                   "DGDG","MGDG",
                   "SQDG","SM","AEG",
                   "LPC","LPE","PG","PA")
class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

df_long2 <- combined %>%
  pivot_longer(
    cols = -c(Compound_Name, Condition),
    names_to  = "Lipid",
    values_to = "Intensity"
  ) %>%
  rename(Sample = Compound_Name) %>%
  mutate(Class = str_extract(Lipid, class_pat)) %>%
  filter(!is.na(Class))

# 2. Helper rounding
round_to_sum_100 <- function(raw_pct_vec, target_sum = 100) {
  n <- length(raw_pct_vec)
  if (n == 1) {
    return(target_sum)
  }
  floored <- floor(raw_pct_vec)
  rem <- raw_pct_vec - floored
  current_sum <- sum(floored)
  diff <- target_sum - current_sum
  if (diff > 0) {
    idx <- order(rem, decreasing = TRUE)[ seq_len(min(diff, n)) ]
    floored[idx] <- floored[idx] + 1
  } else if (diff < 0) {
    idx <- order(rem, decreasing = FALSE)[ seq_len(min(-diff, n)) ]
    floored[idx] <- floored[idx] - 1
  }
  return(floored)
}

# 3. Compute within‐class composition with threshold lumping
compute_within_pct_table <- function(df_long2, condition_name, threshold_pct = 3) {
  df_long2 %>%
    filter(Condition == condition_name) %>%
    group_by(Sample, Class) %>%
    mutate(class_total = sum(Intensity, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(rel_within = if_else(class_total > 0, Intensity / class_total, 0)) %>%
    group_by(Class, Lipid) %>%
    summarise(mean_prop = mean(rel_within, na.rm = TRUE), .groups = "drop") %>%
    group_by(Class) %>%
    group_modify(~ {
      dfc <- .x %>% mutate(mean_prop = if_else(is.na(mean_prop), 0, mean_prop))
      raw_pct <- dfc$mean_prop * 100
      small_mask <- raw_pct < threshold_pct
      sum_small_prop <- sum(dfc$mean_prop[small_mask], na.rm = TRUE)
      if (all(small_mask)) {
        # everything < threshold → single Other
        df2 <- tibble(Lipid = "Other", mean_prop = sum(dfc$mean_prop))
      } else {
        df_big <- dfc[!small_mask, , drop = FALSE]
        if (sum_small_prop > 0) {
          df2 <- bind_rows(df_big,
                           tibble(Lipid = "Other", mean_prop = sum_small_prop))
        } else {
          df2 <- df_big
        }
      }
      raw_pct2 <- df2$mean_prop * 100
      if (nrow(df2) == 1) {
        df2 <- df2 %>% mutate(mean_pct = 100L)
      } else {
        int_pct <- round_to_sum_100(raw_pct2, target_sum = 100)
        df2 <- df2 %>% mutate(mean_pct = int_pct)
      }
      return(df2)
    }) %>%
    ungroup() %>%
    arrange(Class, desc(mean_pct))
}

table_control  <- compute_within_pct_table(df_long2, "Control",  threshold_pct = 3)
table_lowinput <- compute_within_pct_table(df_long2, "LowInput", threshold_pct = 3)

# Check sums:
table_control  %>% group_by(Class) %>% summarise(sum100 = sum(mean_pct))
table_lowinput %>% group_by(Class) %>% summarise(sum100 = sum(mean_pct))

library(dplyr)

# Rename Lipid -> Species for clarity
control_long <- table_control %>%
  select(Class, Species = Lipid, Percentage = mean_pct)

lowinput_long <- table_lowinput %>%
  select(Class, Species = Lipid, Percentage = mean_pct)


# Combine the two 
combined_long <- bind_cols(
  control_long %>% mutate(Condition = "Control"),
  lowinput_long %>% mutate(Condition = "LowInput")
)



# Rename Percentage columns to distinguish
ctrl <- control_long %>%
  rename(Pct_Control = Percentage)

lwp <- lowinput_long %>%
  rename(Pct_LowInput = Percentage)

# Full join by Class & Species
combined <- full_join(ctrl, lwp, by = c("Class", "Species"))

# Flag common species (present in both) and compute a sorting value
combined <- combined %>%
  mutate(
    common = !is.na(Pct_Control) & !is.na(Pct_LowInput),
    # For ordering: compute max percent across conditions, treating NA as 0 for ordering
    max_pct = pmax(replace_na(Pct_Control, 0), replace_na(Pct_LowInput, 0))
  )

combined <- combined %>%
  arrange(Class, desc(common), desc(max_pct), Species) %>%
  dplyr::select(-common, -max_pct)

# See the combined table:
print(combined, n = Inf)





#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
#### SOrghum Lipidomics Database (SOLD)
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


###############################################################################
## Supplementary Figures – Overlap of Control vs Low‑input lipid lists
## Three panels:
##   S1 – All detected lipids               (Venn)
##   S2 – Traditional classes only          (Venn)
##   S3 – Curated “Class” annotation set    (table + bar *optional*)
###############################################################################

# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------
library(vroom)      # fast CSV import
library(dplyr)      # data wrangling
library(tidyr)      # pivoting
library(ggplot2)    # plotting
library(stringr)    # string manipulation
library(forcats)    # factor reordering
library(tidyverse)
library(ggpubr)      # For stat_compare_means
library(scales) 
library(ggVennDiagram)   # nice, ggplot‑based 2‑set Venn
library(readr) 
library(sf)
library(ggVennDiagram) # For Venn diagrams)
library(viridis)

# ------------------------------------------------------------------------------
# 2. Read in the data
# ------------------------------------------------------------------------------
control  <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Control_all_lipids_final_non_normalized.csv")

lowinput <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv")
colnames(control)


# Remove PC(17:0) from control and lowinput (internal standard) from the columns
control  <- control %>% dplyr::select(-`PC(17:0)`)
lowinput <- lowinput %>% dplyr::select(-`PC(17:0)`)

#### CHANGE THE NAME TO COMMONNAME
# Lipid class
lipid_class_info <- vroom::vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/data/lipid_class.csv", show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Class))


# Create a named vector of replacements: names = original, values = new names
name_map <- lipid_class_info %>%
  filter(!is.na(CommonName)) %>%
  select(Lipids, CommonName) %>%
  deframe()  # turns it into a named vector: Lipids -> CommonName

# Function to rename lipid columns in any dataset
rename_lipid_columns <- function(df, name_map) {
  original_names <- colnames(df)
  
  renamed <- original_names %>%
    # Replace only if a match exists
    sapply(function(x) if (x %in% names(name_map)) name_map[[x]] else x)
  
  # Apply renamed vector back
  colnames(df) <- renamed
  return(df)
}

# Apply to both control and lowinput
control  <- rename_lipid_columns(control,  name_map)
lowinput <- rename_lipid_columns(lowinput, name_map)



# Lipid class
lipid_class_info <- vroom::vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/data/lipid_class.csv", show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Class))
colnames(control)
# Lipid class Traditional lipids
traditional_lipid_classes <- c("Glycerolipid", "Glycerophospholipid", "Glycoglycerolipid",
                               "Sphingolipid", "Sterol", "Betaine lipid", "Fatty acid and derivative")
lipid_class_traditional <- lipid_class_info %>% 
  dplyr::filter(Class %in% traditional_lipid_classes)

# Lipid class Traditional lipids
nontraditional_lipid_class <- c("N-acylethanolamine","Terpenoid","Prenol","Tetrapyrrole","Vitamin")

lipid_class_nontraditional <- lipid_class_info %>% 
  dplyr::filter(Class %in% nontraditional_lipid_class)


# Use this updated lipid_class_info
lipid_class_info <- lipid_class_info %>%
  mutate(final_name = ifelse(is.na(CommonName), Lipids, CommonName))

# Redefine traditional and nontraditional sets
lipid_class_traditional <- lipid_class_info %>%
  filter(Class %in% traditional_lipid_classes)

lipid_class_nontraditional <- lipid_class_info %>%
  filter(Class %in% nontraditional_lipid_class)

# Use final_name column now
traditional_lipids     <- lipid_class_traditional$final_name
nontraditional_lipids  <- lipid_class_nontraditional$final_name

# Assuming you renamed control/lowinput already
trad_cols_control     <- intersect(traditional_lipids, colnames(control))
nontrad_cols_control  <- intersect(nontraditional_lipids, colnames(control))

trad_cols_lowinput    <- intersect(traditional_lipids, colnames(lowinput))
nontrad_cols_lowinput <- intersect(nontraditional_lipids, colnames(lowinput))

# Subset
control_traditional     <- control  %>% select(Compound_Name, all_of(trad_cols_control))
control_nontraditional  <- control  %>% select(Compound_Name, all_of(nontrad_cols_control))
lowinput_traditional    <- lowinput %>% select(Compound_Name, all_of(trad_cols_lowinput))
lowinput_nontraditional <- lowinput %>% select(Compound_Name, all_of(nontrad_cols_lowinput))



# 1) First, make sure you have a lookup from final_name → Class:
lookup <- lipid_class_info %>%
  mutate(final_name = ifelse(is.na(CommonName), Lipids, CommonName)) %>%
  select(final_name, Class)

# 2) Grab the lipid‐column names (minus the sample ID) for each condition/set
ctr_trad_cols   <- setdiff(names(control_traditional),    "Compound_Name")
low_trad_cols   <- setdiff(names(lowinput_traditional),   "Compound_Name")
ctr_non_cols    <- setdiff(names(control_nontraditional), "Compound_Name")
low_non_cols    <- setdiff(names(lowinput_nontraditional),"Compound_Name")

# 3) A helper that, for a vector of classes, builds the summary table
make_summary <- function(classes, ctrl_cols, low_cols){
  lapply(classes, function(cl){
    # which lipids in this class?
    members <- lookup %>% filter(Class == cl) %>% pull(final_name)
    in_ctrl <- intersect(ctrl_cols, members)
    in_low  <- intersect(low_cols,  members)
    tibble(
      Class            = cl,
      Shared           = paste0(intersect(in_ctrl, in_low), collapse = "; "),
      Unique_Control   = paste0(setdiff(in_ctrl, in_low), collapse = "; "),
      Unique_Lowinput  = paste0(setdiff(in_low,  in_ctrl), collapse = "; ")
    )
  }) %>% bind_rows()
}

# 4) Build both tables
trad_summary <- make_summary(
  traditional_lipid_classes,
  ctr_trad_cols, low_trad_cols
)

nontrad_summary <- make_summary(
  nontraditional_lipid_class,
  ctr_non_cols, low_non_cols
)

# 5) Inspect
print(trad_summary)
print(nontrad_summary)

# 6) (Optional) write out
write.csv(trad_summary,    "SuppTable_1_traditional_lipid_overlap.csv",    row.names = FALSE)
write.csv(nontrad_summary, "SuppTable_2_nontraditional_lipid_overlap.csv", row.names = FALSE)
