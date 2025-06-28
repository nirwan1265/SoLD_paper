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

length(unique(rownames(A_filtered)))
length(rownames(A_filtered))
################################################################################
####################### OBTAINING LIPID NAMES  #################################
################################################################################

##### Filtering based on lipids
lipid_A <- read.csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/lipid_class_A.csv")
lipid_A$X.Scan. <- as.character(lipid_A$X.Scan.)
lipid_B <- read.csv("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/data/lipid_class_B.csv")
lipid_B$X.Scan. <- as.character(lipid_B$X.Scan.)

library(dplyr)
library(stringr)

##### ── 1. Build a minimal lookup of scan → compound name ───────────────────
lookup_A <- lipid_A %>% 
  dplyr::select(X.Scan., Compound_Name)

lookup_B <- lipid_B %>% 
  dplyr::select(X.Scan., Compound_Name)

##### ── 2. Join annotations into your intensity tables ───────────────────────
A_with_names <- inner_join(A_filtered, lookup_A, by = "X.Scan.")
B_with_names <- inner_join(B_filtered, lookup_B, by = "X.Scan.")

##### ── 3. Remove exact duplicates of (scan, compound, intensities) ─────────
A_unique <- A_with_names %>% 
  distinct(X.Scan., Compound_Name, across(where(is.numeric)), .keep_all = TRUE)

B_unique <- B_with_names %>% 
  distinct(X.Scan., Compound_Name, across(where(is.numeric)), .keep_all = TRUE)

##### ── 4. Identify which columns are your sample intensities ────────────────
int_cols_A <- setdiff(names(A_unique), c("X.Scan.", "Compound_Name"))
int_cols_B <- setdiff(names(B_unique), c("X.Scan.", "Compound_Name"))

##### ── 5. Collapse to **one row per compound** by summing **unique** values ─
final_A <- A_unique %>%
  group_by(Compound_Name) %>%
  summarise(
    across(all_of(int_cols_A), ~ sum(unique(.x), na.rm = TRUE)),
    .groups = "drop"
  )

final_B <- B_unique %>%
  group_by(Compound_Name) %>%
  summarise(
    across(all_of(int_cols_B), ~ sum(unique(.x), na.rm = TRUE)),
    .groups = "drop"
  )

# final_A / final_B now have:
#  • one row per Compound_Name
#  • each sample column = sum of that compound’s *unique* feature intensities
#  • duplicates with identical intensities are only counted once




library(dplyr)
library(stringr)

# 1. Identify which columns are your sample intensities
int_cols_A <- setdiff(colnames(final_A), "Compound_Name")
int_cols_B <- setdiff(colnames(final_B), "Compound_Name")

# 2. For A: strip CE suffix, then group‑and‑sum unique compounds
final_A_cleaned <- final_A %>%
  mutate(
    Compound_Base = str_remove(Compound_Name, " CollisionEnergy:.*$")
  ) %>%
  group_by(Compound_Base) %>%
  summarise(
    across(all_of(int_cols_A), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(Compound_Name = Compound_Base)

# 3. For B: same procedure
final_B_cleaned <- final_B %>%
  mutate(
    Compound_Base = str_remove(Compound_Name, " CollisionEnergy:.*$")
  ) %>%
  group_by(Compound_Base) %>%
  summarise(
    across(all_of(int_cols_B), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(Compound_Name = Compound_Base)

# final_A_cleaned and final_B_cleaned now:
# • have one row per base compound name (CE stripped)
# • each sample column = sum of intensities from all rows that collapsed into that name





library(dplyr)
library(stringr)

# Suppose final_A_cleaned is your CE‑stripped, compound‑level table
# with columns: Compound_Name, Sample1, Sample2, …

# 1. Remove the Massbank prefix (covers both "Massbank: " and "Massbank:XYZ123 ")
A_massbank_fixed <- final_A_cleaned %>%
  mutate(
    Compound_Name = str_remove(Compound_Name,
                               "^Massbank:[^ ]*\\s*")
  )

# 2. Identify your intensity columns again
int_cols <- setdiff(names(A_massbank_fixed), "Compound_Name")

# 3. Collapse duplicates that now share the same name
final_A_nomassbank <- A_massbank_fixed %>%
  group_by(Compound_Name) %>%
  summarise(
    across(all_of(int_cols), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Repeat for B:
B_massbank_fixed <- final_B_cleaned %>%
  mutate(Compound_Name = str_remove(Compound_Name, "^Massbank:[^ ]*\\s*"))

int_cols_B <- setdiff(names(B_massbank_fixed), "Compound_Name")

final_B_nomassbank <- B_massbank_fixed %>%
  group_by(Compound_Name) %>%
  summarise(
    across(all_of(int_cols_B), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )





library(dplyr)
library(stringr)

# 1. Starting from your Massbank‑fixed tables
#    (or you can apply directly to final_A_cleaned / final_B_cleaned)
A_nomassbank <- final_A_nomassbank
B_nomassbank <- final_B_nomassbank

# 2. Remove the MassbankEU: prefix + any non‑space chars + optional space
A_nomassbankeu <- A_nomassbank %>%
  mutate(
    Compound_Name = str_remove(Compound_Name,
                               "^MassbankEU:[^ ]*\\s*")
  ) %>%
  # 3. Re‑collapse in case this created duplicates
  group_by(Compound_Name) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

B_nomassbankeu <- B_nomassbank %>%
  mutate(
    Compound_Name = str_remove(Compound_Name,
                               "^MassbankEU:[^ ]*\\s*")
  ) %>%
  group_by(Compound_Name) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )




library(dplyr)
library(stringr)

# --- for A ---
A_cleaned <- A_nomassbankeu %>%
  mutate(
    # 1) drop "Spectral Match to " at start
    Compound_Name = str_remove(Compound_Name, "^Spectral Match to\\s*"),
    # 2) drop " from NIST14" (or any " from <WORD>") at end
    Compound_Name = str_remove(Compound_Name, "\\s+from\\s+\\S+$")
  ) %>%
  # 3) re‑collapse in case this merging created duplicates
  group_by(Compound_Name) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# --- for B ---
B_cleaned <- B_nomassbankeu %>%
  mutate(
    Compound_Name = str_remove(Compound_Name, "^Spectral Match to\\s*"),
    Compound_Name = str_remove(Compound_Name, "\\s+from\\s+\\S+$")
  ) %>%
  group_by(Compound_Name) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# A_cleaned / B_cleaned are your final tables:
# • no Massbank: or MassbankEU: prefixes
# • no "Spectral Match to " or " from NIST14" wrappers
# • one row per cleaned Compound_Name
# • intensities summed across any scans that collapsed to the same name




library(dplyr)
library(stringr)

# ─── A: remove MoNA: prefixes ─────────────────────────────────────────────────
A_noMoNA <- A_cleaned %>%  # or A_cleaned from your Massbank/NIST cleanup
  mutate(
    # strip any "MoNA:123456 " prefix
    Compound_Name = str_remove(Compound_Name, "^MoNA:[^ ]*\\s*")
  ) %>%
  # re‑collapse in case identical names now merged
  group_by(Compound_Name) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# ─── B: same for B table ─────────────────────────────────────────────────────
B_noMoNA <- B_cleaned %>% 
  mutate(
    Compound_Name = str_remove(Compound_Name, "^MoNA:[^ ]*\\s*")
  ) %>%
  group_by(Compound_Name) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Now A_noMoNA / B_noMoNA:
# • have no MoNA:##### prefixes
# • one row per unique cleaned Compound_Name
# • intensities summed across any rows that collapsed together





library(dplyr)
library(stringr)

# Starting from your MoNA‑ and Massbank‑cleaned tables:
#   A_noMoNA, B_noMoNA

# 1. Remove everything from the first ";" onward 
A_final <- A_noMoNA %>%
  mutate(
    Compound_Name = str_remove(Compound_Name, ";.*$")
  ) %>%
  # 2. Re‑collapse in case this created new duplicates
  group_by(Compound_Name) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

B_final <- B_noMoNA %>%
  mutate(
    Compound_Name = str_remove(Compound_Name, ";.*$")
  ) %>%
  group_by(Compound_Name) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# A_final / B_final now:
# • Have clean compound names (no adduct/formula suffix)
# • One row per compound, intensities summed across all features





library(dplyr)
library(stringr)

# Starting from your fully cleaned tables:
#   A_final  and  B_final

# 1. For A: drop everything from the first '|' onward
A_pipe_cleaned <- A_final %>%
  mutate(
    Compound_Clean = str_replace(Compound_Name, "\\|.*$", "")
  ) %>%
  group_by(Compound_Clean) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(Compound_Name = Compound_Clean)

# 2. For B: same procedure
B_pipe_cleaned <- B_final %>%
  mutate(
    Compound_Clean = str_replace(Compound_Name, "\\|.*$", "")
  ) %>%
  group_by(Compound_Clean) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(Compound_Name = Compound_Clean)

# A_pipe_cleaned and B_pipe_cleaned now:
#  • Have no '|' or trailing synonyms—only the first name remains
#  • One row per unique Compound_Name
#  • Intensities summed across any features that collapsed into that name




library(dplyr)
library(stringr)

# 1. Remove space before "(" in compound names, then re-collapse
A_final2 <- A_pipe_cleaned %>%
  mutate(Compound_Clean = str_replace(Compound_Name, " \\(", "(")) %>%
  group_by(Compound_Clean) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(Compound_Name = Compound_Clean)

B_final2 <- B_pipe_cleaned %>%
  mutate(Compound_Clean = str_replace(Compound_Name, " \\(", "(")) %>%
  group_by(Compound_Clean) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(Compound_Name = Compound_Clean)

# Now PC-DAG (16:0/18:3) becomes PC-DAG(16:0/18:3), etc.




library(dplyr)
library(stringr)

# 1. For A:
A_norm <- A_final2 %>%
  mutate(Compound_Clean = str_replace_all(Compound_Name,
                                          "\\bDAG\\b",   # whole word “DAG”
                                          "DG")) %>%
  group_by(Compound_Clean) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(Compound_Name = Compound_Clean)

# 2. For B:
B_norm <- B_final2 %>%
  mutate(Compound_Clean = str_replace_all(Compound_Name,
                                          "\\bDAG\\b",
                                          "DG")) %>%
  group_by(Compound_Clean) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename(Compound_Name = Compound_Clean)

# A_norm and B_norm now have “DAG” → “DG” everywhere, and any duplicates are merged.






### Name change:
library(dplyr)

# 1. Read your name-change table
name_change <- read.csv(
  "/Users/nirwantandukar/Documents/Research/results/SAP/Lipid Class/lipid_class.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
glimpse(name_change)
# 2. Define a helper to apply to any df with a Compound_Name column
apply_common_names <- function(df) {
  df %>%
    # join on the exact Compound_Name
    left_join(name_change, by = "Compound_Name") %>%
    # whenever CommonName is non‑empty, use it; otherwise keep the original
    mutate(
      Compound_Name = coalesce(
        na_if(CommonName, ""),  # treat "" as NA
        Compound_Name
      )
    ) %>%
    select(-CommonName)  # drop the lookup column
}

# 3. Apply to your two tables
A_final_named <- apply_common_names(A_norm)
B_final_named <- apply_common_names(B_norm)

# Inspect
head(A_final_named$Compound_Name, 20)




library(dplyr)
library(rlang)

# 1) Read your lookup table (with Class, SubClass, Sub-subclass included)
name_change <- read.csv(
  "/Users/nirwantandukar/Documents/Research/results/SAP/Lipid Class/lipid_class.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 2) Enhanced apply_common_names:
apply_common_names <- function(df) {
  # perform a left_join, pulling in CommonName + taxonomy columns
  joined <- df %>%
    left_join(name_change, by = "Compound_Name", suffix = c("", "_new")) %>%
    
    # replace with CommonName if non-empty
    mutate(
      Compound_Name = coalesce(na_if(CommonName, ""), Compound_Name)
    )
  
  # list of taxonomy columns to reconcile
  tax_cols <- c("Class", "SubClass", "Sub-subclass")
  
  # for each taxonomy column, merge the original and _new versions
  for (col in tax_cols) {
    orig <- col
    new  <- paste0(col, "_new")
    
    if (all(c(orig, new) %in% names(joined))) {
      joined <- joined %>%
        mutate(
          !!orig := case_when(
            is.na(.data[[new]])             ~ .data[[orig]],       # only original
            .data[[orig]] == .data[[new]]  ~ .data[[orig]],       # identical
            TRUE                            ~ paste(.data[[orig]],
                                                    .data[[new]],
                                                    sep = ";")         # differ → combine
          )
        )
    }
  }
  
  # drop helper columns
  joined %>%
    dplyr::select(-CommonName, -ends_with("_new"))
}

# 3) Apply to your two final tables
A_final_named <- apply_common_names(A_norm)
B_final_named <- apply_common_names(B_norm)

# 4) Bring those taxonomy columns to just after Compound_Name
A_final_named <- A_final_named %>%
  relocate(Class, SubClass, `Sub-subclass`, .after = Compound_Name)

B_final_named <- B_final_named %>%
  relocate(Class, SubClass, `Sub-subclass`, .after = Compound_Name)





library(dplyr)

A_clean_names <- A_final_named %>%
  rename(Sub_subclass = `Sub-subclass`)  # use the exact dash you see in colnames()

A_clean_names <- A_clean_names %>%
  mutate(Compound_Clean = str_replace_all(Compound_Name, "\\bTAG\\b", "TG")) %>%
  group_by(Compound_Clean, Class, SubClass, Sub_subclass) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups="drop") %>%
  rename(Compound_Name = Compound_Clean)


B_clean_names <- B_final_named %>%
  rename(Sub_subclass = `Sub-subclass`)  # use the exact dash you see in colnames()


B_clean_names <- B_clean_names %>%
  mutate(Compound_Clean = str_replace_all(Compound_Name, "\\bTAG\\b", "TG")) %>%
  group_by(Compound_Clean, Class, SubClass, Sub_subclass) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups="drop") %>%
  rename(Compound_Name = Compound_Clean)



# Change Lyso-PC to LPC
A_clean_names <- A_clean_names %>%
  mutate(Compound_Name = str_replace(Compound_Name, "\\bLyso-PC\\b", "LPC"))
B_clean_names <- B_clean_names %>%
  mutate(Compound_Name = str_replace(Compound_Name, "\\bLyso-PC\\b", "LPC"))






library(dplyr)
library(stringr)
library(tibble)

# 1) Define the headgroup → taxonomy map (as before)
tax_map <- tribble(
  ~Headgroup, ~Class,               ~SubClass,                             ~Sub_subclass,
  "PC",       "Glycerophospholipid","Glycerophosphocholine",               "Phosphatidylcholine",
  "LPC",      "Glycerophospholipid","Glycerophosphocholine",               "Lysophosphatidylcholine",
  "PE",       "Glycerophospholipid","Glycerophosphoethanolamine",           "Phosphatidylethanolamine",
  "LPE",      "Glycerophospholipid","Glycerophosphoethanolamine",           "Lysophosphatidylethanolamine",
  "PA",       "Glycerophospholipid","Glycerophosphate",                     "Phosphatidic acid",
  "PG",       "Glycerophospholipid","Glycerophosphoglycerol",               "Phosphatidylglycerol",
  "PS",       "Glycerophospholipid","Glycerophosphoserine",                 "Phosphatidylserine",
  "DG",       "Glycerolipid",      "Diacylglycerol",                     "DG",
  "TG",       "Glycerolipid",      "Triacylglycerol",                     "TG",
  "MG",       "Glycerolipid",      "Monoacylglycerols",                    "MG",
  "DGDG",     "Glycerolipid",      "Glycosyldiacylglycerol",              "Digalactosyldiacylglycerol",
  "MGDG",     "Glycerolipid",      "Glycosyldiacylglycerol",              "Monogalactosyldiacylglycerol",
  "SQDG",     "Glycerolipid",      "Sulfoquinovosyldiacylglycerol",       "SQDG",
  "DGTS",     "Glycerophospholipid","Glycerophosphocholine",               "DGTS",
  "Cer",      "Sphingolipid",      "Ceramides",                            "Ceramides",
  "SM",       "Sphingolipid",      "Sphingomyelin",                       "Sphingomyelins",
  "AEG",      "Ether lipid",       "Acyl‑ether glycerol",                 "AEG"
)

# 2) Helper to grab headgroup (prefix before any "(")
extract_head <- function(name) {
  str_extract(name, "^[A-Za-z0-9]+")
}

# 3) Fill in taxonomy, preserving any existing values
fill_taxonomy <- function(df) {
  df %>%
    mutate(Headgroup = extract_head(Compound_Name)) %>%
    left_join(
      tax_map,
      by      = "Headgroup",
      suffix  = c("", "_map")   # original columns stay as-is, map columns get _map
    ) %>%
    mutate(
      Class        = coalesce(Class,        Class_map),
      SubClass     = coalesce(SubClass,     SubClass_map),
      Sub_subclass = coalesce(Sub_subclass, Sub_subclass_map)
    ) %>%
    select(-Headgroup, -ends_with("_map"))
}

# 4) Apply to your final cleaned tables
A_filled <- fill_taxonomy(A_clean_names)
B_filled <- fill_taxonomy(B_clean_names)

# 5) If you like, move the taxonomy right after Compound_Name
A_filled <- A_filled %>%
  relocate(Class, SubClass, Sub_subclass, .after = Compound_Name)

B_filled <- B_filled %>%
  relocate(Class, SubClass, Sub_subclass, .after = Compound_Name)

# Now A_filled and B_filled have:
# • Existing Class/SubClass/Sub_subclass kept when present
# • Missing ones filled from tax_map



library(dplyr)

# 1️⃣ Define exactly the classes you want to drop
classes_to_drop <- c(
  "Plasticizer",
  "Psychoactive compound",
  "Functional agent",
  "Functional compound",
  "Herbicide",
  "Lipid binding protein",
  "Organic compound",
  "Phenol derivative",
  "Phenylpropanoid",
  "Surfactant",
  "Thioester"
)

# 2️⃣ Pull out those “junk” compounds for inspection
A_junk <- A_filled %>% 
  filter(Class %in% classes_to_drop)

B_junk <- B_filled %>% 
  filter(Class %in% classes_to_drop)

# 3️⃣ Your cleaned lipid tables (everything NOT in those classes)
A_clean <- A_filled %>% 
  filter(!Class %in% classes_to_drop)

B_clean <- B_filled %>% 
  filter(!Class %in% classes_to_drop)

# 4️⃣ Quick check of what remains
table(A_junk$Class)   # see the junk pulled out
table(A_clean$Class)  # see your lipid classes left

# Now A_junk/B_junk hold all the unwanted classes for manual review,
# and A_clean/B_clean are your lipid‑only tables ready for analysis.

# Save the cleaned data
#write.csv(A_clean, 
#          "data/raw_lipid_intensities/A_cleaned_lipids.csv", 
#          row.names = FALSE)
#write.csv(B_clean,
#          "data/raw_lipid_intensities/B_cleaned_lipids.csv", 
#          row.names = FALSE)


# NOTE: GO TO EXCEL AND CHECK WITH THE LIPID CLASSES AND SUM THE C AND H BONDS

# Read the cleaned data to sum the common lipids now:
A_clean <- vroom("data/raw_lipid_intensities/A_cleaned_lipids.csv")
B_clean <- vroom("data/raw_lipid_intensities/B_cleaned_lipids.csv")




library(dplyr)

A_unique <- A_clean %>%
  # Step 1: Compute rowwise average (across all sample columns)
  rowwise() %>%
  dplyr::mutate(row_avg = mean(c_across(matches("^S\\d")), na.rm = TRUE)) %>%
  ungroup() %>%
  
  # Step 2: Group by the average and collapse names
  dplyr::group_by(row_avg) %>%
  summarise(
    Compound_Name = paste(sort(unique(Compound_Name)), collapse = "::"),
    across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
    across(c(Presummed_name, Class, SubClass, Sub_subclass), ~first(na.omit(.x)))
  ) %>%
  dplyr::select(-row_avg)  # Optional: remove the helper column

# Done
A_unique



B_unique <- B_clean %>%
  rowwise() %>%
  dplyr::mutate(row_avg = mean(c_across(matches("^S\\d")), na.rm = TRUE)) %>%
  ungroup() %>%
  
  dplyr::group_by(row_avg) %>%
  summarise(
    Compound_Name = paste(sort(unique(Compound_Name)), collapse = "::"),
    across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
    across(c(Presummed_name, Class, SubClass, Sub_subclass), ~first(na.omit(.x)))
  ) %>%
  dplyr::select(-row_avg)


B_unique



# Sum the lipids. 


A_summed <- A_unique %>%
  dplyr::group_by(Compound_Name) %>%
  dplyr::summarise(
    across(where(is.numeric), sum, na.rm = TRUE),
    Class = first(Class),
    SubClass = first(SubClass),
    Sub_subclass = first(Sub_subclass),
    .groups = "drop"
  )

colnames(A_summed)
B_summed <- B_unique %>%
  dplyr::group_by(Compound_Name) %>%
  dplyr::summarise(
    across(where(is.numeric), sum, na.rm = TRUE),
    Class = first(Class),
    SubClass = first(SubClass),
    Sub_subclass = first(Sub_subclass),
    .groups = "drop"
  )




# Clean sample names in column names
colnames(A_summed) <- colnames(A_summed) %>%
  stringr::str_remove("^S1_") %>%                 # remove leading S1_
  stringr::str_remove("_Run\\d+$")                # remove trailing _Run<number>


colnames(B_summed) <- colnames(B_summed) %>%
  stringr::str_remove("^S1_") %>%                 # remove leading S1_
  stringr::str_remove("_Run\\d+$")                # remove trailing _Run<number>


# Save the cleaned and summed data
#write.csv(A_summed, "data/summed_lipid_intensities/A_summed_lipids.csv", row.names = FALSE)
#write.csv(B_summed, "data/summed_lipid_intensities/B_summed_lipids.csv", row.names = FALSE)

# Check the multi-names and select the unique names. 
# Remove shady lipids and mark them as delete and delete them
# Also remove the Odd number carbons atoms lipids


# Sum them again:
A_summed <- vroom("data/summed_lipid_intensities/A_summed_lipids.csv")
A_summed_final <- A_summed %>%
  # make sure there’s no stray whitespace in the names
  mutate(Compound_Name = str_trim(Compound_Name)) %>%
  # now collapse any rows that share the same name, summing all their intensities
  group_by(Compound_Name) %>%
  summarise(
    across(where(is.numeric), sum, na.rm = TRUE),
    .groups = "drop"
  )

str(A_summed)


B_summed <- vroom("data/summed_lipid_intensities/B_summed_lipids.csv")
B_summed_final <- B_summed %>%
  group_by(Compound_Name) %>%
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Save 
write.csv(A_summed_final, "data/summed_lipid_intensities/A_summed_lipids_final.csv", row.names = FALSE)
# write.csv(B_summed_final, "data/summed_lipid_intensities/B_summed_lipids_final.csv", row.names = FALSE)
