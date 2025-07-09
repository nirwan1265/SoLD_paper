# ───────────────────────────────────────────────────────────────────────────────
# 0.  Load packages
# ───────────────────────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(stringr)
library(vroom)
library(SpATS)
library(gstat)

# ───────────────────────────────────────────────────────────────────────────────
# 1.  Read & tidy the field map and phenotypes
# ───────────────────────────────────────────────────────────────────────────────
fieldmap_lowinput <- read.csv("data/fieldmap/fieldmap_lowinput.csv", header = FALSE)
nrow(fieldmap_lowinput)
ncol(fieldmap_lowinput)

# phenotype <- as.tibble(t(vroom::vroom(
#   "/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv"
# )))

# 1) Read & transpose ----------------------------------------------------------------
raw <- read.csv(
  "data/summed_lipid_intensities/B_summed_lipids_final_max.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)

# transpose to get lines as rows, compounds as columns
pheno_mat <- as.data.frame(t(raw), stringsAsFactors = FALSE)
# first row of that (the original header) becomes column names
colnames(pheno_mat) <- pheno_mat[1, ]
pheno_mat <- pheno_mat[-1, ]

# these are the exact row-names in your transposed object
annot_rows <- c("Class", "SubClass", "Sub_subclass")

# subset pheno_mat to just those three rows
class_raw <- pheno_mat[annot_rows, , drop = FALSE]

dup_names <- colnames(class_raw)[duplicated(colnames(class_raw))]
dup_names


class_df_lowinput <- class_raw %>% 
  # bring the rownames ("Class", etc) into a column
  rownames_to_column(var = "Annotation") %>% 
  # pivot long: one row per Annotation × Compound
  pivot_longer(
    cols      = -Annotation,
    names_to  = "Compound",
    values_to = "Value"
  ) %>% 
  # pivot back so each Compound has its 3 columns
  pivot_wider(
    names_from  = Annotation,
    values_from = Value
  )

colnames(class_df_lowinput)[1] <- "Lipids"
#write.csv(class_df, "data/lipid_class/lipid_classes_control.csv", row.names = FALSE)

# 3) isolate the *numeric* portion (the actual lines × compounds) --------------------
data_rows <- setdiff(rownames(pheno_mat), annot_rows)
pheno_num <- pheno_mat[data_rows, , drop = FALSE]
# convert those to numeric
pheno_num[] <- lapply(pheno_num, as.numeric)

# add back your line IDs
pheno_num <- pheno_num %>%
  as_tibble(rownames = "LineRaw") %>%
  mutate(
    LineID = str_extract(LineRaw, "\\d+")
  ) %>%
  filter(!is.na(LineID))

# pivot the 10×12 map into long format
map_long <- fieldmap_lowinput %>%
  mutate(row = row_number()) %>%
  pivot_longer(
    cols      = starts_with("V"),
    names_to  = "col",
    values_to = "LineRaw"
  ) %>%
  mutate(
    col    = as.integer(str_remove(col, "^V")),
    LineID = case_when(
      LineRaw == "X"        ~ NA_character_,
      TRUE                  ~ str_remove(LineRaw, "^PI_")
    )
  ) %>%
  filter(!is.na(LineID))

# rename & clean the phenotype table
pheno <- pheno_num

# pheno <- phenotype %>%
#   rename(LineID = Compound_Name) %>%
#   mutate(LineID = str_extract(LineID, "\\d+"))
# tail(pheno)

# Assuming filtered_data has columns: LineRaw, LineID, row, col, PlotID, plus all your lipid names
meta_cols <- c("LineRaw", "LineID", "row", "col", "PlotID")

# Make sure you don’t accidentally include any of those in your trait list:
# trait_names <- setdiff(names(filtered_data), meta_cols)
# print(trait_names)


# ───────────────────────────────────────────────────────────────────────────────
# 2.  Prepare loop over all lipid traits
# ───────────────────────────────────────────────────────────────────────────────
trait_names <- setdiff(names(pheno), c("LineRaw","LineID"))

# merge map + pheno
filtered_data <- pheno %>%
  inner_join(map_long, by = "LineID") %>%
  mutate(
    row = as.integer(row),
    col = as.integer(col)
  )


# start a BLUP‐accumulator with one row per LineID
blup_df <- data.frame(LineID = pheno$LineID, stringsAsFactors = FALSE)

fitted_df <- filtered_data %>%
  # keep identifiers for each plot
  dplyr::select(LineID, row, col) %>%
  dplyr::mutate(PlotID = row_number()) %>%
  dplyr::select(PlotID, everything())  # PlotID, LineID, row, col

# ───────────────────────────────────────────────────────────────────────────────
# 3.  Loop: fit SpATS, pull BLUPs, join into blup_df; 
#    also save diagnostics for DG(20:1) only
# ───────────────────────────────────────────────────────────────────────────────
for(tr in trait_names) {
  message("Fitting trait: ", tr)
  
  m <- SpATS(
    response            = tr,
    spatial             = ~ SAP(col, row, nseg = c(8,2), degree = 3, pord = 2), # use 10% x 2 of rows and columns
    genotype            = "LineID",
    data                = filtered_data,
    control             = list(tolerance = 1e-3, maxit = 500),
    genotype.as.random  = TRUE
  )
  
  # extract the BLUP vector
  geno_names  <- m$terms$geno$geno_names
  geno_blups  <- m$coeff[geno_names]
  blups_i     <- data.frame(LineID = geno_names,
                            BLUP   = as.numeric(geno_blups),
                            stringsAsFactors = FALSE)
  
  # join into the master table
  #col_i       <- paste0("BLUP_", tr)
  col_i       <- tr
  blup_df     <- left_join(
    blup_df,
    blups_i %>% rename(!!col_i := BLUP),
    by = "LineID"
  )
  
  
  # ── fitted plot‐values ───────────────────────────────────────────────────────
  fit <- m$fitted  # numeric vector, one per row of filtered_data
  fitted_df <- fitted_df %>%
    mutate(!!tr := fit)
  
  
  # if this is DG(20:1), save the diagnostic plot + variogram
  if (tr == "DG(20:1)") {
    safe <- gsub("[^[:alnum:]]", "_", tr)
    
    png(paste0(safe, "_SpATS_diagnostics.png"), width = 12, height = 12, res = 300, units = "in", bg = "white")
    plot(m)
    dev.off()
    
    png(paste0(safe, "_variogram.png"),width = 12, height = 12, res = 300, units = "in", bg = "white")
    plot(SpATS::variogram(m))
    dev.off()
  }
  
  # png(paste0(safe, "_variogram.png"),
  #     width = 12, height = 12, res = 300, units = "in", bg = "white")
  # var.m0 <- SpATS::variogram(m)
  # plot(var.m0)
  #dev.off()
}

# # variogram of residuals
varobj <- variogram(residuals(m) ~ 1,
                    locations = ~ col + row,
                    data = filtered_data)
png(paste0("model_variogram_residual.png"),
    width = 6, height = 6, res = 300, units = "in", bg = "white")
plot(varobj)
dev.off()


# 1. Add LineRaw to blup_df
blup_df2 <- blup_df %>%
  dplyr::mutate(LineRaw = paste0("PI", LineID)) %>%
  dplyr::select(LineRaw, everything()) %>%
  dplyr::select(-LineID)

# 2. Build a little data.frame for the fitted phenotypes
fitted_df <- fitted_df %>%
  dplyr::mutate(LineRaw = paste0("PI", LineID)) %>%
  dplyr::select(LineRaw, everything()) %>%
  dplyr::select(-LineID)

# 3. Write them out
write.csv(blup_df2,
          "lowinput_all_lipids_BLUPs_max.csv",
          row.names = FALSE)

write.csv(fitted_df,
          "lowinput_all_lipids_fitted_phenotype_non_normalized_max.csv",
          row.names = FALSE)

message("Done!  All BLUPs in all_lipids_BLUPs.csv; DG(20:1) plots saved.")




# ───────────────────────────────────────────────────────────────────────────────
# ───────────────────────────────────────────────────────────────────────────────
# FOR CONTROL
# ───────────────────────────────────────────────────────────────────────────────
# ───────────────────────────────────────────────────────────────────────────────

# ───────────────────────────────────────────────────────────────────────────────
# 1.  Read & tidy the field map and phenotypes
# ───────────────────────────────────────────────────────────────────────────────
fieldmap_control <- read.csv("data/fieldmap/fieldmap_control.csv", header = FALSE)
nrow(fieldmap_control)
ncol(fieldmap_control)

# phenotype <- as.tibble(t(vroom::vroom(
#   "/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/control_all_lipids_final_non_normalized.csv"
# )))

# 1) Read & transpose ----------------------------------------------------------------
raw <- read.csv(
  "data/summed_lipid_intensities/A_summed_lipids_final.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)

# transpose to get lines as rows, compounds as columns
pheno_mat <- as.data.frame(t(raw), stringsAsFactors = FALSE)
# first row of that (the original header) becomes column names
colnames(pheno_mat) <- pheno_mat[1, ]
pheno_mat <- pheno_mat[-1, ]

# these are the exact row-names in your transposed object
annot_rows <- c("Class", "SubClass", "Sub_subclass")

# # subset pheno_mat to just those three rows
# class_raw <- pheno_mat[annot_rows, , drop = FALSE]
# 
# dup_names <- colnames(class_raw)[duplicated(colnames(class_raw))]
# dup_names


class_df_control <- class_raw %>% 
  # bring the rownames ("Class", etc) into a column
  rownames_to_column(var = "Annotation") %>% 
  # pivot long: one row per Annotation × Compound
  pivot_longer(
    cols      = -Annotation,
    names_to  = "Compound",
    values_to = "Value"
  ) %>% 
  # pivot back so each Compound has its 3 columns
  pivot_wider(
    names_from  = Annotation,
    values_from = Value
  )

colnames(class_df_control)[1] <- "Lipids"

# bind lipid class
class_df <- bind_rows(
  class_df_lowinput %>% mutate(Source = "Lowinput"),
  class_df_control %>% mutate(Source = "Control")
)
#write.csv(class_df, "data/lipid_class/lipid_classes_control.csv", row.names = FALSE)

# 3) isolate the *numeric* portion (the actual lines × compounds) --------------------
data_rows <- setdiff(rownames(pheno_mat), annot_rows)
pheno_num <- pheno_mat[data_rows, , drop = FALSE]
# convert those to numeric
pheno_num[] <- lapply(pheno_num, as.numeric)

# add back your line IDs
pheno_num <- pheno_num %>%
  as_tibble(rownames = "LineRaw") %>%
  mutate(
    LineID = str_extract(LineRaw, "\\d+")
  ) %>%
  filter(!is.na(LineID))

# pivot the 10×12 map into long format
map_long <- fieldmap_control %>%
  mutate(row = row_number()) %>%
  pivot_longer(
    cols      = starts_with("V"),
    names_to  = "col",
    values_to = "LineRaw"
  ) %>%
  mutate(
    col = as.integer(str_remove(col, "^V")),
    # Preserve the full CHECK1_#, CHECK2_#, or PI##### label:
    Genotype = case_when(
      LineRaw == "X"                 ~ NA_character_,
      TRUE                           ~ LineRaw
    )
  ) %>%
  filter(!is.na(Genotype))


# rename & clean the phenotype table
pheno <- pheno_num

# pheno <- phenotype %>%
#   rename(LineID = Compound_Name) %>%
#   mutate(LineID = str_extract(LineID, "\\d+"))
# tail(pheno)

# ───────────────────────────────────────────────────────────────────────────────
# 2.  Prepare loop over all lipid traits
# ───────────────────────────────────────────────────────────────────────────────
trait_names <- setdiff(names(pheno), c("LineRaw","LineID"))

# merge map + pheno
filtered_data <- pheno %>%
  inner_join(map_long, by = "LineRaw") %>%
  mutate(
    row      = as.integer(row),
    col      = as.integer(col),
    Genotype = factor(Genotype)
  )


# start a BLUP‐accumulator with one row per LineID
#blup_df <- data.frame(LineID = pheno$LineID, stringsAsFactors = FALSE)
blup_df <- data.frame(Genotype = filtered_data$Genotype, stringsAsFactors = FALSE)

fitted_df <- filtered_data %>%
  # keep identifiers for each plot
  dplyr::select(LineID, row, col) %>%
  dplyr::mutate(PlotID = row_number()) %>%
  dplyr::select(PlotID, everything())  # PlotID, LineID, row, col

# ───────────────────────────────────────────────────────────────────────────────
# 3.  Loop: fit SpATS, pull BLUPs, join into blup_df; 
#    also save diagnostics for DG(20:1) only
# ───────────────────────────────────────────────────────────────────────────────
for(tr in trait_names) {
  message("Fitting trait: ", tr)
  
  m <- SpATS(
    response            = tr,
    spatial             = ~ SAP(col, row, nseg = c(4,10), degree = 3, pord = 2), # use 10% x 2 of rows and columns
    genotype            = "Genotype",
    data                = filtered_data,
    control             = list(tolerance = 1e-3, maxit = 500),
    genotype.as.random  = TRUE
  )
  
  # extract the BLUP vector
  geno_names  <- m$terms$geno$geno_names
  geno_blups  <- m$coeff[geno_names]
  blups_i     <- data.frame(Genotype = geno_names,
                            BLUP   = as.numeric(geno_blups),
                            stringsAsFactors = FALSE)
  
  # join into the master table
  #col_i       <- paste0("BLUP_", tr)
  col_i       <- tr
  blup_df     <- left_join(
    blup_df,
    blups_i %>% rename(!!tr := BLUP),
    by = "Genotype"
  )
  
  
  # ── fitted plot‐values ───────────────────────────────────────────────────────
  fit <- m$fitted  # numeric vector, one per row of filtered_data
  fitted_df <- fitted_df %>%
    mutate(!!tr := fit)
  
  
  # if this is DG(20:1), save the diagnostic plot + variogram
  if (tr == "DG(20:1)") {
    safe <- gsub("[^[:alnum:]]", "_", tr)
    
    png(paste0(safe, "_SpATS_diagnostics.png"), width = 12, height = 12, res = 300, units = "in", bg = "white")
    plot(m)
    dev.off()
    
    png(paste0(safe, "_variogram.png"),width = 12, height = 12, res = 300, units = "in", bg = "white")
    plot(SpATS::variogram(m))
    dev.off()
  }
  
  # png(paste0(safe, "_variogram.png"),
  #     width = 12, height = 12, res = 300, units = "in", bg = "white")
  # var.m0 <- SpATS::variogram(m)
  # plot(var.m0)
  #dev.off()
}

# # variogram of residuals
varobj <- variogram(residuals(m) ~ 1,
                    locations = ~ col + row,
                    data = filtered_data)
png(paste0("control_model_variogram_residual.png"),
    width = 6, height = 6, res = 300, units = "in", bg = "white")
plot(varobj)
dev.off()


# 1. Add LineRaw to blup_df
blup_df2 <- blup_df %>%
  dplyr::mutate(LineRaw = paste0("PI", Genotype)) %>%
  dplyr::select(LineRaw, everything()) %>%
  dplyr::select(-Genotype)

# 2. Build a little data.frame for the fitted phenotypes
fitted_df <- fitted_df %>%
  dplyr::mutate(LineRaw = paste0("PI", LineID)) %>%
  dplyr::select(LineRaw, everything()) %>%
  dplyr::select(-LineID)

# 3. Write them out
write.csv(blup_df2,
          "control_all_lipids_BLUPs.csv",
          row.names = FALSE)

write.csv(fitted_df,
          "control_all_lipids_fitted_phenotype_non_normalized.csv",
          row.names = FALSE)

message("Done!  All BLUPs in all_lipids_BLUPs.csv; DG(20:1) plots saved.")
