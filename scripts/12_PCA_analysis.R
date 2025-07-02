
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###############################################################################
## INDIVIDUAL LIPIDS
###############################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



###############################################################################
## 1)  LOAD PACKAGES
###############################################################################

library(vroom)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(broom)      
library(patchwork)  
library(scales)
library(cowplot)
library(ggh4x)


###############################################################################
## 2) READ RAW INTENSITY TABLES AND CLASS             
###############################################################################

# Read the raw files 
control <- vroom(
  "data/SPATS_fitted/non_normalized_intensities/control_all_lipids_fitted_phenotype_non_normalized.csv",
  show_col_types = FALSE,
  col_types = cols(.default = "c")
) %>%
  dplyr::select(-2:-4) %>%                # drop the 3 metadata cols you don’t need
  mutate(Condition = "Control") %>%
  rename(Compound_Name = 1)

lowinput <- vroom(
  "data/SPATS_fitted/non_normalized_intensities/lowinput_all_lipids_fitted_phenotype_non_normalized.csv",
  show_col_types = FALSE,
  col_types = cols(.default = "c")
) %>%
  dplyr::select(-2:-4) %>%
  mutate(Condition = "LowInput") %>%
  rename(Compound_Name = 1)


###############################################################################
## 3)  SUBSET TRADITIONAL LIPID CLASSES
###############################################################################

# Valid classes
valid_classes <- c("TG","DG","MG",
                   "PC","PE","PI",
                   "DGDG","MGDG",
                   "SQDG","SM","AEG",
                   "LPC","LPE","PG","PA","PS")

valid_classes <- c("TG","DG",
                   "PC","PE",
                   "DGDG","MGDG",
                   "SQDG"
                   )


# build a regex like "^(TG|DG|MG|…)"
class_pattern <- paste0("^(", paste(valid_classes, collapse = "|"), ")")

# find lipid columns in control
lipid_cols_ctrl <- grep(class_pattern, names(control), value = TRUE)

# subset control to samples + Condition + those lipids
control_sub <- control %>%
  dplyr::select(Compound_Name, Condition, all_of(lipid_cols_ctrl))

# Convert everything to numeric starting from 3rd column
control_sub <- control_sub %>%
  mutate(across(all_of(lipid_cols_ctrl), as.numeric))


# do the same for lowinput
lipid_cols_low <- grep(class_pattern, names(lowinput), value = TRUE)
lowinput_sub <- lowinput %>%
  dplyr::select(Compound_Name, Condition, all_of(lipid_cols_low))
lowinput_sub <- lowinput_sub %>%
  mutate(across(all_of(lipid_cols_low), as.numeric))



###############################################################################
## 4)  SELECT ONLY SHARED LIPIDS
###############################################################################

# 1. Use the intersected lipid list (shared lipids)
shared_lipids <- intersect(lipid_cols_ctrl, lipid_cols_low)

# 2. Subset and bind together
combined_df <- bind_rows(
  control %>% dplyr::select(Compound_Name, Condition, all_of(shared_lipids)),
  lowinput %>% dplyr::select(Compound_Name, Condition, all_of(shared_lipids))
)

# 3. Convert lipid columns to numeric
combined_df <- combined_df %>%
  mutate(across(all_of(shared_lipids), as.numeric))

# Remove rows with string count less than 5 from Compound_Name
combined_df <- combined_df %>%
  filter(nchar(Compound_Name) >= 5)

# 4. Drop samples with any NA (if any lipid is missing)
pca_input <- combined_df %>%
  drop_na()


# 1) Create the pca_input with a unique Sample_ID
pca_input <- combined_df %>%
  drop_na() %>%                                            # drop any rows with missing lipids
  mutate(Sample_ID = paste(Compound_Name, Condition, sep = "_"))

# 2) Build a pure-numeric matrix with those rownames
lipid_mat <- pca_input %>%
  dplyr::select(Sample_ID, all_of(shared_lipids)) %>%            # keep Sample_ID + all lipids
  as.data.frame() %>%                                     # convert tibble → data.frame
  column_to_rownames("Sample_ID")                         # turn Sample_ID into rownames

# 3) (Optionally) scale & center
#lipid_mat_scaled <- scale(lipid_mat, center = TRUE, scale = TRUE)

# 4) Run PCA
pca_res <- prcomp(lipid_mat, center = TRUE, scale. = TRUE)

# 5) Extract scores and re‐attach metadata
scores_df <- as.data.frame(pca_res$x[,1:2]) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(
    pca_input %>% select(Compound_Name, Condition, Sample_ID),
    by = "Sample_ID"
  )

# 3) Extract scores and merge in metadata
scores_df <- data.frame(pca_res$x[,1:2], Sample_ID = rownames(pca_res$x))
scores_df <- merge(
  scores_df,
  pca_input %>% dplyr::select(Sample_ID, Compound_Name, Condition),
  by = "Sample_ID",
  all.x = TRUE
)

# 4) Prepare loadings as arrows

load_df <- data.frame(pca_res$rotation[,1:2], Class = rownames(pca_res$rotation))
scale_factor <- 5
# compute a data‐driven scale factor
sf <- max(abs(scores_df$PC1), abs(scores_df$PC2)) /
  max(sqrt(load_df$PC1^2 + load_df$PC2^2))
# now multiply loadings by, say, 0.8*sf
load_df <- load_df %>% mutate(
  a1 = PC1 * (0.8 * sf),
  a2 = PC2 * (0.8 * sf)
)
load_df <- load_df %>%
  mutate(a1 = PC1 * scale_factor,
         a2 = PC2 * scale_factor)

# 5) Plot biplot
quartz()
ggplot() +
  # Samples + ellipses
  geom_point(data = scores_df,
             aes(PC1, PC2, color = Condition),
             size = 3, alpha = 0.8) +
  stat_ellipse(data = scores_df,
               aes(PC1, PC2, fill = Condition),
               geom = "polygon", alpha = 0.2, color = NA) +
  # Loading arrows
  geom_segment(data = load_df,
               aes(x = 0, y = 0, xend = a1, yend = a2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "navy", size = 0.8) +
  geom_text(data = load_df,
            aes(x = a1 * 1.05, y = a2 * 1.05, label = Class),
            color = "navy", size = 3) +
  # Axes, labels, theme
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    title = "Biplot: Samples & Lipid-Class Loadings",
    x = paste0("PC1 (", round(100 * pca_res$sdev[1]^2 / sum(pca_res$sdev^2),1), "%)"),
    y = paste0("PC2 (", round(100 * pca_res$sdev[2]^2 / sum(pca_res$sdev^2),1), "%)")
  ) +
  scale_color_manual(values = c(Control = "#F8766D", LowInput = "#00BFC4")) +
  scale_fill_manual(values  = c(Control = "#F8766D", LowInput = "#00BFC4")) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "right",
        plot.title     = element_text(face = "bold"))




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###############################################################################
## SUMMED LIPIDS
###############################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# ╔══════════════════════════════════════════════════════════════════╗
# ║ 1) READ RAW INTENSITY TABLES AND CLASS                           ║
# ╚══════════════════════════════════════════════════════════════════╝

# Read the raw files 
control  <- vroom("data/SPATS_fitted/non_normalized_intensities/control_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
colnames(control)[1] <- "Compound_Name"  

lowinput  <- vroom("data/SPATS_fitted/non_normalized_intensities/lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
colnames(lowinput)[1] <- "Compound_Name"  

# Valid classes
valid_classes <- c("TG","DG","MG",
                   "PC","PE","PI",
                   "DGDG","MGDG",
                   "SQDG","SM","AEG",
                   "LPC","LPE","PG","PA","PS")

class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

raw_filt <- raw_all %>%
  filter(str_detect(Compound_Name, class_pat))

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 2) TIC normalisation + log10 per lipid                           ║
# ╚══════════════════════════════════════════════════════════════════╝

# Long format
long_prep_global <- function(df) {
  df_long <- df %>%
    pivot_longer(-c(Compound_Name, Condition),
                 names_to = "Lipid", values_to = "Intensity") %>%
    rename(Sample = Compound_Name)
  
  # 1) TIC per sample → relative abundance
  df_long <- df_long %>%
    group_by(Sample) %>%
    mutate(TIC       = sum(Intensity, na.rm = TRUE),
           rel_abund = Intensity / TIC) %>%
    ungroup()
  
  # 2) small pseudo-count per sample → log10
  df_long <- df_long %>%
    group_by(Sample) %>%
    mutate(minpos = ifelse(
      all(is.na(rel_abund) | rel_abund <= 0),
      NA_real_, min(rel_abund[rel_abund > 0], na.rm = TRUE)),
      eps    = ifelse(is.na(minpos), 0, minpos * 0.5)) %>%
    ungroup() %>%
    mutate(log_rel = log10(rel_abund + eps))
  
  # 3) summarise mean log10(relative abundance) per lipid class
  df_long <- df_long %>%
    mutate(Class = str_extract(Lipid, class_pat)) %>%
    filter(!is.na(Class)) %>%
    group_by(Sample, Condition, Class) %>%
    summarise(class_log = mean(log_rel, na.rm = TRUE), .groups = "drop")
  
  return(df_long)
}

# combine Control + LowInput, then run prep
combined  <- bind_rows(control  %>% mutate(Condition = "Control"),
                       lowinput %>% mutate(Condition = "LowInput"))
long_all  <- long_prep_global(combined)

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 3)  PIVOT WIDE                                                   ║
# ╚══════════════════════════════════════════════════════════════════╝

# Wide format
wide_log <- long_all %>%
  pivot_wider(names_from = Class, values_from = class_log,
              values_fill = NA_real_)

# Removing  AEG and PG 
#wide_log <- wide_log[,-c(3,13)]
wide_log <- wide_log[,-c(3,8)]



# 1) pivot into “long” so we have one row per Sample/Condition/Class
lc <- wide_log %>%
  pivot_longer(
    cols      = -c(Sample, Condition),
    names_to  = "Class",
    values_to = "Abundance"
  )


# 2) Pivot wide so each row = one Sample×Condition, each col = one Class
wide <- lc %>%
  pivot_wider(
    names_from  = Class,
    values_from = Abundance
  )

# 1) Prepare sample info
sample_info <- wide %>%
  dplyr::select(Sample, Condition) %>%
  mutate(Sample_ID = paste(Sample, Condition, sep = "_"))

# 2) Prepare matrix for PCA
mat <- wide %>%
  dplyr::select(-Sample, -Condition)

mat <- as.data.frame(wide %>% dplyr::select(-Sample, -Condition))
rownames(mat) <- sample_info$Sample_ID

# 3) Run PCA (on scaled data)
p <- prcomp(mat, center = TRUE, scale. = TRUE)

# 4) PCA scores (Sample positions)
scores_df <- as.data.frame(p$x[, 1:2]) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(sample_info, by = "Sample_ID")  # attach Condition

# 5) Loadings (Lipid-class arrows)
load_df <- as.data.frame(p$rotation[, 1:2]) %>%
  rownames_to_column("Class") %>%
  mutate(a1 = PC1 * 5,  # scale for visibility
         a2 = PC2 * 5)

# 6) Plot
quartz()
ggplot() +
  # sample points + ellipses
  geom_point(data = scores_df,
             aes(PC1, PC2, color = Condition),
             size = 2.5, alpha = 0.7) +
  stat_ellipse(data = scores_df,
               aes(PC1, PC2, fill = Condition),
               geom = "polygon", alpha = 0.2, color = NA) +
  # loading arrows
  geom_segment(data = load_df,
               aes(x = 0, y = 0, xend = a1, yend = a2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "navy", size = 0.8) +
  # loading labels
  geom_text(data = load_df,
            aes(x = a1 * 1.05, y = a2 * 1.05, label = Class),
            color = "navy", size = 3) +
  # axis lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  # labels
  labs(
    title = "Biplot: Samples & Lipid-Class Loadings",
    x = paste0("PC1 (", round(100 * p$sdev[1]^2 / sum(p$sdev^2), 1), "%)"),
    y = paste0("PC2 (", round(100 * p$sdev[2]^2 / sum(p$sdev^2), 1), "%)")
  ) +
  # colors
  scale_color_manual(values = c(Control = "#F8766D", LowInput = "#00BFC4")) +
  scale_fill_manual(values  = c(Control = "#F8766D", LowInput = "#00BFC4")) +
  coord_fixed() +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###############################################################################
## RATIOS OF THE LIPIDS
###############################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Ratios
# 2) self-join on Sample+Condition, then keep only one direction of each pair
ratios_long <- lc %>%
  inner_join(long_classes, 
             by     = c("Sample","Condition"),
             suffix = c(".num", ".den")) %>%
  # only keep each unordered pair once
  filter(Class.num < Class.den) %>%
  transmute(
    Sample,
    Condition,
    RatioName = paste(Class.num, Class.den, sep = "/"),
    Ratio     = Abundance.num / Abundance.den
  )

# 3) inspect    
ratios_long


# Remove rows with string count less than 3 from Sample
ratios_long <- ratios_long %>% 
  filter(nchar(Sample) >= 5)


unique(ratios_long$RatioName)


library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# 1) pivot ratios_long into a wide matrix (rows = each sample×condition, cols = each RatioName)
wide_ratios <- ratios_long %>%
  # make a unique row‐name combining Sample + Condition
  mutate(Rep = paste(Sample, Condition, sep = "_")) %>%
  dplyr::select(Rep, RatioName, Ratio) %>%
  pivot_wider(names_from  = RatioName,
              values_from = Ratio) %>%
  column_to_rownames("Rep")


# 1) build metadata: each row = one Sample×Condition
sample_info_ratios <- ratios_long %>% 
  distinct(Sample, Condition) %>% 
  mutate(Rep = paste(Sample, Condition, sep = "_"))




# 2) run PCA (autoscaled)
pca_ratios <- prcomp(wide_ratios, center = TRUE, scale. = TRUE)


# 3) extract scores and re‐attach metadata
scores_df <- as.data.frame(pca_ratios$x[,1:2], check.names = FALSE) %>%
  setNames(c("PC1","PC2")) %>%
  tibble::rownames_to_column("Rep") %>%
  inner_join(sample_info_ratios, by = "Rep")

# 4) extract loadings and compute a data‐driven arrow scale
load_df <- as.data.frame(pca_ratios$rotation[,1:2], check.names = FALSE) %>%
  setNames(c("PC1","PC2")) %>%
  tibble::rownames_to_column("RatioName")

# compute scale so that 80% of the longest arrow matches the score‐cloud radius
max_score <- max(abs(c(scores_df$PC1, scores_df$PC2)))
max_load  <- max(sqrt(load_df$PC1^2 + load_df$PC2^2))
arrow_scale <- 0.8 * max_score / max_load

load_df <- load_df %>%
  mutate(a1 = PC1 * arrow_scale,
         a2 = PC2 * arrow_scale)

# 5) biplot
quartz()
ggplot() +
  # samples + ellipses
  geom_point(data = scores_df,
             aes(PC1, PC2, color = Condition),
             size = 2.5, alpha = 0.7) +
  stat_ellipse(data = scores_df,
               aes(PC1, PC2, fill = Condition),
               geom = "polygon", alpha = 0.2, color = NA) +
  # loading arrows
  geom_segment(data = load_df,
               aes(x = 0, y = 0, xend = a1, yend = a2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "navy", size = 0.8) +
  # loading labels
  geom_text(data = load_df,
            aes(x = a1 * 1.05, y = a2 * 1.05, label = RatioName),
            color = "navy", size = 3) +
  # axis lines & fixed aspect
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  coord_fixed() +
  # axis labels with % variance
  labs(
    title = "Biplot: Samples & Lipid Ratios Loadings",
    x = paste0("PC1 (", round(100 * pca_ratios$sdev[1]^2 /
                                sum(pca_ratios$sdev^2), 1), "%)"),
    y = paste0("PC2 (", round(100 * pca_ratios$sdev[2]^2 /
                                sum(pca_ratios$sdev^2), 1), "%)")
  ) +
  scale_color_manual(values = c(Control = "#F8766D",
                                LowInput = "#00BFC4")) +
  scale_fill_manual(values = c(Control = "#F8766D",
                               LowInput = "#00BFC4")) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title     = element_text(face = "bold")
  )


