# ╔══════════════════════════════════════════════════════════════════╗
# ║ 0)  PACKAGES                                                    ║
# ╚══════════════════════════════════════════════════════════════════╝
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

# keep the Sample and Condition columns around for plotting
sample_info <- wide %>% dplyr::select(Sample, Condition)

# drop them to get only the numeric matrix for PCA
mat <- wide %>% dplyr::select(-Sample, -Condition)

# 3) Run PCA (autoscale on classes)
p <- prcomp(scale(mat), center = TRUE, scale. = TRUE)

# 4) Extract PC scores and re-attach Sample + Condition
scores <- as.data.frame(p$x) %>%
  bind_cols(sample_info)

# 5a) Plot samples colored by Condition
quartz()
ggplot(scores, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Condition), alpha = 0.2, geom = "polygon") +
  labs(
    title = "PCA of Samples (by class-level abundance)",
    x = "PC1", y = "PC2"
  ) +
  theme_bw()

# 5b) If you want to see how the *classes* themselves cluster,
#     just transpose the matrix and repeat:
p2 <- prcomp(scale(t(mat)), center = TRUE, scale. = TRUE)
cls_scores <- as.data.frame(p2$x)
cls_scores$Class <- rownames(cls_scores)

quartz()
ggplot(cls_scores, aes(PC1, PC2, label = Class)) +
  geom_text(size = 3) +
  labs(
    title = "PCA of Lipid Classes (across all samples)",
    x = "PC1", y = "PC2"
  ) +
  theme_bw()





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

# 2) run PCA (autoscaled)
pca_ratios <- prcomp(wide_ratios, center = TRUE, scale. = TRUE)

# 3) extract scores + metadata
scores_ratios <- as.data.frame(pca_ratios$x) %>%
  rownames_to_column("Rep") %>%
  separate(Rep, into = c("Sample", "Condition"), sep = "_")

# 4) plot PC1 vs PC2
quartz()
ggplot(scores_ratios, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(alpha = 0.6, size = 2) +
  stat_ellipse(level = 0.68) +
  labs(title = "PCA of all pairwise lipid ratios",
       x = paste0("PC1 (", round(summary(pca_ratios)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca_ratios)$importance[2,2]*100,1), "%)")) +
  theme_minimal(base_size = 12)




library(ggrepel)

# 1) pull out the PC1/PC2 loadings for each ratio
loadings_ratios <- as.data.frame(pca_ratios$rotation[, 1:2]) %>%
  tibble::rownames_to_column("Ratio") %>%
  setNames(c("Ratio", "PC1", "PC2"))

# 2) plot them
quartz()
ggplot(loadings_ratios, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_text_repel(aes(label = Ratio),
                  size = 3,
                  box.padding   = 0.3,
                  point.padding = 0.2) +
  labs(
    title = "PCA of Lipid Ratios (loadings)",
    x     = paste0("PC1 (", round(summary(pca_ratios)$importance[2,1]*100,1), "%)"),
    y     = paste0("PC2 (", round(summary(pca_ratios)$importance[2,2]*100,1), "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )




###############################################################################
## INDIVIDUAL-SPECIES  PCA  (add this right after you create `wide_log`)
###############################################################################

## ---------------------------------------------------------------------------
## 1. build a “long” table that keeps every lipid species separately
##    – the only change is that we DON’T collapse to Class.
## ---------------------------------------------------------------------------
species_long <- combined %>%                         # <- same   raw table
  filter(str_detect(Compound_Name, class_pat)) %>%   # keep only valid classes
  pivot_longer(-c(Compound_Name, Condition),
               names_to  = "Sample",
               values_to = "Intensity") %>%
  group_by(Sample) %>%
  mutate(TIC      = sum(Intensity, na.rm = TRUE),
         rel      = Intensity / TIC,
         eps      = 0.5 * min(rel[rel > 0], na.rm = TRUE),
         log_rel  = log10(rel + eps)) %>%
  ungroup() %>%
  dplyr::select(Sample, Condition,
         Variable = Compound_Name,   # <- each lipid *species*
         Value    = log_rel)

## ---------------------------------------------------------------------------
## 2. pivot wide  →  run PCA  →  plot
##    (identical code structure you already used for classes & ratios)
## ---------------------------------------------------------------------------
wide_sp <- species_long %>%
  pivot_wider(names_from = Variable, values_from = Value)

meta_sp <- wide_sp %>% dplyr::select(Sample, Condition)
X_sp    <- wide_sp %>% dplyr::select(-Sample, -Condition)

# drop all-NA or zero-variance columns
ok_sp   <- sapply(X_sp, function(z) !all(is.na(z)) && var(z, na.rm = TRUE) > 0)
X_sp    <- X_sp[, ok_sp, drop = FALSE]

if(ncol(X_sp) >= 2){
  pca_sp <- prcomp(X_sp, center = TRUE, scale. = TRUE)
  
  scores_sp <- as.data.frame(pca_sp$x[,1:2]) %>%
    bind_cols(meta_sp)
  
  pct_sp <- round(summary(pca_sp)$importance[2,1:2] * 100, 1)
  
  library(ggplot2)
  library(cowplot)
  
  p_species <- ggplot(scores_sp,
                      aes(PC1, PC2, colour = Condition)) +
    geom_point(size = 2.5) +
    stat_ellipse(aes(fill = Condition),
                 geom = "polygon", alpha = 0.2, show.legend = FALSE) +
    labs(title = "PCA: individual lipid species",
         x = paste0("PC1 (", pct_sp[1], "%)"),
         y = paste0("PC2 (", pct_sp[2], "%)")) +
    theme_bw(base_size = 11)
  
  ## ---- Combine all three panels -------------------------------------------
  plot_grid(p_species,          # newly-made species plot
            last_plot(),        # your class-level PCA (already displayed)
            pca_ratios_plot,    # your ratio-level PCA (replace with object name)
            ncol = 3, align = "h")
} else {
  message("Species matrix collapsed to < 2 usable variables – no species PCA produced.")
}
