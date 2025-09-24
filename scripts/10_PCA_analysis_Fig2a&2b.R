
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
library(ggtext)
library(ggpp)
library(textshape)
library(viridis)
library(RColorBrewer)
library(tibble)

###############################################################################
## 2) READ RAW INTENSITY TABLES AND CLASS             
###############################################################################

# Read the raw files 
control <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  select(-c(2,3,4)) %>%
  rename(Compound_Name = 1) %>%
  mutate(Condition = "Control")

lowinput <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  select(-c(2,3,4)) %>%
  rename(Compound_Name = 1) %>%
  mutate(Condition = "LowInput")

colnames(control)

###############################################################################
## 3)  SUBSET TRADITIONAL LIPID CLASSES
###############################################################################

# Valid classes
valid_classes <- c("TG","DG","MG",
                   "PC","PE",
                   "DGDG","MGDG",
                   "SQDG",
                   "LPC","LPE","PG","PA","PS")



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
#shared_lipids <- intersect(lipid_cols_ctrl, lipid_cols_low)
shared_lipids <- intersect(colnames(control), colnames(lowinput))
shared_lipids <- shared_lipids[!shared_lipids %in% c("Compound_Name", "Condition")]

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


# 3) scale
# build matrix of shared lipids (your code up to pca_input is fine)
mat <- lipid_mat %>% dplyr::select(all_of(shared_lipids)) %>% as.matrix()

# 1) per-sample TIC normalization
tic <- rowSums(mat, na.rm = TRUE)
rel <- sweep(mat, 1, tic, "/")

# 2) per-sample pseudocount (half the smallest positive per sample)
minpos <- apply(rel, 1, function(x) if (any(x > 0, na.rm=TRUE)) min(x[x > 0], na.rm=TRUE) else NA_real_)
eps    <- 0.5 * minpos
rel_eps <- rel
for (i in seq_len(nrow(rel_eps))) {
  if (!is.na(eps[i]) && eps[i] > 0) rel_eps[i, ] <- rel_eps[i, ] + eps[i]
}

# 3) log10 transform
log_rel <- log10(rel_eps)


# 4) Run PCA
pca_res <- prcomp(log_rel, center = TRUE, scale. = TRUE)

# 3) Extract scores and merge in metadata
scores_df <- data.frame(pca_res$x[,1:2], Sample_ID = rownames(pca_res$x))
scores_df <- merge(
  scores_df,
  pca_input %>% dplyr::select(Sample_ID, Compound_Name, Condition),
  by = "Sample_ID",
  all.x = TRUE
)
scores_df_individual <- scores_df

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
  mutate(
    LipidClass = str_extract(Class, class_pattern)   # pull off the prefix
  )

# 3) build a named colour vector
# build a Paired palette of exactly as many colors as you have classes
n_classes <- length(valid_classes)
qual_cols <- brewer.pal(n = max(3, min(12, n_classes)), name = "Paired")
# Change TG to black
qual_cols[which(valid_classes == "TG")] <- "#000000"  # black for TG

# if you have more than 12 classes, you can interpolate:
if(n_classes > length(qual_cols)){
  qual_cols <- colorRampPalette(qual_cols)(n_classes)
}
names(qual_cols) <- valid_classes



# Theme
plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(
      size   = 14,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x   = element_text(
      size = 16,      # X‐axis title size
      face = "bold"
    ),
    axis.title.y   = element_text(
      size = 16,      # Y‐axis title size
      face = "bold"
    ),
    axis.text.x    = element_text(
      size = 16,      # X‐axis tick label size
      color = "black"
    ),
    axis.text.y    = element_text(
      size = 16,      # Y‐axis tick label size
      color = "black"
    ),
    axis.line      = element_line(color = "black"),
    panel.grid     = element_blank(),
    
    legend.position      = c(0.95, 0.95),
    legend.justification = c("right","top"),
    legend.background    = element_rect(fill="white", color="grey70", size=0.4),
    legend.direction     = "vertical",
    legend.spacing.y     = unit(0.2,"cm"),
    legend.title         = element_blank(),
    legend.text          = element_text(size=16),
    
    plot.margin    = margin(15, 15, 15, 15)
  )

# 5) Plot biplot

# With loadings
individual_PCA <- ggplot() +
  # samples
  geom_point(data = scores_df,
             aes(PC1, PC2, color = Condition),
             size = 3, alpha = 0.8) +
  stat_ellipse(data = scores_df,
               aes(PC1, PC2, fill = Condition),
               geom = "polygon", alpha = 0.2, colour = NA) +
  # loading arrows, coloured by lipid class
  geom_segment(data = load_df,
               aes(x = 0, y = 0, xend = a1, yend = a2, color = LipidClass),
               arrow = arrow(length = unit(0.2, "cm")),
               size = 0.8) +
  geom_text(data = load_df,
            aes(x = a1 * 1.05, y = a2 * 1.05, label = Class, color = LipidClass),
            size = 2.5) +


  scale_fill_manual(name   = "Condition",
                    values = c(Control = "#440154FF", LowInput = "#FDE725FF")
                    #guide  = guide_legend(override.aes = list(shape = 22, size = 4, alpha = 0.))
  ) +
  # guides & palettes
  scale_color_manual(
    name   = "Lipid Class",
    values = qual_cols
  ) +
  # keep your sample‐point colours as before

  # axes, theme
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    title = "Biplot: Individual Lipid Species",
    x = paste0("PC1 (", round(100 * pca_res$sdev[1]^2 / sum(pca_res$sdev^2),1), "%)"),
    y = paste0("PC2 (", round(100 * pca_res$sdev[2]^2 / sum(pca_res$sdev^2),1), "%)")
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title     = element_text(face = "bold")
  )

quartz()
individual_PCA

# Without loadings. 
individual_PCA <- ggplot() +
  # Samples + ellipses
  geom_point(data = scores_df,
             aes(PC1, PC2, color = Condition),
             size = 2.5, alpha = 0.8) +
  stat_ellipse(data = scores_df,
               aes(PC1, PC2, fill = Condition),
               geom = "polygon", alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    x = paste0("PC1 (", round(100 * pca_res$sdev[1]^2 / sum(pca_res$sdev^2),1), "%)"),
    y = paste0("PC2 (", round(100 * pca_res$sdev[2]^2 / sum(pca_res$sdev^2),1), "%)")
  ) +
  scale_color_manual(values = c(Control = "#440154FF", LowInput = "#FDE725FF")) +
  scale_fill_manual(values  = c(Control = "#440154FF", LowInput = "#FDE725FF")) +
  
  coord_fixed(xlim = c(min(scores_df$PC1)*1.05, max(scores_df$PC1)*1.05),
              ylim = c(min(scores_df$PC2)*1.05, max(scores_df$PC2)*1.05)) +
  
  theme_bw(base_size = 16) +  # base theme with box and grid
  
  theme(
    panel.grid.major = element_line(color = "grey85", size = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
    
    plot.title = element_text(face = "bold"),
    
    legend.position      = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(
      fill     = "white",
      color    = "grey70",
      size     = 0.4,
      linetype = "solid"
    ),
    legend.direction     = "vertical",
    legend.spacing.y     = unit(0.2, "cm")
  )  +

  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
   
  ) +
  plot_theme 

quartz()
individual_PCA

# Save the plot
#ggsave("fig/main/Fig2a_individual_indv_lipid_PCA_plot.png", individual_PCA, width = 6, height = 6, dpi = 300, bg = "white")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###############################################################################
## SUMMED LIPIDS
###############################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# ╔══════════════════════════════════════════════════════════════════╗
# ║ 1) READ RAW INTENSITY TABLES AND CLASS                           ║
# ╚══════════════════════════════════════════════════════════════════╝

class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

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
# Remove rows with string count less than 5 from Sample
long_all <- long_all %>% 
  filter(nchar(Sample) >= 5)
# ╔══════════════════════════════════════════════════════════════════╗
# ║ 3)  PIVOT WIDE                                                   ║
# ╚══════════════════════════════════════════════════════════════════╝

# Wide format
wide_log <- long_all %>%
  pivot_wider(names_from = Class, values_from = class_log,
              values_fill = NA_real_)

# Removing  AEG and PG 
#wide_log <- wide_log[,-c(3,13)]
#wide_log <- wide_log[,-c(3,8)]


long_classes <- wide_log %>%
  pivot_longer(
    cols      = -c(Sample, Condition),
    names_to  = "Class",
    values_to = "Abundance"
  )

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


# 3) Scale
# 1) per-sample TIC normalization
tic <- rowSums(mat, na.rm = TRUE)
rel <- sweep(mat, 1, tic, "/")

# 2) per-sample pseudocount (half the smallest positive per sample)
minpos <- apply(rel, 1, function(x) if (any(x > 0, na.rm=TRUE)) min(x[x > 0], na.rm=TRUE) else NA_real_)
eps    <- 0.5 * minpos
rel_eps <- rel
for (i in seq_len(nrow(rel_eps))) {
  if (!is.na(eps[i]) && eps[i] > 0) rel_eps[i, ] <- rel_eps[i, ] + eps[i]
}

# 3) log10 transform
log_rel <- log10(rel_eps)





# 3) Run PCA (on scaled data)
p <- prcomp(log_rel, center = T, scale. = T)

# 4) PCA scores (Sample positions)
scores_df <- as.data.frame(p$x[, 1:2]) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(sample_info, by = "Sample_ID")  # attach Condition
scores_df_summed <- scores_df

# 5) Loadings (Lipid-class arrows)
load_df <- as.data.frame(p$rotation[, 1:2]) %>%
  rownames_to_column("Class") #%>%
  # mutate(a1 = PC1 * 5,  # scale for visibility
  #        a2 = PC2 * 5)

scale_factor <- 5
# compute a data‐driven scale factor
sf <- max(abs(scores_df$PC1), abs(scores_df$PC2)) /
  max(sqrt(load_df$PC1^2 + load_df$PC2^2))
# now multiply loadings by, say, 0.8*sf
load_df <- load_df %>% mutate(
  a1 = PC1 * (0.8 * sf),
  a2 = PC2 * (0.8 * sf)
)


# 6) Plot
# left side legends
plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(
      size   = 14,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x   = element_text(
      size = 16,      # X‐axis title size
      face = "bold"
    ),
    axis.title.y   = element_text(
      size = 16,      # Y‐axis title size
      face = "bold"
    ),
    axis.text.x    = element_text(
      size = 16,      # X‐axis tick label size
      color = "black"
    ),
    axis.text.y    = element_text(
      size = 16,      # Y‐axis tick label size
      color = "black"
    ),
    axis.line      = element_line(color = "black"),
    panel.grid     = element_blank(),
    
    legend.position      = c(0.3, 0.95),
    legend.justification = c("right","top"),
    legend.background    = element_rect(fill="white", color="grey70", size=0.4),
    legend.direction     = "vertical",
    legend.spacing.y     = unit(0.2,"cm"),
    legend.title         = element_blank(),
    legend.text          = element_text(size=16),
    
    plot.margin    = margin(15, 15, 15, 15)
  )
quartz()

summed_pca <- ggplot() +
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
               color = "black", size = 0.8) +
  
  # loading labels
  geom_text(data = load_df,
            aes(x = a1 * 1.05, y = a2 * 1.05, label = Class),
            color = "black", size = 3) +
  
  # axis lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  
  # axis labels
  labs(
    x = paste0("PC1 (", round(100 * p$sdev[1]^2 / sum(p$sdev^2), 1), "%)"),
    y = paste0("PC2 (", round(100 * p$sdev[2]^2 / sum(p$sdev^2), 1), "%)")
  ) +
  
  # condition colors
  scale_color_manual(values = c(Control = "#440154FF", LowInput = "#FDE725FF")) +
  scale_fill_manual(values  = c(Control = "#440154FF", LowInput = "#FDE725FF")) +
  
  coord_fixed() +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85", size = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
    
    plot.title = element_text(face = "bold"),
    
    #legend.position      = c(0.95, 0.95),
    #legend.justification = c("right", "top"),
    legend.background    = element_rect(
      fill     = "white",
      color    = "grey70",
      size     = 0.4,
      linetype = "solid"
    ),
    legend.direction     = "vertical",
    legend.spacing.y     = unit(0.2, "cm")
  ) +
  plot_theme


summed_pca

# Save
#ggsave("fig/main/Fig2b_summed_lipid_PCA_biplot.png", summed_pca,width = 8, height = 8, dpi = 300, bg = "white")










### 3D PCA
# ---- 3D PCA with loadings (plotly) ----
library(dplyr)
library(tibble)
library(plotly)
library(htmlwidgets)

# variance explained for axis titles
var_exp <- (p$sdev^2) / sum(p$sdev^2)

# scores (PC1:PC3)
scores3d <- as.data.frame(p$x[, 1:3]) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(sample_info, by = "Sample_ID")

# loadings (PC1:PC3)
load3d <- as.data.frame(p$rotation[, 1:3]) %>%
  rownames_to_column("Class")

# scale loadings to sit nicely inside the score cloud
score_rad <- max(sqrt(rowSums(scores3d[,c("PC1","PC2","PC3")]^2)))
load_rad  <- max(sqrt(rowSums(load3d[,c("PC1","PC2","PC3")]^2)))
sf <- 0.8 * (score_rad / load_rad)  # 0.8 keeps a little margin

load3d <- load3d %>%
  mutate(a1 = PC1 * sf, a2 = PC2 * sf, a3 = PC3 * sf)

pal <- c(Control = "#440154FF", LowInput = "#FDE725FF")

fig <- plot_ly()

# sample points
fig <- fig %>%
  add_markers(
    data = scores3d,
    x = ~PC1, y = ~PC2, z = ~PC3,
    color = ~Condition, colors = pal,
    marker = list(size = 4, opacity = 0.75),
    hoverinfo = "text",
    text = ~paste0(
      "Sample: ", Sample_ID,
      "<br>PC1: ", round(PC1, 2),
      "<br>PC2: ", round(PC2, 2),
      "<br>PC3: ", round(PC3, 2),
      "<br>Condition: ", Condition
    ),
    showlegend = TRUE,
    name = "Samples"
  )

# loading arrows (one trace per arrow)
for (i in seq_len(nrow(load3d))) {
  fig <- fig %>%
    add_trace(
      type = "scatter3d", mode = "lines",
      x = c(0, load3d$a1[i]),
      y = c(0, load3d$a2[i]),
      z = c(0, load3d$a3[i]),
      line = list(color = "grey20", width = 4),
      hoverinfo = "text",
      text = paste0(load3d$Class[i]),
      showlegend = FALSE
    ) %>%
    add_trace(
      type = "scatter3d", mode = "text",
      x = load3d$a1[i] * 1.05,
      y = load3d$a2[i] * 1.05,
      z = load3d$a3[i] * 1.05,
      text = load3d$Class[i],
      textposition = "middle right",
      textfont = list(color = "black", size = 10),
      showlegend = FALSE
    )
}

fig <- fig %>%
  layout(
    title = list(text = "PCA (PC1–PC3) with loading vectors"),
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(100*var_exp[1],1), "%)"),
                   zeroline = TRUE, showspikes = FALSE),
      yaxis = list(title = paste0("PC2 (", round(100*var_exp[2],1), "%)"),
                   zeroline = TRUE, showspikes = FALSE),
      zaxis = list(title = paste0("PC3 (", round(100*var_exp[3],1), "%)"),
                   zeroline = TRUE, showspikes = FALSE)),
    legend = list(orientation = "h", x = 0.5, y = 1.05, xanchor = "center")
  )

# fig  # interactive 3D
#
# # Save an interactive HTML (great for supplement)
# saveWidget(as_widget(fig), "PCA_3D_loadings.html", selfcontained = TRUE)





# # palette (same as before)
# pal <- c(Control = "#440154FF", LowInput = "#FDE725FF")
#
# # general biplot for any PC pair
# plot_pca_pair <- function(p, sample_info, pcs = c(1, 2), arrow_scale = 0.8) {
#   # scores on selected PCs
#   sc <- as.data.frame(p$x[, pcs, drop = FALSE]) |>
#     tibble::rownames_to_column("Sample_ID") |>
#     dplyr::left_join(sample_info, by = "Sample_ID")
#   names(sc)[2:3] <- c("PCx", "PCy")
#
#   # loadings on selected PCs
#   ld <- as.data.frame(p$rotation[, pcs, drop = FALSE]) |>
#     tibble::rownames_to_column("Class")
#   names(ld)[2:3] <- c("PCx", "PCy")
#
#   # scale loadings to fit score cloud
#   score_rad <- max(sqrt(sc$PCx^2 + sc$PCy^2), na.rm = TRUE)
#   load_rad  <- max(sqrt(ld$PCx^2 + ld$PCy^2), na.rm = TRUE)
#   sf <- arrow_scale * score_rad / load_rad
#   ld <- mutate(ld, a1 = PCx * sf, a2 = PCy * sf)
#
#   # axis labels with variance explained
#   ve <- (p$sdev^2) / sum(p$sdev^2)
#   xlab <- sprintf("PC%d (%.1f%%)", pcs[1], 100 * ve[pcs[1]])
#   ylab <- sprintf("PC%d (%.1f%%)", pcs[2], 100 * ve[pcs[2]])
#
#   ggplot() +
#     geom_point(data = sc, aes(PCx, PCy, color = Condition), size = 2.5, alpha = 0.75) +
#     stat_ellipse(data = sc, aes(PCx, PCy, fill = Condition),
#                  geom = "polygon", alpha = 0.18, color = NA) +
#     geom_segment(data = ld, aes(x = 0, y = 0, xend = a1, yend = a2),
#                  arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.6) +
#     geom_text(data = ld, aes(x = a1 * 1.05, y = a2 * 1.05, label = Class),
#               color = "black", size = 3) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
#     labs(x = xlab, y = ylab) +
#     scale_color_manual(values = pal) +
#     scale_fill_manual(values  = pal) +
#     coord_fixed() +
#     theme_bw(base_size = 14) +
#     theme(
#       panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
#       panel.grid.minor = element_blank(),
#       legend.position  = "top",
#       legend.title     = element_blank()
#     )
# }
#
# # Make the two plots you asked for
# p13 <- plot_pca_pair(p, sample_info, pcs = c(1, 3))  # PC1–PC3
# p23 <- plot_pca_pair(p, sample_info, pcs = c(2, 3))  # PC2–PC3
#
# # (optional) show side-by-side if you use patchwork
# # library(patchwork)
# quartz()
# p13 | p23



# combined <- plot_grid(
#   individual_PCA +
#     theme(plot.margin = margin(5, 5, 5, 5)),   # tiny inner margins
#   summed_pca +
#     theme(plot.margin = margin(5, 5, 5, 5)),
#   ncol        = 2,
#   labels      = c("(c)", "(d)"),                # optional panel tags
#   align       = "hv",                           # align both horizontally & vertically
#   axis        = "tblr",                         # align top/bottom/left/right
#   rel_widths  = c(1, 1)                         # both panels equal width
# )
#
# quartz()
# combined






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
unique(ratios_long$RatioName)

# Remove rows with string count less than 3 from Sample
ratios_long <- ratios_long %>% 
  filter(nchar(Sample) >= 5)

unique(ratios_long$RatioName)

# Subset the ratios that you got from OPLS-DA from 12_OPLS_analysis.R
only_vips <- c("MG/SQDG",   "PS/SQDG"   ,"MG/MGDG"  ,"PG/SQDG"   ,"DG/MG"     ,"LPC/MG"   , "DGDG/MG"  , "MGDG/PS" , 
"PE/SQDG",   "LPC/PS"    ,"DG/PS"   ,  "LPE/SQDG" , "DGDG/PS"   ,"PC/PS"  ,   "DGDG/SQDG" ,"LPC/LPE" , 
"PE/PS" ,    "PG/PS"     ,"SQDG/TG"  , "MGDG/PG"   ,"MG/PG"     ,"PA/PS"   ,  "MG/PC"     ,"MGDG/PE"  ,
"PC/SQDG",   "LPE/MGDG"  ,"DGDG/MGDG" ,"MG/PA"     ,"PS/TG"     ,"LPE/PS"   , "MG/PE"     ,"MGDG/TG"  ,
"DG/LPE"  ,  "LPC/PE"    ,"LPE/PA"    ,"DG/TG"     ,"LPC/PG"    ,"PA/TG"    , "DG/SQDG"   ,"DG/PG"    ,
"MGDG/PC"  , "LPC/TG"    ,"DG/PE"     ,"DGDG/PE"   ,"MG/PS"     ,"DGDG/LPE" , "PA/PG"     ,"PA/PE"    ,
"LPE/MG"    ,"DGDG/PG"   ,"PC/PE" )

ratios_long <- ratios_long %>%
  filter(RatioName %in% only_vips)

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
scale_factor <- 5
# compute a data‐driven scale factor
sf <- max(abs(scores_df$PC1), abs(scores_df$PC2)) /
  max(sqrt(load_df$PC1^2 + load_df$PC2^2))
# now multiply loadings by, say, 0.8*sf
load_df <- load_df %>% mutate(
  a1 = PC1 * (0.8 * sf),
  a2 = PC2 * (0.8 * sf)
)

groups_lipids <- data.frame()

sulfolipid_adjustments   <- c("PC/SQDG","PE/SQDG","PG/SQDG","PS/SQDG","DGDG/SQDG","DG/SQDG","MG/SQDG")
galactolipid_dynamics <- c("DGDG/MG","DGDG/MGDG","DGDG/PG","DGDG/PE","DGDG/PS","MG/MGDG","MGDG/PC","MGDG/PE","MGDG/PG", "MGDG/PS")
phospholipid_homeostasis <- c("PC/PE","PC/PS","PA/PG","PA/PS","PA/PE","PE/PS","PG/PS")

turn_diacylglycerol <- c("DG/MG","DG/PC","DG/PE","DG/PG","DG/PS","DG/TG","MG/PG","MG/PE","MG/PC","MG/PA","MG/PS")
turnover_lysophospho    <- c("LPC/LPE","LPE/PA","LPC/PE","LPC/PS","LPC/MG","LPE/PS","LPC/PG","DGDG/LPE","LPE/MG","LPE/MGDG","LPE/SQDG","DG/LPE")
turn_triacylglycerol <- c("LPC/TG","MGDG/TG","PA/TG","PS/TG","SQDG/TG")


all_combine <- c(
  sulfolipid_adjustments,
  galactolipid_dynamics,
  phospholipid_homeostasis,
  turnover_lysophospho,
  turn_diacylglycerol,
  turn_triacylglycerol)


load_df <- load_df %>%
  mutate(
    Category = case_when(
      RatioName %in% sulfolipid_adjustments   ~ "Sulfolipid",
      RatioName %in% galactolipid_dynamics    ~ "Galactolipid",
      RatioName %in% phospholipid_homeostasis ~ "Phospholipid",
      RatioName %in% turnover_lysophospho      ~ "Lysophospholipid",
      RatioName %in% turn_triacylglycerol           ~ "Triacylglycerol",
      RatioName %in% turn_diacylglycerol ~ "Mono/Diacylglycerols",
      TRUE                                     ~ "Other"
    )
  )

# pick a qualititative palette (one colour per category)
cat_cols <- c(
  "Sulfolipid"   = "#1b9e77",
  "Galactolipid"    = "#d95f02",
  "Phospholipid" = "#000000",
  "Lysophospholipid" = "#e7298a",
  "Triacylglycerol" = "#7570b3",
  "Mono/Diacylglycerols"           = "#66a61e",
  "Other"                    = "grey70"
)


quartz()
ratio_pca <- ggplot() +
  # sample‐scores + ellipses
  geom_point(data = scores_df,
             aes(PC1, PC2, color = Condition),
             size = 2.5, alpha = 0.7, show.legend = FALSE) +    # hide from colour legend
  stat_ellipse(data = scores_df,
               aes(PC1, PC2, fill = Condition),
               geom = "polygon", alpha = 0.2, colour = NA) +
  
  # arrows coloured by your ratio‐Category
  geom_segment(data = load_df,
               aes(x = 0, y = 0,
                   xend = a1, yend = a2,
                   color = Category),
               arrow = arrow(length = unit(0.2, "cm")),
               size = 0.8) +
  geom_text(data = load_df,
            aes(x = a1 * 1.05, y = a2 * 1.05,
                label = RatioName,
                color = Category),
            size = 3) +
  
  # axes
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  coord_fixed() +
  
  # axes titles with % variance
  labs(title = "Biplot: Samples & Lipid Ratios Loadings",
       x = paste0("PC1 (", round(100 * pca_ratios$sdev[1]^2 /
                                 sum(pca_ratios$sdev^2), 1), "%)"),
       y = paste0("PC2 (", round(100 * pca_ratios$sdev[2]^2 /
                                 sum(pca_ratios$sdev^2), 1), "%)")) +
  
  # 1) Category colours only in the colour legend
  scale_color_manual(
    name   = "Category",
    values = cat_cols,
    breaks = names(cat_cols)      # only these show up in the legend
  ) +
  # 2) keep your sample fill legend for ellipses
  scale_fill_manual(
    name   = "Condition",
    values = c(Control  = "#440154FF",
               LowInput = "#FDE725FF")
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold")
  )

ratio_pca





# Save
# ggsave("figures/ratios_lipid_PCA_biplot.png", ratio_pca, width = 8, height = 6, dpi = 300, bg = "white")


# function to plot PCA biplot for a given category
plot_pca_loadings <- function(cat_name, scores_df, load_df, pca_ratios, cat_cols) {
  load_sub <- load_df %>% filter(Category == cat_name)
  
  ggplot() +
    # sample scores + ellipses
    geom_point(data = scores_df,
               aes(PC1, PC2, color = Condition),
               size = 2.5, alpha = 0.7, show.legend = FALSE) +
    stat_ellipse(data = scores_df,
                 aes(PC1, PC2, fill = Condition),
                 geom = "polygon", alpha = 0.2, colour = NA) +
    
    # loadings (arrows + labels)
    geom_segment(data = load_sub,
                 aes(x = 0, y = 0,
                     xend = a1, yend = a2),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = cat_cols[cat_name], size = 0.8) +
    geom_text(data = load_sub,
              aes(x = a1 * 1.05, y = a2 * 1.05, label = RatioName),
              color = cat_cols[cat_name], size = 3) +
    
    # axes
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    coord_fixed() +
    
    # variance explained
    labs(
      title = paste0("Biplot: ", cat_name, " Ratios"),
      x = paste0("PC1 (", round(100 * pca_ratios$sdev[1]^2 / sum(pca_ratios$sdev^2), 1), "%)"),
      y = paste0("PC2 (", round(100 * pca_ratios$sdev[2]^2 / sum(pca_ratios$sdev^2), 1), "%)")
    ) +
    scale_fill_manual(
      name   = "Condition",
      values = c(Control  = "#440154FF", LowInput = "#FDE725FF")
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      plot.title      = element_text(face = "bold")
    ) +
    plot_theme
}

quartz()
p_sulfo  <- plot_pca_loadings("Sulfolipid", scores_df, load_df, pca_ratios, cat_cols)
p_galac  <- plot_pca_loadings("Galactolipid", scores_df, load_df, pca_ratios, cat_cols)
p_phosph <- plot_pca_loadings("Phospholipid", scores_df, load_df, pca_ratios, cat_cols)
p_lyso   <- plot_pca_loadings("Lysophospholipid", scores_df, load_df, pca_ratios, cat_cols)
p_mdg    <- plot_pca_loadings("Mono/Diacylglycerols", scores_df, load_df, pca_ratios, cat_cols)
p_tg     <- plot_pca_loadings("Triacylglycerol", scores_df, load_df, pca_ratios, cat_cols)

library(patchwork)

quartz()
(p_sulfo | p_galac) / (p_phosph | p_lyso) / (p_mdg | p_tg)
