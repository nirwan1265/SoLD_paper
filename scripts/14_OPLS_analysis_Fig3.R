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
library(ropls)
library(ggbreak)
library(ggplot2)
library(dplyr)
library(grid)
library(knitr) 
library(kableExtra)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(patchwork)
library(ggplot2)
library(ggbreak)
library(broom); library(effsize)
library(lme4); library(lmerTest)


# ╔══════════════════════════════════════════════════════════════════╗
# ║ 1) READ RAW INTENSITY TABLES AND CLASS                           ║
# ╚══════════════════════════════════════════════════════════════════╝

# Read the raw files 
control <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  select(-c(2,3,4))
colnames(control)[1] <- "Compound_Name"  # Rename the first column to "Compound_Name"
dim(control)

lowinput <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  select(-c(2,3,4))
colnames(lowinput)[1] <- "Compound_Name"  # Rename the first column to "Compound_Name"

colnames(control)

# Valid classes
valid_classes <- c("DGDG","DGTS","MGDG",
                   "TG","DG","MG",
                   "PC","PE",
                   "SQDG",
                   "LPC","LPE","PG","PA","PS"
                   )

class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 2) TIC normalisation + log10 per lipid                           ║
# ╚══════════════════════════════════════════════════════════════════╝

# Long format
long_prep_global <- function(df) {
  df_long <- df %>%
    pivot_longer(-c(Compound_Name, Condition),
                 names_to = "Lipid", values_to = "Intensity") %>%
    dplyr::rename(Sample = Compound_Name)
  
  # 1) TIC per sample → relative abundance
  df_long <- df_long %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(TIC       = sum(Intensity, na.rm = TRUE),
           rel_abund = Intensity / TIC) %>%
    ungroup()
  
  # 2) small pseudo-count per sample → log10
  df_long <- df_long %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(minpos = ifelse(
      all(is.na(rel_abund) | rel_abund <= 0),
      NA_real_, min(rel_abund[rel_abund > 0], na.rm = TRUE)),
      eps    = ifelse(is.na(minpos), 0, minpos * 0.5)) %>%
    ungroup() %>%
    dplyr::mutate(log_rel = log10(rel_abund + eps))
  
  # 3) summarise mean log10(relative abundance) per lipid class
  df_long <- df_long %>%
    dplyr::mutate(Class = str_extract(Lipid, class_pat)) %>%
    filter(!is.na(Class)) %>%
    dplyr::group_by(Sample, Condition, Class) %>%
    dplyr::summarise(class_log = mean(log_rel, na.rm = TRUE), .groups = "drop")
  
  return(df_long)
}

# combine Control + LowInput, then run prep
combined  <- bind_rows(control  %>% dplyr::mutate(Condition = "Control"),
                       lowinput %>% dplyr::mutate(Condition = "LowInput"))
colnames(combined)
long_all  <- long_prep_global(combined)
str(control)

# Remove na:
long_all <- long_all %>%
  filter(!is.na(class_log))
table(is.na(long_all))

# ╔══════════════════════════════════════════════════════════════════╗
# ║ 3)  PIVOT WIDE                                                   ║
# ╚══════════════════════════════════════════════════════════════════╝

# Wide format
wide_log <- long_all %>%
  pivot_wider(names_from = Class, values_from = class_log,
              values_fill = NA_real_)



# 1) pivot into “long” so we have one row per Sample/Condition/Class
long_classes <- wide_log %>%
  pivot_longer(
    cols      = -c(Sample, Condition),
    names_to  = "Class",
    values_to = "Abundance"
  )

# 2) self-join on Sample+Condition, then keep only one direction of each pair
ratios_long <- long_classes %>%
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
unique(ratios_long$RatioName)

# Remove rows with string count less than 3 from Sample
ratios_long <- ratios_long %>% 
  filter(nchar(Sample) >= 5)

table(ratios_long$Sample)


# find the “full” count
full_n <- ratios_long %>% 
  group_by(Sample) %>% 
  tally() %>% 
  pull(n) %>% 
  max()

# keep only those with the full complement
ratios_clean <- ratios_long %>% 
  group_by(Sample) %>% 
  filter(n() == full_n) %>% 
  ungroup()


# 1) Pivot your ratios into a (n_samples × n_features) matrix
wide_ratios <- ratios_clean %>% 
  # make one row per Sample×Condition
  pivot_wider(
    id_cols   = c(Sample, Condition),
    names_from = RatioName,
    values_from = Ratio
  ) %>% 
  arrange(Sample, Condition)

# keep the metadata
meta <- wide_ratios %>% dplyr::select(Sample, Condition)

# pull out the numeric feature matrix
lipid_mat <- wide_ratios %>% 
  dplyr::select(-Sample, -Condition) %>% 
  as.matrix()
dim(lipid_mat)
# 2) build your class vector: a factor of your two Conditions
class_vec <- factor(meta$Condition)

# sanity check
stopifnot(nrow(lipid_mat) == length(class_vec))

# 3) fit your OPLS
set.seed(123)  # for reproducibility of any CV‐based choices
quartz()
opls_mod <- opls(
  lipid_mat, 
  class_vec,
  predI  = 1,       # one predictive component
  orthoI = NA,      # have opls choose orthogonal comps by CV
  scaleC = "standard"  # mean‐center & unit‐variance scale
)

# 4) take a look
str(opls_mod)

scores_df <- data.frame(
  Sample    = rownames(opls_mod@scoreMN),
  t1        = opls_mod@scoreMN[, 1],
  to1       = opls_mod@orthoScoreMN[, 1],
  Condition = class_vec,
  stringsAsFactors = FALSE
)


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

# 2) Make the OPLS-DA score plot
opls_plot <- ggplot(scores_df, aes(x = t1, y = to1, color = Condition)) +
  
  # dashed center lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  
  # samples + ellipses
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(fill = Condition),
               geom   = "polygon",
               type   = "norm",
               level  = 0.95,
               alpha  = 0.2,
               colour = NA) +
  
  # color/fill palette (same as PCA plot)
  scale_color_manual(
    values = c(Control = "#440154FF", LowInput = "#FDE725FF")
  ) +
  scale_fill_manual(
    values = c(Control = "#440154FF", LowInput = "#FDE725FF")
  ) +
  
  # fixed aspect and expand limits a bit
  coord_fixed(
    xlim = c(min(scores_df$t1) * 1.05, max(scores_df$t1) * 1.05),
    ylim = c(min(scores_df$to1) * 1.05, max(scores_df$to1) * 1.05)
  ) +
  
  # labels with % of X‐variance on t1
  labs(
    #title = "OPLS-DA Scores",
    x     = paste0(
      "t1 (", 
      round(opls_mod@summaryDF$`R2X(cum)` * 100, 1), 
      "% of X variance)"
    ),
    y     = "to1"
  ) +
  
  # styling
  theme_bw(base_size = 14) +
  theme_classic(base_size = 14) +
  
  # legend inside the panel
  theme(
    legend.position      = c(0.95, 0.95),        # 95% from left, 95% from bottom
    legend.justification = c("right", "top"),     # anchor legend at its top‐right
    legend.background    = element_rect(
      fill   = "white",
      color  = "grey70",
      size   = 0.4,
      linetype = "solid"
    ),
    legend.direction     = "vertical",
    legend.spacing.y     = unit(0.2, "cm")
  ) +
  plot_theme

quartz()
print(opls_plot)

# Save the plot
#ggsave("fig/main/Fig2c_OPLS_scores_plot.png", plot = opls_plot, width = 10, height = 6, dpi = 300, bg = "white")

# 2) Run the permutation test (200 label-swaps)
perm_res <- opls(
  lipid_mat, class_vec,
  predI = 1,
  permI = 500,
  scaleC = "standard"
)

# 3) Pull out the observed metrics
sumDF   <- opls_mod@summaryDF
obsR2X  <- sumDF$`R2X(cum)`  # e.g. 0.514
obsR2Y  <- sumDF$`R2Y(cum)`  # e.g. 0.903
obsQ2   <- sumDF$`Q2(cum)`   # e.g. 0.901

# 4) Extract and clean the permuted distributions
perm_mat   <- as.data.frame(perm_res@suppLs$permMN)
colnames(perm_mat) <- c("R2Xcum","R2Ycum","Q2cum","RMSEE","pre","ort","unused")
perm_only  <- perm_mat[-1, ]    # drop row 1 (the real model)

# 5) Compute “exact” one-sided p-values
nPerm     <- nrow(perm_only)
pR2Y_ex   <- (sum(perm_only$R2Ycum >= obsR2Y) + 1) / (nPerm + 1)
pQ2_ex    <- (sum(perm_only$Q2cum  >= obsQ2 ) + 1) / (nPerm + 1)

# 6) Tidy the permuted metrics for plotting
plot_df <- perm_only %>%
  select(R2Ycum, Q2cum) %>%
  pivot_longer(everything(),
               names_to  = "metric",
               values_to = "value") %>%
  mutate(metric = recode(metric, R2Ycum = "R²Y", Q2cum = "Q²"))

# 7) Build a small table of observed values + p-values
obs_df <- tibble(
  metric = c("R²Y","Q²"),
  obs    = c(obsR2Y, obsQ2),
  p_ex   = c(pR2Y_ex, pQ2_ex)
)

# 8) Plot!
# # 2) Build the plot with the break
# perm_plot_broken_nature <- ggplot(plot_df, aes(x = value)) + 
#   
#   # Histograms (fill mapped here!)
#   geom_histogram(aes(fill = metric),
#                  position = "identity",
#                  alpha    = 0.6,
#                  bins     = 30,
#                  color    = "white") + 
#   
#   # Dashed lines (manually colored)
#   geom_vline(data = obs_df,
#              aes(xintercept = obs),
#              linetype = "dashed",
#              size     = 0.8,
#              color    = c("#440154FF", "#FDE725FF")) + 
#   
#   # Labels (also manually colored)
#   geom_text(data = obs_df,
#             aes(x = obs, y = Inf,
#                 label = sprintf("%s = %.3f\np = %.3f", metric, obs, p_ex)),
#             hjust = -0.1,
#             vjust = c(1.2, 3.5),
#             size  = 3,
#             color = c("#440154FF", "#FDE725FF")) +
#   
#   # R²X arrow and label
#   geom_segment(aes(x = obsR2X, xend = obsR2X, y = 0, yend = -5),
#                arrow = arrow(length = unit(0.2, "cm")),
#                color = "darkgreen") +
#   annotate("text",
#            x     = obsR2X,
#            y     = -8,
#            label = sprintf("R²X = %.3f", obsR2X),
#            color = "darkgreen",
#            size  = 3,
#            hjust = 0.5) +
#   
#   # Manual fill legend
#   scale_fill_manual(
#     values = c("R²Y" = "#440154FF", "Q²" = "#FDE725FF"),
#     labels = c("Perm R²Y", "Perm Q²")
#   ) +
#   
#   labs(x = "Metric value", y = "Frequency") +
#   
#   scale_x_break(c(0.05, 0.90), scales = 0.5) +
#   
#   guides(
#     fill = guide_legend(
#       override.aes = list(alpha = 0.6),
#       title = NULL
#     )
#   ) +
#   theme(
#     legend.position      = "top",
#     legend.justification = "center",
#     legend.background    = element_rect(
#       fill     = "white",
#       color    = "grey70",
#       size     = 0.4,
#       linetype = "solid"
#     ),
#     legend.direction     = "vertical",
#     legend.spacing.y     = unit(0.2, "cm")
#   ) +
#   plot_theme
# 
# 
# 
# # 3) Draw it
# quartz()
# print(perm_plot_broken_nature)

# Summarize permuted distributions
perm_summary <- plot_df %>%
  group_by(metric) %>%
  summarise(
    mean_perm = mean(value, na.rm = TRUE),
    sd_perm   = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

# Merge with observed
bar_df <- perm_summary %>%
  left_join(obs_df, by = "metric")

pal <- c("R²Y"="#440154FF","Q²"="#FDE725FF")

make_hist_zoom <- function(metric_name, bw = 0.002, col = "#440154FF") {
  df  <- dplyr::filter(plot_df, metric == metric_name)
  obs <- dplyr::filter(obs_df,  metric == metric_name)$obs
  xl  <- min(df$value, na.rm = TRUE)
  xh  <- max(df$value, na.rm = TRUE)
  pad <- 0.01
  
  ggplot(df, aes(value)) +
    geom_histogram(binwidth = bw, fill = scales::alpha(col, 0.85), colour = "white") +
    coord_cartesian(xlim = c(xl - pad, xh + pad)) +
    # arrow + label pointing to the (off-panel) observed value near 1
    annotate("segment", x = xh + pad*0.9, xend = xh + pad*0.4,
             y = Inf, yend = Inf, colour = col,
             arrow = arrow(length = unit(0.18,"cm")), lineend = "round") +
    annotate("text", x = xh + pad*0.4, y = Inf,
             label = sprintf("obs = %.3f\np = %.3f", obs,
                             dplyr::filter(obs_df, metric == metric_name)$p_ex),
             colour = col, hjust = 1, vjust = 1.2, size = 4) +
    labs(x = paste0(metric_name, " (permutation)"), y = "Count") +
    theme_bw(24) +
    plot_theme
}

p_r2y_zoom <- make_hist_zoom("R²Y", bw = 0.002, col = "#440154FF")
p_q2_zoom  <- make_hist_zoom("Q²",  bw = 0.002, col = "#FDE725FF")

quartz()
perm_plot <- (p_r2y_zoom | p_q2_zoom)
perm_plot

# Save the plot
#ggsave("fig/supp/SuppFig6_OPLS_permutation_plot.png", plot = perm_plot, width = 12, height = 8, dpi = 300, bg = "white")



# Extract VIP scores
vip_vec <- slot(opls_mod, "vipVn")     # numeric vector with names = lipid IDs

vip_df <- tibble(
  Lipid = names(vip_vec),
  VIP   = as.numeric(vip_vec)
)

# 2.  Keep VIP > 1.3  (and sort descending)
vip_hits <- vip_df %>%
  filter(VIP > 1) %>%
  arrange(desc(VIP))

# How many?
cat("Number of discriminatory lipids (VIP > 1):", nrow(vip_hits), "\n")
print(vip_hits, n = Inf)




library(dplyr)
library(forcats)
library(ggplot2)

# --- your buckets ---
mem_sulfolipid   <- c("PC_SQDG","PE_SQDG","PG_SQDG","PS_SQDG","DGDG_SQDG","DG_SQDG","MG_SQDG")
mem_galactolipid <- c("DGDG_MG","DGDG_MGDG","DGDG_PG","DGDG_PE","DGDG_PS","MG_MGDG","MGDG_PC","MGDG_PE","MGDG_PG", "MGDG_PS")
mem_phospholipid <- c("PC_PE","PC_PS","PA_PG","PA_PS","PA_PE","PE_PS","PG_PS")

turn_diacylglycerol <- c("DG_MG","DG_PC","DG_PE","DG_PG","DG_PS","DG_TG","MG_PG","MG_PE","MG_PC","MG_PA","MG_PS")
turn_lysophospho    <- c("LPC_LPE","LPE_PA","LPC_PE","LPC_PS","LPC_MG","LPE_PS","LPC_PG","DGDG_LPE","LPE_MG","LPE_MGDG","LPE_SQDG","DG_LPE")
turn_triacylglycerol <- c("LPC_TG","MGDG_TG","PA_TG","PS_TG","SQDG_TG")




group_map <- c(
  setNames(rep("Sulfolipid adjustments",   length(mem_sulfolipid)),   mem_sulfolipid),
  setNames(rep("Galactolipid dynamics",    length(mem_galactolipid)), mem_galactolipid),
  setNames(rep("Phospholipid homeostasis", length(mem_phospholipid)), mem_phospholipid),
  setNames(rep("DG/MG turnover / signaling",  length(turn_diacylglycerol)), turn_diacylglycerol),
  setNames(rep("Lyso-phospho turnover",    length(turn_lysophospho)), turn_lysophospho),
  setNames(rep("TG turnover / signaling",  length(turn_triacylglycerol)), turn_triacylglycerol)
)

vip_df <- vip_hits %>%
  mutate(Key   = gsub("/", "_", Lipid),
         Group = unname(group_map[Key]),
         Group = ifelse(is.na(Group), "Other", Group)) %>%
  mutate(Group = factor(
    Group,
    levels = c("Sulfolipid adjustments","Galactolipid dynamics",
               "Phospholipid homeostasis","DG/MG turnover / signaling",
               "Lyso-phospho turnover","TG turnover / signaling","Other")
  ))

# reorder within each facet without extra deps
vip_df <- vip_df %>%
  group_by(Group) %>%
  mutate(Lipid_f = fct_reorder(Lipid, VIP, .desc = FALSE)) %>%
  ungroup()

p_vip_facets <- ggplot(vip_df, aes(x = Lipid_f, y = VIP, fill = Group)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2, show.legend = FALSE) +
  coord_flip() +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_fill_viridis_d(end = 0.92) +
  labs(x = NULL, y = "VIP score", title = "OPLS-DA VIP by pathway group") +
  theme_bw(base_size = 11) +
  theme(
    strip.placement    = "outside",
    strip.text.y.left  = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )
quartz()
p_vip_facets






# --- your buckets ---

group_map <- c(
  setNames(rep("Sulfolipid",   length(mem_sulfolipid)),   mem_sulfolipid),
  setNames(rep("Galactolipid",    length(mem_galactolipid)), mem_galactolipid),
  setNames(rep("Phospholipid", length(mem_phospholipid)), mem_phospholipid),
  setNames(rep("Mono/Diacylglycerol",  length(turn_diacylglycerol)), turn_diacylglycerol),
  setNames(rep("Lyso-phospho",    length(turn_lysophospho)), turn_lysophospho),
  setNames(rep("Triacylglycerol",  length(turn_triacylglycerol)), turn_triacylglycerol)
)

# ---------------------------
# 1) VIP table with groups + per-group order
# ---------------------------
vip_df <- vip_hits %>%
  mutate(Key   = gsub("/", "_", Lipid),
         Group = unname(group_map[Key]),
         Group = ifelse(is.na(Group), "Other", Group)) %>%
  mutate(Group = factor(
    Group,
    levels = c("Sulfolipid","Galactolipid",
               "Phospholipid","Mono/Diacylglycerol",
               "Lyso-phospho","Triacylglycerol","Other")
  )) %>%
  group_by(Group) %>%
  mutate(Lipid_f = fct_reorder(Lipid, VIP, .desc = FALSE)) %>%
  ungroup()

# LOOKUP (single source of truth for the box side)
lookup <- vip_df %>% select(Lipid, Group, Lipid_f)

# ---------------------------
# 2) Build ratio_tbl from wide_log ONLY for VIP ratios
#    (wide_log is log-space, so log-ratio = num - den)
# ---------------------------
build_ratio_tbl <- function(wide_log, vip_df) {
  ratios <- unique(vip_df$Lipid)
  parts  <- tibble(ratio = ratios) %>%
    separate(ratio, into = c("num","den"), sep = "/", remove = FALSE) %>%
    mutate(col = gsub("/", "_", ratio))
  
  # optional sanity check
  need <- unique(c(parts$num, parts$den))
  miss <- setdiff(need, names(wide_log))
  if (length(miss)) warning("Missing in wide_log: ", paste(miss, collapse = ", "))
  
  out <- wide_log %>% select(Sample, Condition)
  for (i in seq_len(nrow(parts))) {
    n  <- parts$num[i]; d <- parts$den[i]; cn <- parts$col[i]
    out[[cn]] <- wide_log[[n]] - wide_log[[d]]
  }
  out %>% select(Sample, Condition, all_of(parts$col))
}

ratio_tbl <- build_ratio_tbl(wide_log, vip_df)

# ---------------------------
# 3) Long data for distributions, stamped with VIP order/group
# ---------------------------
keep_keys  <- gsub("/", "_", vip_df$Lipid)

ratio_long <- ratio_tbl %>%
  pivot_longer(cols = all_of(keep_keys), names_to = "Key", values_to = "Value") %>%
  mutate(Lipid = gsub("_", "/", Key)) %>%
  left_join(lookup, by = "Lipid") %>%     # <- DO NOT recompute groups; trust VIP
  mutate(Condition = factor(Condition, levels = c("Control","LowInput"))) %>%
  filter(!is.na(Group))                   # should be none, but just in case



# ---------------------------
# 4) Plots
# ---------------------------

# Left: VIP bars (horizontal)
p_vip <- ggplot(vip_df, aes(x = VIP, y = Lipid_f, fill = Group)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2, show.legend = FALSE) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  scale_fill_viridis_d(end = 0.92) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "VIP score", y = NULL, title = "OPLS-DA VIP by pathway group") +
  theme_bw(base_size = 24) +
  theme(
    strip.placement    = "outside",
    strip.text.y.left  = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

# Right: LEFT-ALIGNED distributions (ratios on y; dodge by Condition)
p_boxes <- ggplot(ratio_long, aes(x = Lipid_f, y = Value, fill = Condition)) +
  geom_violin(position = position_dodge(width = 0.72),
              width = 0.9, alpha = 0.25, colour = NA, trim = TRUE) +
  geom_boxplot(position = position_dodge(width = 0.72),
               width = 0.55, outlier.shape = NA, colour = "black", alpha = 0.95) +
  coord_flip() +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c(Control = "#440154FF", LowInput = "#FDE725FF")) +
  labs(x = NULL, y = "Z-score") +
  theme_bw(base_size = 24) +
  theme(
    # keep horizontal separators between ratios (after flip these are major.y)
    #panel.grid.major.y = element_line(colour = "grey90"),
    #panel.grid.minor.y = element_blank(),
    # optional faint vertical grid
    #panel.grid.major.x = element_line(colour = "grey95"),
    
    # nuke y-axis text/ticks (the ratio names)
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    
    # remove row facet strips on the right panel
    strip.text.y.left  = element_blank(),
    strip.background.y = element_blank(),
    
    legend.position = "bottom",
    #strip.placement = "outside",
    # Remove all grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 


# Stitch side-by-side
quartz()
p_out <- (p_vip | p_boxes) + plot_layout(widths = c(1, 2), guides = "collect")

# theme the *combined* plot
p_out <- p_out & theme(
  legend.position       = c(0.5, 1),   # top-center
  legend.justification  = c(0.5, 1),
  legend.direction      = "horizontal",
  legend.box.background = element_rect(color = "grey70", fill = "white", size = 0.6),
  legend.background     = element_blank(),
  legend.key            = element_rect(fill = NA, colour = NA),
  legend.margin         = margin(4, 6, 4, 6)
)
p_out

# High-res PNG (journal friendly)
ggsave("VIP_by_group_boxes.png", plot = p_out,
       width = 12, height = 24, units = "in", dpi = 300, bg = "white")












# ---------------------------
### STATISTISTICS:
# ---------------------------


# ── 1) z-scores for plotting ONLY (within each Condition × Ratio) ─────────────
ratio_long2 <- ratio_long

ratio_long <- ratio_long %>%
  group_by(Lipid, Condition) %>%
  mutate(Z = as.numeric(scale(Value))) %>%
  ungroup()

# (use `Z` instead of `Value` in your violin/box plot if you want standardized viz)

# ── 2) robust stats on raw log-ratios (Value) ─────────────────────────────────

# helper that returns Wilcoxon, HL shift, AUC (prob. superiority), Cliff’s Δ,
# and jackknife sign stability of the median difference
one_ratio_tests <- function(df, ratio, alternative = "two.sided") {
  xi <- df %>% filter(Condition == "LowInput", Lipid == ratio) %>% pull(Value)
  x0 <- df %>% filter(Condition == "Control",  Lipid == ratio) %>% pull(Value)
  
  if (length(xi) < 3 || length(x0) < 3) return(
    tibble(Ratio = ratio, n_LI = length(xi), n_C = length(x0),
           median_LI = median(xi, na.rm=TRUE), median_C = median(x0, na.rm=TRUE),
           HL_shift = NA_real_, HL_low = NA_real_, HL_high = NA_real_,
           p_wilcox = NA_real_, AUC = NA_real_, cliffs_delta = NA_real_,
           jackknife_stability = NA_real_)
  )
  
  wt <- wilcox.test(xi, x0, alternative = alternative,
                    conf.int = TRUE, exact = FALSE)  # HL shift + CI
  
  # U and AUC (probability an LI value > a C value; for 'less' interpret accordingly)
  U   <- unname(wt$statistic)
  auc <- as.numeric(U) / (length(xi) * length(x0))
  
  cd  <- tryCatch(effsize::cliff.delta(xi, x0)$estimate, error = \(e) NA_real_)
  
  # jackknife sign stability on medians
  d_full <- sign(median(xi) - median(x0))
  N      <- length(xi) + length(x0)
  keep   <- logical(N)
  for (i in seq_len(N)) {
    xi2 <- xi; x0_2 <- x0
    if (i <= length(xi)) xi2 <- xi2[-i] else x0_2 <- x0_2[-(i - length(xi))]
    keep[i] <- sign(median(xi2) - median(x0_2)) == d_full
  }
  stab <- mean(keep)
  
  tibble(
    Ratio   = ratio,
    n_LI    = length(xi),
    n_C     = length(x0),
    median_LI = median(xi),
    median_C  = median(x0),
    HL_shift   = unname(wt$estimate),   # robust median (LI − C) on log10 scale
    HL_low     = wt$conf.int[1],
    HL_high    = wt$conf.int[2],
    p_wilcox   = wt$p.value,
    AUC        = auc,                   # P(LI > C)
    cliffs_delta = cd,                  # −1..1
    jackknife_stability = stab          # 0..1
  )
}

# which ratios to test?
#   - all:            unique(ratio_long$Lipid)
#   - only VIP ones:  vip_df$Lipid
ratios_to_test <- unique(ratio_long$Lipid)

# If you want ONE-SIDED direction for a few known cases:
greater_in_LI <- c("PS/SQDG","PG/SQDG","PE/SQDG","PC/SQDG","DGDG/SQDG")  # expect LI > C
less_in_LI    <- c("SQDG/TG")                                            # expect LI < C
dir_vec <- setNames(rep("two.sided", length(ratios_to_test)), ratios_to_test)
dir_vec[intersect(ratios_to_test, greater_in_LI)] <- "greater"
dir_vec[intersect(ratios_to_test, less_in_LI)]    <- "less"

# run the tests
ratio_stats <- map_dfr(
  ratios_to_test,
  \(r) one_ratio_tests(ratio_long, r, alternative = dir_vec[[r]])
) %>%
  mutate(p_adj_BH = p.adjust(p_wilcox, method = "BH")) %>%
  arrange(p_adj_BH)

# pretty columns
ratio_stats_out <- ratio_stats %>%
  mutate(
    effect_log10 = HL_shift,
    effect_fc    = 10^HL_shift,                    # fold-change of ratios LI/C
    AUC_pct      = scales::percent(AUC, accuracy = 0.1)
  ) %>%
  select(Ratio, n_C, n_LI, median_C, median_LI,
         effect_log10, effect_fc,
         HL_low, HL_high,
         p_wilcox, p_adj_BH, AUC_pct, cliffs_delta, jackknife_stability)

# save + view
write.csv(ratio_stats_out, "table/supp/SuppTable_ratio_stats_all.csv")
print(ratio_stats_out, n = 30)




colnames(ratio_stats_out)













# Convert to latex format:
# kable(vip_hits, format = "latex", booktabs = TRUE, linesep = "", escape = FALSE) %>%
#   kable_styling(latex_options = c("striped", "hold_position")) %>%
#   column_spec(1, bold = TRUE) %>%
#   column_spec(2, width = "3cm") %>%
#   row_spec(0, bold = TRUE, background = "#f7f7f7") %>%
#   row_spec(nrow(vip_hits), bold = TRUE, background = "#f7f7f7") %>%
#   row_spec(1:nrow(vip_hits)-1, background = "#f7f7f7") %>%
#   row_spec(1:nrow(vip_hits)-1, color = "black") %>%
#   row_spec(nrow(vip_hits), color = "black") %>%
#   row_spec(0, color = "black") %>%
#   row_spec(0, background = "#f7f7f7") %>%
#   row_spec(1:nrow(vip_hits)-1, background = "#f7f7f7") %>%
#   row_spec(1:nrow(vip_hits)-1, color = "black") %>%
#   row_spec(nrow(vip_hits), background = "#f7f7f7") %>%
#   row_spec(nrow(vip_hits), color = "black") %>%
#   row_spec(0, background = "#f7f7f7") %>%
#   row_spec(0, color = "black") %>%
#   row_spec(1:nrow(vip_hits)-1, background = "#f7f7f7") %>%
#   row_spec(1:nrow(vip_hits)-1, color = "black") %>%
#   row_spec(nrow(vip_hits), background = "#f7f7f7") %>%
#   row_spec(nrow(vip_hits), color = "black") %>%
#   row_spec(0, background = "#f7f7f7") %>%
#   row_spec(0, color = "black") %>%
#   row_spec(1:nrow(vip_hits)-1, background = "#f7f7f7")
# 

# pg_cv <- raw_peak_tbl %>%     # your raw species-level table
#   filter(Class == "PG") %>% 
#   group_by(SampleGroup) %>%   # QC or pooled injections
#   summarise(CV = sd(Intensity) / mean(Intensity)) 


# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ################################################################################
# #### SOrghum Lipidomics Database (SOLD)
# ################################################################################
# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 
# 
# ###############################################################################
# ## Supplementary Figures – Overlap of Control vs Low‑input lipid lists
# ## Three panels:
# ##   S1 – All detected lipids               (Venn)
# ##   S2 – Traditional classes only          (Venn)
# ##   S3 – Curated “Class” annotation set    (table + bar *optional*)
# ###############################################################################
# 
# # ------------------------------------------------------------------------------
# # 1. Load libraries
# # ------------------------------------------------------------------------------
# library(vroom)      # fast CSV import
# library(dplyr)      # data wrangling
# library(tidyr)      # pivoting
# library(ggplot2)    # plotting
# library(stringr)    # string manipulation
# library(forcats)    # factor reordering
# library(tidyverse)
# library(ggpubr)      # For stat_compare_means
# library(scales) 
# library(ggVennDiagram)   # nice, ggplot‑based 2‑set Venn
# library(readr) 
# library(sf)
# library(ggVennDiagram) # For Venn diagrams)
# library(viridis)
# library(ropls)    # OPLS‑DA
# library(limma)
# library(knitr)
# library(kableExtra)
# 
# # ------------------------------------------------------------------------------
# # 2. Read in the data
# # ------------------------------------------------------------------------------
# #control  <- vroom("data/SPATS_fitted/BLUP_GWAS_phenotype/control_all_lipids_BLUPs_rename.csv")
# control  <- vroom("data/SPATS_fitted/non_normalized_intensities/control_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
# colnames(control)[1] <- "Compound_Name"  
# 
# # Select only partial string match to columns with strings TG, DG, PC, PE, PA, PS, DGDG, SQDG, MGDG, MG
# control  <- control %>% dplyr::select(Compound_Name, contains("TG("), contains("DG("), contains("PC("),
#                                       contains("PE("), contains("PA("), contains("PS("),
#                                       contains("DGDG("), contains("SQDG("), contains("MGDG("), contains("MG("))
# 
# # Remove rows with strin gcount less than 5 from Compound_Name
# control <- control %>% 
#   filter(nchar(Compound_Name) >= 5)
# 
# # Lowinput
# #lowinput  <- vroom("data/SPATS_fitted/BLUP_GWAS_phenotype/lowinput_all_lipids_BLUPs_rename.csv")
# lowinput  <- vroom("data/SPATS_fitted/non_normalized_intensities/lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
# colnames(lowinput)[1] <- "Compound_Name"  
# 
# 
# # Select only partial string match to columns with strings TG, DG, PC, PE, PA, PS, DGDG, SQDG, MGDG, MG
# lowinput  <- lowinput %>% dplyr::select(Compound_Name, contains("TG("), contains("DG("), contains("PC("),
#                                           contains("PE("), contains("PA("), contains("PS("),
#                                           contains("DGDG("), contains("SQDG("), contains("MGDG("), contains("MG("))
# 
# 
# # Remove rows with strin gcount less than 5 from Compound_Name
# lowinput <- lowinput %>% 
#   filter(nchar(Compound_Name) >= 5)
# 
# 
# 
# #### CHANGE THE NAME TO COMMONNAME
# # Lipid class
# # lipid_class_info <- vroom::vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/data/lipid_class.csv", show_col_types = FALSE) %>%
# #   dplyr::filter(!is.na(Class))
# lipid_class_info <- vroom::vroom("data/lipid_class/lipid_classes.csv", show_col_types = FALSE) %>% dplyr::distinct(Lipids, .keep_all = TRUE)  # keep unique Lipids with their CommonName
# 
# # Create a named vector of replacements: names = original, values = new names
# name_map <- lipid_class_info %>%
#   dplyr::filter(!is.na(CommonName)) %>%
#   dplyr::select(Lipids, CommonName) %>%
#   deframe()  # turns it into a named vector: Lipids -> CommonName
# 
# # Function to rename lipid columns in any dataset
# rename_lipid_columns <- function(df, name_map) {
#   original_names <- colnames(df)
#   
#   renamed <- original_names %>%
#     # Replace only if a match exists
#     sapply(function(x) if (x %in% names(name_map)) name_map[[x]] else x)
#   
#   # Apply renamed vector back
#   colnames(df) <- renamed
#   return(df)
# }
# 
# # Apply to both control and lowinput
# control  <- rename_lipid_columns(control,  name_map)
# lowinput <- rename_lipid_columns(lowinput, name_map)
# 
# 
# 
# # 0.  Ensure identical lipid columns in both data frames
# shared_lipids <- intersect(names(control)[-1],   # exclude Compound_Name
#                            names(lowinput)[-1])
# 
# ctrl  <- control  %>% dplyr::select(Compound_Name, all_of(shared_lipids))
# lowin <- lowinput %>% dplyr::select(Compound_Name, all_of(shared_lipids))
# 
# 
# # 1.  Assemble one long table with a Condition flag
# full_tbl <- bind_rows(
#   ctrl  %>% mutate(Condition = "Control"),
#   lowin %>% mutate(Condition = "LowInput")
# )
# 
# 
# # 2.  Build lipid_mat  (samples × lipids)  and  class_vec
# 
# lipid_mat <- full_tbl %>%
#   dplyr::select(-Compound_Name, -Condition) %>%   # leave only numeric columns
#   as.matrix()
# 
# class_vec <- factor(full_tbl$Condition, levels = c("Control", "LowInput"))
# 
# 
# # 3.  (Optional) Pareto‑scale or autoscale before OPLS‑DA
# #  -- Pareto scaling: mean‑centre, divide by sqrt(sd)
# lipid_mat <- scale(lipid_mat, center = TRUE, scale = sqrt(apply(lipid_mat, 2, sd)))
# 
# # Check dimensions
# dim(lipid_mat)   #  ➜  samples × lipids
# table(class_vec) #  ➜  sample counts per condition
# lipid_mat[1:5,1:5]
# str(lipid_mat)
# 
# # 1. Pull out only the numeric matrix
# #    (drop any non-numeric columns if you have e.g. lipid names in the first column)
# lipid_numeric <- as.matrix(lipid_mat)  
# lipid_numeric <- t(lipid_numeric)
# dim(lipid_numeric)
# # if lipid_mat has a "Compound_Name" column, do:
# # lipid_numeric <- as.matrix(lipid_mat %>% select(-Compound_Name))
# 
# 
# # 2. Build a simple factor/vector of batch labels
# batch_vec <- c(
#   rep("Control", 394),
#   rep("LowInput", 363)
# )
# 
# # check it matches the number of columns
# stopifnot(length(batch_vec) == ncol(lipid_numeric))
# 
# # 2.  Run OPLS‑DA  -----------------------------------------------------
# set.seed(123)
# quartz()
# opls_mod <- opls(lipid_mat,
#                  class_vec,
#                  predI    = 1,       # 1 predictive component (default)
#                  orthoI   = NA,      # ropls chooses # orthogonal comps via CV
#                  scaleC   = "standard")  # autoscale (mean‑center, unit variance)
# 
# # 3.  Model diagnostics  ----------------------------------------------
# summary(opls_mod)      # R2X, R2Y, Q2, # components
# 
# quartz()
# plot(opls_mod, type = "overview")      # score + loading plots
# quartz()
# plot(opls_mod, type = "correlation")   # VIP vs correlation
# 
# # 4.  Permutation test (over‑fit check) -------------------------------
# perm_res <- opls(lipid_mat, class_vec,
#                  permI = 200,         # 200 permuted models
#                  predI = 1)
# 
# # look at p‑value in perm_res@summaryDF
# perm_res@summaryDF
# 
# 
# # 1.  Extract VIP scores
# 
# vip_vec <- slot(opls_mod, "vipVn")     # numeric vector with names = lipid IDs
# 
# vip_df <- tibble(
#   Lipid = names(vip_vec),
#   VIP   = as.numeric(vip_vec)
# )
# 
# # 2.  Keep VIP > 1.3  (and sort descending)
# 
# vip_hits <- vip_df %>%
#   filter(VIP > 1.2) %>%
#   arrange(desc(VIP))
# 
# # How many?
# cat("Number of discriminatory lipids (VIP > 1.3):", nrow(vip_hits), "\n")
# 
# 
# # 3.  Annotate with lipid class
# lookup_tbl <-
#   lipid_class_info %>%
#   dplyr::select(Lipid = Lipids, SubClass) %>%
#   mutate(Lipid = str_trim(Lipid)) %>%
#   distinct(Lipid, .keep_all = TRUE)
# 
# 
# # Compare with vip hits
# vip_hits <- vip_hits %>%
#   left_join(lookup_tbl, by = "Lipid") %>%          # add Class column
#   mutate(SubClass = ifelse(is.na(SubClass), "Unk", SubClass))
# 
# 
# # 4A.  Save as Supplementary Table (CSV)
# write.csv(vip_hits,
#           "VIP_lipids_OPLSDA_Control_vs_LowInput.csv",
#           row.names = FALSE)
# 
# # 4B.  Publication bar plot  (top 30 VIPs)
# vip_top30 <- vip_hits %>% slice_head(n = 30)
# 
# quartz()
# ggplot(vip_top30,
#        aes(x = reorder(Lipid, VIP),
#            y = VIP,
#            fill = SubClass)) +
#   geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
#   coord_flip() +
#   geom_hline(yintercept = 1.3, linetype = "dashed", colour = "red") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
#   scale_fill_viridis_d(option = "D", direction = 1) +
#   labs(
#     x     = NULL,
#     y     = "VIP score",
#     title = "Top 30 lipids driving Control vs Low-input separation\n(OPLS-DA, VIP > 1.3)"
#   ) +
#   theme_bw(base_size = 9) +
#   theme(
#     axis.line          = element_line(colour = "black"),
#     panel.grid.major.y = element_blank(),
#     panel.grid.major.x = element_line(colour = "grey90"),
#     panel.grid.minor   = element_blank(),
#     legend.position    = "right",
#     legend.title       = element_blank(),
#     plot.title         = element_text(hjust = 0.5, face = "plain", size = 10)
#   )













### VOLCANO
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# — 1) bring in your ratio matrix and VIP vector —  
# assume `lipid_mat` is an R matrix or data.frame with 710 rows (samples) × 78 columns (ratios),
# with rownames(lipid_mat) = sample IDs, and
# vip_vec is a named numeric vector length 78, names(vip_vec) = the same ratios.

# turn the matrix into a tibble, preserving rownames
ratios_df <- as_tibble(lipid_mat, rownames = "SampleID")

# add a Condition column (first 394 control, then 316 lowinput)
ratios_df <- ratios_df %>%
  mutate(Condition = rep(c("Control","LowInput"), times = c(394,316)))

# — 2) pivot longer so each row is one measurement —  
long_df <- ratios_df %>%
  pivot_longer(
    -c(SampleID, Condition),
    names_to  = "Ratio",
    values_to = "Value"
  )

# — 3) keep only those Ratios that have VIP scores —  
vip_df <- enframe(vip_vec, name="Ratio", value="VIP")
long_df <- long_df %>%
  inner_join(vip_df, by="Ratio")

# — 4) compute statistics for each Ratio —  
stats_df <- long_df %>%
  group_by(Ratio, VIP) %>%
  summarise(
    # mean per condition
    mean_ctrl = mean(Value[Condition=="Control"], na.rm=TRUE),
    mean_low  = mean(Value[Condition=="LowInput"], na.rm=TRUE),
    # log2 fold‐change
    log2FC    = log2(mean_low / mean_ctrl),
    # p‐value from two‐sample t‐test
    pval      = t.test(
      Value[Condition=="LowInput"],
      Value[Condition=="Control"]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    negLog10P = -log10(pval),
    highVIP   = VIP > 1
  )

# — 5) volcano plot —  
quartz()
ggplot(stats_df, aes(x = log2FC, y = negLog10P)) +
  geom_point(aes(color = highVIP), alpha = 0.7, size=2) +
  scale_color_manual(
    values = c("FALSE"="gray70", "TRUE"="firebrick"),
    labels = c("VIP ≤ 1","VIP > 1"),
    name   = "High VIP"
  ) +
  geom_vline(xintercept = c(-0.05,0.05), linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  labs(
    x = "log₂ fold‐change (LowInput / Control)",
    y = "-log₁₀(p-value)",
    title = "Volcano of lipid‐ratio changes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill="white", color="black")
  )



library(ggrepel)

# define your thresholds
fc_cut   <- 0.001
p_cut    <- 0.001
vip_cut  <- 1.0

stats_df <- stats_df %>%
  mutate(
    sig = (abs(log2FC) >= fc_cut & pval < p_cut),
    vip_hi = VIP >= vip_cut
  )
quartz()
ggplot(stats_df, aes(x=log2FC, y=negLog10P)) +
  geom_point(aes(color = vip_hi), alpha=0.7, size=2) +
  scale_color_manual(
    values = c("FALSE"="gray70","TRUE"="firebrick"),
    labels = c(paste0("VIP < ",vip_cut), paste0("VIP ≥ ",vip_cut)),
    name   = "High VIP"
  ) +
  geom_vline(xintercept = c(-fc_cut, fc_cut), linetype="dashed", color="black") +
  geom_hline(yintercept = -log10(p_cut),    linetype="dashed", color="black") +
  geom_text_repel(
    data = filter(stats_df, sig & vip_hi),
    aes(label = Ratio),
    size = 3,
    max.overlaps = 20
  ) +
  labs(
    x = "log₂ fold-change (LowInput / Control)",
    y = "-log₁₀(p-value)",
    title = "Volcano of lipid-ratio changes"
  ) +
  theme_minimal(base_size=14) +
  theme(
    legend.position = c(0.8,0.8),
    legend.background = element_rect(fill="white", color="black")
  )
