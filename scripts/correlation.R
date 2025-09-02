
###############################################################################
# 1.  LIBRARIES
###############################################################################

suppressPackageStartupMessages({
  library(vroom);   library(dplyr);   library(tidyr);  library(stringr)
  library(ggplot2); library(viridis); library(gridExtra); library(grid)
})


################################################################################
################################################################################
#### FIGURE 1A
################################################################################
################################################################################


################################################################################
# 1.  LOAD DATA  (Control & LowInput)
################################################################################

control <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  dplyr::select(-c(2,3,4)) %>%
  dplyr::rename(Compound_Name = LineRaw)

lowinput <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  dplyr::select(-c(2,3,4)) %>%
  dplyr::rename(Compound_Name = LineRaw) 


################################################################################
# 2.  SELECT THE CLASSES 
################################################################################

valid_classes <- c(
  "TG",
  
  "DGDG","MGDG","SQDG",
  
  "PC","PE","PI","PG","PA","PS","LPC","LPE",
  
  "DG","MG"
)
#class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")
class_pat <- paste0("^(", paste(valid_classes, collapse = "|"), ")")

reshape_plate <- function(df, label) {
  df %>%
    pivot_longer(-Compound_Name, names_to = "Lipid", values_to = "Intensity") %>%
    # now using our reordered, “longest‐first” regex
    mutate(Class = str_extract(Lipid, class_pat)) %>%
    filter(!is.na(Class)) %>%
    group_by(Sample = Compound_Name, Class) %>%
    summarise(sum_int = sum(Intensity, na.rm = TRUE), .groups="drop") %>%
    mutate(Condition = label)
}

ctrl_long <- reshape_plate(control,  "Control")
low_long  <- reshape_plate(lowinput, "LowInput")

unique(ctrl_long$Class)
unique(low_long$Class)


# assumes you already have `ctrl_long` and `low_long` as in your script
# columns: Sample, Class, sum_int, Condition

library(dplyr)
library(tidyr)
library(ggplot2)

## 1) keep only samples present in BOTH conditions
common_ids <- intersect(unique(ctrl_long$Sample), unique(low_long$Sample))

ctrl_c <- ctrl_long %>%
  filter(Sample %in% common_ids) %>%
  select(Sample, Class, sum_int) %>%
  rename(sum_ctrl = sum_int)

low_c <- low_long %>%
  filter(Sample %in% common_ids) %>%
  select(Sample, Class, sum_int) %>%
  rename(sum_low = sum_int)

## 2) pivot to wide (samples x classes), fill missing with 0
ctrl_wide <- ctrl_c %>%
  group_by(Sample, Class) %>% summarise(sum_ctrl = sum(sum_ctrl), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Class, values_from = sum_ctrl, values_fill = 0) %>%
  arrange(Sample)

low_wide <- low_c %>%
  group_by(Sample, Class) %>% summarise(sum_low = sum(sum_low), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Class, values_from = sum_low, values_fill = 0) %>%
  arrange(Sample)

## 3) make sure rows (samples) line up
stopifnot(identical(ctrl_wide$Sample, low_wide$Sample))

## 4) log-transform (raw-sum on log scale)
ctrl_mat <- as.matrix(ctrl_wide[, -1])
low_mat  <- as.matrix(low_wide [, -1])

#ctrl_log <- log10(ctrl_mat + 1)
#low_log  <- log10(low_mat  + 1)
ctrl_log <- ctrl_mat
low_log  <- low_mat
## 5) compute class-by-class correlations (control vs low-input)
# cor(X, Y) => matrix of corr between columns of X and columns of Y
Cmat <- cor(ctrl_log, low_log, use = "pairwise.complete.obs", method = "pearson")

## 6) long format for ggplot
cor_long <- as.data.frame(as.table(Cmat))
colnames(cor_long) <- c("Class_ctrl", "Class_low", "r")

## 7) optional: order classes in a nice sequence
class_order <- c("TG","DG","MG","PC","PE","PG","PI","PA","PS","LPC","LPE","MGDG","DGDG","SQDG")
class_order <- intersect(class_order, colnames(ctrl_log))  # keep only those that exist
cor_long$Class_ctrl <- factor(cor_long$Class_ctrl, levels = class_order)
cor_long$Class_low  <- factor(cor_long$Class_low,  levels = class_order)

## 8) heatmap
p_heat <- ggplot(cor_long, aes(x = Class_low, y = Class_ctrl, fill = r)) +
  geom_tile(color = "grey90", linewidth = 0.2, na.rm = FALSE) +
  # emphasize the diagonal (same class vs same class)
  geom_tile(
    data = subset(cor_long, as.character(Class_ctrl) == as.character(Class_low)),
    fill = NA, color = "black", linewidth = 0.5
  ) +
  # optional: print numbers
  geom_text(aes(label = ifelse(is.na(r), "", sprintf("%.2f", r))),
            size = 3, color = "black") +
  scale_fill_gradient2(
    name = "Pearson r",
    limits = c(-1, 1), midpoint = 0,
    low = "#2c7bb6", mid = "white", high = "#d7191c",
    na.value = "grey85"
  ) +
  coord_equal() +
  labs(
    title = "Control vs Low-input: class-by-class correlation (log10(sum intensity + 1))",
    x = "Lowinput",
    y = "Control"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    plot.title   = element_text(face = "bold")
  )

quartz()
print(p_heat)

## 9) if you also want the *same-class* r's as a quick table:
same_class_r <- cor_long %>%
  filter(as.character(Class_ctrl) == as.character(Class_low)) %>%
  arrange(desc(r))
print(same_class_r)
