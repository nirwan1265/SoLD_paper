#################################################################################
###### LOAD THE LIBRARIES
#################################################################################

### Libraries
library(vroom)
library(dplyr)
library(caret)
library(mlr)
library(tuneRanger)
library(ggplot2)
library(treeshap)
library(shapviz)
library(patchwork)
library(stringr)
library(viridis) 
library(tidyr)
# Clear the list
rm(list=ls())

################################################################################
###### LOAD THE LIPIDS, PCA, AND PHENOTYPES
################################################################################

## Get the lipids, Log10 transformed
lipids <- vroom("data/SPATS_fitted/log10_median_centered/control_all_lipids_fitted_phenotype_log_centered.csv")
dim(lipids)
colnames(lipids)[1] <- "Line"

# 2) rename that long column to “GABA”
colnames(lipids)[
  colnames(lipids) == ".beta.-Phenyl-.gamma.-aminobutyric acid"
] <- "GABA"

# lipids <- lipids %>%
#   dplyr::select(Line,
#                 starts_with("TG("),
#                 starts_with("PC("),
#                 starts_with("LPC"),
#                 starts_with("SQDG"),
#                 starts_with("MGDG"),
#                 starts_with("DGDG"),
#                 # “MG(” will catch monogalactosyl diacylglycerol if it’s written as MG(...)
#                 starts_with("MG("),
#                 # “DG(” will catch diacylglycerols
#                 starts_with("DG(")
#   )

### Get the phenotypes in the field
# Our data: Plant height and Flowering Time
pheno <- vroom("data/phenotypes/control_field_phenotypes.csv") %>% dplyr::select(c(1,3))

# Plant Height and Stem Diameter
pheno <- vroom("data/phenotypes/plantheight_diameter_SAP.csv") %>% dplyr::select(c(1,2))

# Grain carotenoids
pheno <- vroom("data/phenotypes/grain_carotenoid_Clara_Cruet_Burgos.csv") %>% dplyr::select(c(1,7))

# Yield Traits
pheno <- vroom("data/phenotypes/yield_traits_Richard_E_Boyles.csv") %>% dplyr::select(c(1,4))

colnames(pheno)[1] <- "Line"
colnames(pheno)[2] <- "FlowerTime"

### Drop any rows with NA
keep_rows <- !is.na(pheno$FlowerTime)
pheno    <- pheno[keep_rows, ]
#pheno$FlowerTime <- log10(pheno$FlowerTime + 1)  # Log10 transform 
pheno$FlowerTime <- pheno$FlowerTime  # Log10 transform 

### Combine the phenotypes and lipids
dat <- lipids %>% inner_join(pheno, by = "Line")
dim(dat)


### Get the PCA - first 3 6
PC <- vroom("table/PCA_SAP.csv") %>% dplyr::select(c(1:4))
colnames(PC)[1] <- "Line"
dim(PC)
str(PC)

### Align and parse the phenotype as numeric  
common <- intersect(dat$Line, PC$Line)

dat2   <- dat %>% 
  filter(Line %in% common) %>% 
  arrange(Line) %>% 
  mutate(FlowerTime = as.numeric(FlowerTime))

PC2    <- PC  %>% 
  filter(Line %in% common) %>% 
  arrange(Line) %>% 
  select(-Line)           


################################################################################
###### POPULATION STRUCTURE ADJUSTMENT
################################################################################

### Build a pure‑numeric data.frame for residualization
df_resid <- data.frame(
  FlowerTime = dat2$FlowerTime,
  PC2   
)

### Residualise both the phenotype and the lipids on PCs so RF doesn’t learn population structure instead of biology:
adj_formula <- as.formula(paste("~", paste(colnames(PC2), collapse = " + ")))

### Residualize the phenotype on PCs  
rFT <- resid( lm(FlowerTime ~ ., data = df_resid) )

# Sanity check  
length(rFT)       # should equal nrow(df_resid)
summary(rFT)      # non‑zero residuals with mean ≈0

### Extract the lipid matrix (numeric only: drop Line & FlowerTime)
lipid_mat <- dat2 %>%
  select(-Line, -FlowerTime) %>%
  as.matrix()

### Build the design matrix X from your PCs only (no Line column)
#adj_formula still = "~ EV1 + EV2 + EV3 + EV4 + EV5"
X <- model.matrix(adj_formula, data = PC2)

### Compute the projection matrix P
P <- X %*% solve(crossprod(X)) %*% t(X)

### Residualize each lipid:  
lipids_adj_mat <- lipid_mat - (P %*% lipid_mat)

### Convert back to a tibble and restore column names
lipids_adj <- as_tibble(lipids_adj_mat)
colnames(lipids_adj) <- colnames(lipid_mat)

### Re‑join the Line IDs so you can keep track
lipids_adj <- bind_cols(Line = dat2$Line, lipids_adj)



################################################################################
###### STRATIFIED TRAIN/TEST SPLIT ON PCS, TUNING, AND MODEL FITTING
################################################################################

### Split the training and testing 80-20
set.seed(22)
lipids_num <- lipids_adj %>% select(-Line)
clusters  <- kmeans(PC2, centers = 6)$cluster
train_idx <- createDataPartition(clusters, p = 0.8, list = FALSE)

train_X <- lipids_num[train_idx, ]
train_y <- rFT      [train_idx]

test_X  <- lipids_num[-train_idx, ]
test_y  <- rFT      [-train_idx]


### Create an mlr regression task for tuning
train_df   <- data.frame(train_X, FlowerTime = train_y)
regr_task  <- makeRegrTask(data = train_df, target = "FlowerTime")


### Hyperparameter tuning with tuneRanger
tune_res <- tuneRanger(
  task            = regr_task,
  measure         = list(rmse),
  num.trees       = 1000,
  tune.parameters = c("mtry", "min.node.size", "sample.fraction"),
  num.threads     = 10
)

# SANITY CHECK
print(tune_res$recommended.pars)


################################################################################
###### 5‑FOLD CV ON THE TRAINING SET USING MLR
################################################################################

###Build the mlr learner with the tuned parameters
rf_learner <- makeLearner(
  "regr.ranger",
  par.vals = list(
    num.trees       = 1000,
    mtry            = tune_res$recommended.pars$mtry,
    min.node.size   = tune_res$recommended.pars$min.node.size,
    sample.fraction = tune_res$recommended.pars$sample.fraction
  ),
  predict.type = "response"
)

### Define the CV
rdesc <- makeResampleDesc("CV", iters = 5L)

### Run the resampling
cv_res <- resample(
  learner    = rf_learner,
  task       = regr_task,
  resampling = rdesc,
  measures   = list(rmse, mae, rsq),
  show.info  = TRUE
)

### Get the raw per‑fold metrics
df_meas <- as.data.frame(cv_res$measures.test)

### Add a simple iteration index
df_meas$Iteration <- seq_len(nrow(df_meas))

### Pivot to long form and clean up names
df_long <- pivot_longer(
  df_meas,
  cols = c("rmse", "mae", "rsq"),
  names_to = "Metric",
  values_to = "Value"
)

### Change to metrics labels
df_long$Metric <- recode(df_long$Metric,
                         "rmse.test.rmse" = "RMSE",
                         "mae.test.mean"  = "MAE",
                         "rsq.test.mean"  = "R²"
)

### Compute the mean for each metric
mean_df <- df_long %>%
  group_by(Metric) %>%
  summarize(Mean = mean(Value))

### Plot

# Pre-set theme
nature_theme <- theme_minimal(base_size = 16) +
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
    
    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(
      size = 16        # legend label size
    ),
    
    plot.margin    = margin(15, 15, 15, 15)
  )

mid_iter <- mean(df_long$Iteration)

cv_plot <- ggplot(df_long, aes(x = Iteration, y = Value, color = Metric)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  # dashed horizontal line at each mean
  geom_hline(data = mean_df,
             aes(yintercept = Mean, color = Metric),
             linetype = "dashed", size = 0.8) +
  # centered text at the middle iteration
  geom_text(data = mean_df,
            aes(x = mid_iter, y = Mean,
                label = sprintf("%s mean = %.3f", Metric, Mean),
                color = Metric),
            angle = 30,      # tilt a bit
            hjust = 0.5,     # center horizontally
            vjust = -0.5,    # nudge above the line
            size = 4) +
  scale_x_continuous(breaks = df_long$Iteration,
                     limits = c(min(df_long$Iteration), max(df_long$Iteration) + 1)) +
  scale_color_viridis_d(option = "D", end = 0.8) + 
  labs(
    title = "5 Fold CV Metrics by Iteration (RF)",
    x     = "CV iteration",
    y     = "Metric value"
  ) +
  nature_theme


quartz()
print(cv_plot)

# Save the plot
# ggsave("Fig5a_CV_RF_metrics.png",cv_plot, width = 8, height = 6, dpi = 300, 
#        units = "in", bg = "white")

################################################################################
###### FINAL MODEL FITTING AND EVALUATION
################################################################################

### Fit the final RF model with tuned parameters
rf_model <- ranger(
  x               = train_X,
  y               = train_y,
  num.trees       = 1000,
  mtry            = tune_res$recommended.pars$mtry,
  min.node.size   = tune_res$recommended.pars$min.node.size,
  sample.fraction = tune_res$recommended.pars$sample.fraction,
  importance      = "none",
  num.threads     = 10,
  seed            = 42
)

### Evaluate on the test set
preds <- predict(rf_model, data = test_X)$predictions
rmse  <- sqrt(mean((preds - test_y)^2))
print(paste0("Test RMSE: ", round(rmse, 4)))


### Make sure you’re using exactly the same test set for both:
length(preds)    # how many predictions did you get?
length(test_y)   # how many observations are in your test set?

# 1) RMSE (you already have)
rmse <- sqrt(mean((preds - test_y)^2))

# 2) MAE (mean absolute error)
mae  <- mean(abs(preds - test_y))

# 3) Pearson correlation & R²
pearson_r <- cor(preds, test_y)
r2        <- pearson_r^2

# 4) Bias (mean error)
bias <- mean(preds - test_y)

# 5) normalized RMSE (NRMSE) as % of observed SD
nrmse <- rmse / sd(test_y) * 100

# 6) summary table
metrics <- tibble(
  RMSE       = rmse,
  MAE        = mae,
  PearsonR   = pearson_r,
  R2         = r2,
  Bias       = bias,
  NRMSE_pct  = nrmse
)
print(metrics)


################################################################################
###### TreeSHAP ANALYSIS
################################################################################

### Extract feature matrix
X <- lipids_adj %>% select(-Line) %>% as.data.frame()

### Convert your ranger forest into a "unified" model
um <- unify(rf_model,X)

### treeshap wants an “unwrapped” ranger forest + the same data.frame
ts <- treeshap(um, X)
# ts$shaps is an (n_samples × n_features) matrix of exact SHAP values

### Mean absolute SHAP per lipid
mean_abs <- colMeans(abs(ts$shaps))
shap_rank <- tibble(
  Lipid       = names(mean_abs),
  MeanAbsSHAP = mean_abs
) %>%
  arrange(desc(MeanAbsSHAP))

### View top 20
top20 <- shap_rank %>% slice_head(n = 20)
print(top20)

### Convert to shapviz and plot
sv <- shapviz(ts, X = X)   # TreeSHAP under the hood

### ts is your treeshap result, ts$shaps is an (n_samples × n_features) matrix
shap_mat <- ts$shaps

### Build a global‐importance table by mean(|SHAP|)
global_rank <- tibble(
  Feature = colnames(shap_mat),
  MeanAbs = colMeans(abs(shap_mat))
) %>% 
  arrange(desc(MeanAbs))

### Grab the top 20
top50 <- global_rank %>% slice_head(n = 24)
print(top50, n = 24)

write.csv(global_rank, "table/SuppTable_lipid_SHAP_Flowering_time_Ruben_Rellan_Alvarez.csv", row.names = FALSE)


# 3) Bar chart of global importance
quartz()
ggplot(top50, aes(x = reorder(Feature, MeanAbs), y = MeanAbs)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x     = NULL,
    y     = "Mean |SHAP| (days)",
    title = "Top 50 Lipids by TreeSHAP Importance"
  ) +
  theme_minimal()

# 4) Prepare data for a beeswarm (only those top 20 features)
shap_long <- as_tibble(shap_mat) %>%
  mutate(Sample = row_number()) %>%
  pivot_longer(-Sample, names_to = "Feature", values_to = "SHAP") %>%
  filter(Feature %in% top50$Feature)

# 5) Beeswarm‐style scatter
quartz()
ggplot(shap_long, aes(x = SHAP, y = Feature)) +
  geom_jitter(height = 0.2, size = 1, alpha = 0.4) +
  labs(
    x     = "SHAP value (days)",
    y     = NULL,
    title = "Beeswarm of Top 50 Lipid SHAP Values"
  ) +
  theme_minimal()



### Residual distributions (raw vs PC‑residualized phenotype and an example lipid)
df_ft <- tibble(
  raw   = dat2$FlowerTime,
  resid = rFT
)

df_ft_long <- tibble(
  raw   = dat2$FlowerTime,
  resid = rFT
) %>%
  pivot_longer(
    cols      = c(raw, resid),
    names_to  = "Type",
    values_to = "Value"
  ) %>%
  # give each level a short code that we can re‐label
  mutate(Type = recode(Type,
                       raw   = "Raw FT",
                       resid = "Residual FT"))


residuals <- ggplot(df_ft_long, aes(x = Value, fill = Type)) +
  geom_density(alpha = 0.6, color = NA) +
  facet_wrap(
    ~Type,
    scales   = "free_x",
    labeller = labeller(Type = c(
      "Raw FT"      = "Flowering time (days)",
      "Residual FT" = "PC‑residualized FT (days)"
    ))
  ) +
  scale_fill_manual(
    values = c(
      "Raw FT"      = "#440154FF",
      "Residual FT" = "#FDE725FF"
    )
  ) +
  # label the axis so it's not blank
  labs(
    x = "Phenotypic value",
    y = "Density"
  ) +
  nature_theme +
  theme(legend.position = "none")

quartz()
print(residuals)

# Save the plot 
# ggsave("Fig5b_Residuals_phenotype.png", residuals, width = 8, height = 6, dpi = 300, 
#        units = "in", bg = "white")


# Some lipid (first column of lipids_adj after Line)
lipid_name <- colnames(lipids_adj)[251]

df_lip <- tibble(
  raw   = lipid_mat[,1],
  resid = lipids_adj[[lipid_name]]
)
# assume lipid_name is already set, and df_lip has columns `raw` and `resid`
df_lip_long <- df_lip %>%
  pivot_longer(
    cols      = c(raw, resid),
    names_to  = "Type",
    values_to = "Value"
  ) %>%
  mutate(
    Type = recode(
      Type,
      raw   = paste0("Raw ", lipid_name),
      resid = paste0("PC residualized ", lipid_name)
    )
  )


lipid_residuals <- ggplot(df_lip_long, aes(x = Value, fill = Type)) +
  geom_density(alpha = 0.6, color = NA) +
  facet_wrap(~ Type, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c(
    "Raw Zeaxanthin"              = "#440154FF",
    "PC residualized Zeaxanthin"  = "#FDE725FF"
  )) +
  labs(
    x = "Abundance",
    y = "Density"
  ) +
  nature_theme +
  theme(
    strip.text      = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

quartz()
print(lipid_residuals)

# Save the lipid residuals plot
# ggsave("Fig5c_Residuals_lipid.png", lipid_residuals, width = 8, height = 6, dpi = 300, 
#        units = "in", bg = "white")


### Observed vs Predicted FT
df_pred <- tibble(obs = test_y, pred = preds)
p_obs_pred <- ggplot(df_pred, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  labs(title = "Observed vs Predicted Flowering Time",
       x = "Observed residual FT",
       y = "RF predicted residual FT") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = sprintf("RMSE = %.3f\nR² = %.3f", metrics$RMSE, metrics$R2)) +
  nature_theme
quartz()
p_obs_pred

# Save the plot
# ggsave("Fig3d_obs_vs_pred.png", p_obs_pred, width = 6, height = 6, dpi = 300, units = "in", bg = "white")

### OOB error vs number of trees
ntree_seq <- c(100, 250, 500, 750, 1000, 1500)
oob_rmse   <- sapply(ntree_seq, function(nt) {
  m <- ranger(x = train_X, y = train_y,
              num.trees = nt,
              mtry = tune_res$recommended.pars$mtry,
              min.node.size = tune_res$recommended.pars$min.node.size,
              sample.fraction = tune_res$recommended.pars$sample.fraction,
              importance = "none",
              num.threads = 10,
              seed = 42)
  sqrt(m$prediction.error)
})
df_oob <- tibble(ntree = ntree_seq, OOB_RMSE = oob_rmse)
p_oob <- ggplot(df_oob, aes(x = ntree, y = OOB_RMSE)) +
  geom_line() + geom_point() +
  labs(title = "OOB RMSE vs Number of Trees",
       x = "Number of Trees", y = "OOB RMSE") +
  nature_theme

quartz()
p_oob

# save the plot
# ggsave("SuppFig_Hyperparameter_oob_vs_trees.png", width = 8, height = 6, dpi = 300)

### Hyperparameter tuning surface (mtry × min.node.size)
tuning_df <- tune_res$results

# Uncomment this
# p_tune <- ggplot(tuning_df, aes(x = factor(mtry), y = factor(min.node.size), fill = rmse)) +
#   geom_tile() +
#   scale_fill_viridis_c(name = "OOB RMSE") +
#   labs(title = "Hyperparameter Tuning Landscape",
#        x = "mtry", y = "min.node.size") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme_minimal(base_size = 16) +
#   nature_theme
# quartz()
# p_tune

# Save the plot
# ggsave("SuppFig_Hyperparameter_Tuning_surface.png", width = 16, height = 6, dpi = 300)



# ###Global SHAP importance & Beeswarm for top 50
top50 <- shap_rank %>% slice_head(n = 24)
# p_shap_bar <- ggplot(top50, aes(x = reorder(Lipid, MeanAbsSHAP), y = MeanAbsSHAP)) +
#   geom_col(fill = "steelblue") +
#   coord_flip() +
#   labs(title = "Top 20 Lipids by Mean |SHAP|",
#        x = NULL, y = "Mean |SHAP| (days)") +
#   nature_theme
# quartz()
# p_shap_bar
# ggsave("fig6_shap_bar.png", width = 6, height = 5, dpi = 300)
# 
# # Beeswarm
# shap_long <- as_tibble(ts$shaps) %>%
#   mutate(Sample = row_number()) %>%
#   pivot_longer(-Sample, names_to = "Lipid", values_to = "SHAP") %>%
#   filter(Lipid %in% top20$Lipid)
# 
# p_shap_bee <- ggplot(shap_long, aes(x = SHAP, y = reorder(Lipid, -SHAP))) +
#   geom_jitter(height = 0.2, alpha = 0.4, size = 1) +
#   labs(title = "Beeswarm of Top 50 Lipid SHAP Values",
#        x = "SHAP value (days)", y = NULL) +
#   theme_minimal()
# quartz()
# p_shap_bee


# prepare the bar chart (right panel)
p_shap_bar <- ggplot(top50, aes(
  x = reorder(Lipid, MeanAbsSHAP),
  y = MeanAbsSHAP
)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 50 Lipids by Mean |SHAP|",
    x     = NULL,
    y     = "Mean |SHAP| (days)"
  ) +
  nature_theme


# 1) compute cumulative‐importance and select features covering 80%
shap_cum <- shap_rank %>%
  arrange(desc(MeanAbsSHAP)) %>%
  mutate(cum_pct = cumsum(MeanAbsSHAP) / sum(MeanAbsSHAP))

# 1) pick the top‐80% lipids as a character vector
selected_lipids <- shap_cum %>% 
  filter(cum_pct <= 0.7183362) %>% 
  pull(Lipid)

# 2) reverse it once, to make the top feature end up at the top of a flipped axis
levels_order <- rev(selected_lipids)

# 2) bar chart of those lipids
p_shap_bar <- shap_rank %>% 
  filter(Lipid %in% levels_order) %>% 
  mutate(Lipid = factor(Lipid, levels = levels_order)) %>% 
  ggplot(aes(x = Lipid, y = MeanAbsSHAP)) +
  geom_col(fill = "#440154FF") +
  coord_flip() +
  labs(
    #title = "Top Lipids Covering 80% of SHAP Importance",
    x     = NULL,
    y     = "Mean |SHAP|"
  ) + 
  nature_theme +
  theme(
    axis.text.y  = element_blank(),  # remove the lipid names
    axis.ticks.y = element_blank()   # remove the tick marks
  )


# 3) beeswarm, using the exact same levels (so the top‐item is at the top)
shap_long <- as_tibble(ts$shaps) %>%
  mutate(Sample = row_number()) %>%
  pivot_longer(-Sample, names_to = "Lipid", values_to = "SHAP") %>%
  left_join(
    train_X %>% 
      as_tibble() %>%
      mutate(Sample = row_number()) %>%
      pivot_longer(-Sample, names_to = "Lipid", values_to = "FeatureValue"),
    by = c("Sample","Lipid")
  ) %>%
  filter(Lipid %in% selected_lipids) %>%
  mutate(Lipid = factor(Lipid, levels = rev(selected_lipids)))  # reverse so top is at top

# Beeswarm
shap_long2 <- as_tibble(ts$shaps) %>%
  mutate(Sample = row_number()) %>%
  pivot_longer(-Sample, names_to = "Lipid", values_to = "SHAP") %>%
  left_join(
    train_X %>% 
      as_tibble() %>% 
      mutate(Sample = row_number()) %>% 
      pivot_longer(-Sample, names_to = "Lipid", values_to = "FeatureValue"),
    by = c("Sample","Lipid")
  ) %>%
  filter(Lipid %in% levels_order) %>%
  mutate(Lipid = factor(Lipid, levels = levels_order))

p_shap_bee <- ggplot(shap_long2, aes(
  x     = SHAP, 
  y     = Lipid, 
  color = FeatureValue
)) +
  geom_jitter(height = 0.2, size = 1, alpha = 0.6) +
  scale_color_gradient(low  = "blue", high = "red", name = "Feature\nvalue") +
  labs(
    #title = "SHAP Beeswarm of Lipids Covering 80% Importance",
    x     = "SHAP value",
    y     = NULL
  ) +
  nature_theme +
  theme(
    legend.position = "right",
    #plot.title      = element_text(hjust = 0.5),
    #axis.text.y     = element_text(size = 8)
  )


# 4) combine side by side, giving more room on the left
quartz()
(p_shap_bee + p_shap_bar) +
  plot_layout(ncol = 2, widths = c(2, 1))


# Save the plot
x <- (p_shap_bee + p_shap_bar)

# ggsave("Fig5f_SHAP_beeswarm.png", x, width = 24, height = 12, dpi = 300, units = "in", bg = "white")
ggsave("SuppFig_SHAP_StemDiameter.png", x, width = 32, height = 24, dpi = 300, units = "in", bg = "white")

# 5) stitch side‑by‑side
quartz()
(p_shap_bee + p_shap_bar) +
  plot_layout(ncol = 2, widths = c(2,1))





### Cut off for SHAP importance
shap_cum <- shap_rank %>% 
  arrange(desc(MeanAbsSHAP)) %>% 
  mutate(cum_pct = cumsum(MeanAbsSHAP) / sum(MeanAbsSHAP))

# figure out which rank first exceeds 80%
cutoff_idx <- which(shap_cum$cum_pct >= 0.80)[1]

quartz()
cutoff_shap <- ggplot(shap_cum, aes(x = seq_along(MeanAbsSHAP), y = cum_pct)) +
  geom_line(size = 1) +
  # vertical line at the cutoff rank
  geom_vline(xintercept = cutoff_idx, linetype = "dashed", color = "red") +
  # horizontal line at 80%
  geom_hline(yintercept = 0.80,    linetype = "dashed", color = "red") +
  # annotation text at the intersection
  annotate(
    "text",
    x    = cutoff_idx + 5,      # nudge to the right
    y    = 0.80 + 0.03,         # nudge upward
    label = sprintf("rank = %d\n80%% cum. imp.", cutoff_idx),
    color = "red",
    hjust = 0,
    size  = 4
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(
    x     = "Feature rank (by mean |SHAP|)",
    y     = "Cumulative importance",
    title = "Pareto plot of SHAP importance"
  ) +
  nature_theme

#save the plot
ggsave("SuppFig_SHAP_cumulative_importance.png", cutoff_shap, width = 8, height = 6, dpi = 300, units = "in", bg = "white")
