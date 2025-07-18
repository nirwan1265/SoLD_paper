library(vroom)
library(dplyr)
library(caret)
library(mlr)
library(tuneRanger)
rm(list=ls())

## GET THE LIPID, Log10 transformed
lipids <- vroom("data/SPATS_fitted/log10_median_centered/control_all_lipids_fitted_phenotype_log_centered.csv")
colnames(lipids)[1] <- "Line"
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


### GET THE PHENOTYPES, FLOWERING TIME AND BLOCK AND ROW
pheno <- vroom("data/phenotypes/control_field_phenotypes.csv") %>% dplyr::select(c(1,2))
colnames(pheno)[1] <- "Line"
colnames(pheno)[2] <- "FlowerTime"
# Add this: drop any rows where FlowerTime is NA
keep_rows <- !is.na(pheno$FlowerTime)
pheno    <- pheno[keep_rows, ]
pheno$FlowerTime <- log10(pheno$FlowerTime + 1)  # Log10 transform FlowerTime

## Combine
dat <- lipids %>% inner_join(pheno, by = "Line")
dim(dat)

## GET THE PCA
PC <- vroom("table/PCA_SAP.csv") 
colnames(PC)[1] <- "Line"
dim(PC)
str(PC)


# 1) Align and parse FlowerTime as numeric  
common <- intersect(dat$Line, PC$Line)
dat2   <- dat %>% 
  filter(Line %in% common) %>% 
  arrange(Line) %>% 
  mutate(FlowerTime = as.numeric(FlowerTime))
PC2    <- PC  %>% 
  filter(Line %in% common) %>% 
  arrange(Line) %>% 
  select(-Line)           # ✂ DROP the ID column


# 2) Build a pure‑numeric data.frame for residualization
df_resid <- data.frame(
  FlowerTime = dat2$FlowerTime,
  PC2      # columns EV1…EV5 only
)


### Population structure adjustment
## Residualise both the phenotype and the lipids on PCs so RF doesn’t learn population structure instead of biology:
adj_formula <- as.formula(paste("~", paste(colnames(PC2), collapse = " + ")))

# 3) Residualize FlowerTime on PCs  
rFT <- resid( lm(FlowerTime ~ ., data = df_resid) )

# 4) Sanity check  
length(rFT)       # should equal nrow(df_resid)
summary(rFT)      # non‑zero residuals with mean ≈0



# ─────────────────────────────────────────────────────────────────────────────
# ASSUMING you’ve already done:
#  • aligned ‘dat’ and ‘PC’ on the same set of Lines
#  • parsed FlowerTime to numeric in dat2
#  • created PC2 by filtering+sorting PC and then dropping Line
#  • created dat2 by filtering+sorting dat and converting FlowerTime
# ─────────────────────────────────────────────────────────────────────────────

# 1) Extract the lipid matrix (numeric only: drop Line & FlowerTime)
lipid_mat <- dat2 %>%
  select(-Line, -FlowerTime) %>%
  as.matrix()

# 2) Build the design matrix X from your PCs only (no Line column)
#adj_formula still = "~ EV1 + EV2 + EV3 + EV4 + EV5"
X <- model.matrix(adj_formula, data = PC2)

# 3) Compute the projection matrix P
P <- X %*% solve(crossprod(X)) %*% t(X)

# 4) Residualize each lipid:  
#    lipids_adj_mat[i,j] = lipid_mat[i,j] - (P %*% lipid_mat)[i,j]
lipids_adj_mat <- lipid_mat - (P %*% lipid_mat)

# 5) Convert back to a tibble and restore column names
lipids_adj <- as_tibble(lipids_adj_mat)
colnames(lipids_adj) <- colnames(lipid_mat)

# 6) (Optional) re‑attach the Line IDs so you can keep track
lipids_adj <- bind_cols(Line = dat2$Line, lipids_adj)

# Now lipids_adj has the same number of rows as dat2 (no ID column was used
# in computing the projection), and each column is the PC‑residualized lipid.




# ─────────────────────────────────────────────────────────────────────────────
# Assume you already have:
#  • lipids_adj: a data.frame or matrix of PC‑residualized lipids (rows = lines)
#  • rFT       : a numeric vector of PC‑residualized flowering times (same order)
#  • PC2      : a data.frame of numeric PCs (no Line column), same rows/order as lipids_adj
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# 1) Stratified train/test split on PCs
# ─────────────────────────────────────────────────────────────────────────────
set.seed(42)
lipids_num <- lipids_adj %>% select(-Line)
clusters  <- kmeans(PC2, centers = 6)$cluster
train_idx <- createDataPartition(clusters, p = 0.8, list = FALSE)

train_X <- lipids_num[train_idx, ]
train_y <- rFT      [train_idx]

test_X  <- lipids_num[-train_idx, ]
test_y  <- rFT      [-train_idx]

# ─────────────────────────────────────────────────────────────────────────────
# 2) Create an mlr regression task for tuning
# ─────────────────────────────────────────────────────────────────────────────
train_df   <- data.frame(train_X, FlowerTime = train_y)
regr_task  <- makeRegrTask(data = train_df, target = "FlowerTime")

# ─────────────────────────────────────────────────────────────────────────────
# 3) Hyperparameter tuning with tuneRanger
#    - up to 10 threads
#    - optimize RMSE
# ─────────────────────────────────────────────────────────────────────────────
# rm(rmse) if you want to re run. 
#rm(rmse) 
tune_res <- tuneRanger(
  task            = regr_task,
  measure         = list(rmse),
  num.trees       = 1000,
  tune.parameters = c("mtry", "min.node.size", "sample.fraction"),
  num.threads     = 10
)

# View recommended hyperparameters:
print(tune_res$recommended.pars)
# $mtry
# $min.node.size
# $sample.fraction

# ─────────────────────────────────────────────────────────────────────────────
# 4) Fit the final RF model with tuned parameters
# ─────────────────────────────────────────────────────────────────────────────
rf_model <- ranger(
  x               = train_X,
  y               = train_y,
  num.trees       = 1000,
  mtry            = tune_res$recommended.pars$mtry,
  #mtry            = 5,
  min.node.size   = tune_res$recommended.pars$min.node.size,
  #min.node.size   = 10,
  sample.fraction = tune_res$recommended.pars$sample.fraction,
  importance      = "none",
  num.threads     = 10,
  seed            = 42
)

# ─────────────────────────────────────────────────────────────────────────────
# 5) Evaluate on the test set
# ─────────────────────────────────────────────────────────────────────────────
preds <- predict(rf_model, data = test_X)$predictions
rmse  <- sqrt(mean((preds - test_y)^2))
print(paste0("Test RMSE: ", round(rmse, 4)))


# 1) Make sure you’re using exactly the same test set for both:
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





# 1) Install & load
library(treeshap)
library(dplyr)

# 2) Compute TreeSHAP values
# rf_model is your trained ranger object
# lipids_adj is your PC‐residualized lipid tibble with a “Line” column

# Extract feature matrix (no Line column)
X <- lipids_adj %>% select(-Line) %>% as.data.frame()

# 1) Convert your ranger forest into a "unified" model
um <- unify(rf_model,X)

# treeshap wants an “unwrapped” ranger forest + the same data.frame
ts <- treeshap(um, X)

# ts$shaps is an (n_samples × n_features) matrix of exact SHAP values

# 3) Global ranking: mean absolute SHAP per lipid
mean_abs <- colMeans(abs(ts$shaps))
shap_rank <- tibble(
  Lipid       = names(mean_abs),
  MeanAbsSHAP = mean_abs
) %>%
  arrange(desc(MeanAbsSHAP))

# View top 20
top20 <- shap_rank %>% slice_head(n = 20)
print(top20)

# 4) (Optional) Convert to shapviz and plot
library(shapviz)
sv <- shapviz(ts, X = X)   # TreeSHAP under the hood

library(dplyr)
library(tidyr)
library(ggplot2)

library(dplyr)
library(tidyr)
library(ggplot2)

# ts is your treeshap result, ts$shaps is an (n_samples × n_features) matrix
shap_mat <- ts$shaps

# 1) Build a global‐importance table by mean(|SHAP|)
global_rank <- tibble(
  Feature = colnames(shap_mat),
  MeanAbs = colMeans(abs(shap_mat))
) %>% 
  arrange(desc(MeanAbs))

# 2) Grab the top 20
top50 <- global_rank %>% slice_head(n = 50)
print(top50, n = 50)

#write.csv(top50, "table/top50_lipid_SHAP_floweringtime.csv", row.names = FALSE)
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




# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: Residual distributions (raw vs PC‑residualized FT and example lipid)
# ─────────────────────────────────────────────────────────────────────────────
library(ggplot2)

# A: Flowering time
df_ft <- tibble(
  raw   = dat2$FlowerTime,
  resid = rFT
)
p1 <- ggplot(df_ft, aes(x = raw)) +
  geom_density(fill = "gray80") +
  labs(title = "Raw Flowering Time", x = "log10(FT+1)", y = "Density") +
  theme_minimal()
p2 <- ggplot(df_ft, aes(x = resid)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  labs(title = "PC‑residualized Flowering Time", x = "Residual FT", y = "Density") +
  theme_minimal()

# B: Example lipid (first column of lipids_adj after Line)
lipid_name <- colnames(lipids_adj)[2]
df_lip <- tibble(
  raw   = lipid_mat[,1],
  resid = lipids_adj[[lipid_name]]
)
p3 <- ggplot(df_lip, aes(x = raw)) +
  geom_density(fill = "gray80") +
  labs(title = paste0("Raw ", lipid_name), x = lipid_name, y = "Density") +
  theme_minimal()
p4 <- ggplot(df_lip, aes(x = resid)) +
  geom_density(fill = "tomato", alpha = 0.6) +
  labs(title = paste0("PC‑residualized ", lipid_name), x = "Residual abundance", y = "Density") +
  theme_minimal()

# Arrange and save
library(patchwork)
quartz()
(p1 + p2) / (p3 + p4)

ggsave("fig2_residual_distributions.png", width = 8, height = 8, dpi = 300)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3: Observed vs Predicted FT
# ─────────────────────────────────────────────────────────────────────────────
df_pred <- tibble(obs = test_y, pred = preds)
p_obs_pred <- ggplot(df_pred, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  labs(title = "Observed vs Predicted Flowering Time",
       x = "Observed residual FT",
       y = "RF predicted residual FT") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = sprintf("RMSE = %.3f\nR² = %.3f", metrics$RMSE, metrics$R2)) +
  theme_minimal()
p_obs_pred
ggsave("fig3_obs_vs_pred.png", width = 6, height = 5, dpi = 300)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4: OOB error vs number of trees
# ─────────────────────────────────────────────────────────────────────────────
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
       x = "Number of Trees", y = "OOB RMSE (days)") +
  theme_minimal()
p_oob
ggsave("fig4_oob_vs_trees.png", width = 6, height = 4, dpi = 300)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 5: Hyperparameter tuning surface (mtry × min.node.size)
# ─────────────────────────────────────────────────────────────────────────────
tuning_df <- tune_res$results
p_tune <- ggplot(tuning_df, aes(x = factor(mtry), y = factor(min.node.size), fill = rmse)) +
  geom_tile() +
  scale_fill_viridis_c(name = "OOB RMSE") +
  labs(title = "Hyperparameter Tuning Landscape",
       x = "mtry", y = "min.node.size") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_tune
ggsave("fig5_tuning_surface.png", width = 6, height = 5, dpi = 300)


# ─────────────────────────────────────────────────────────────────────────────
# Figures 6 & 7: Global SHAP importance & Beeswarm for top 20
# ─────────────────────────────────────────────────────────────────────────────
# Global bar chart (top 20)
top20 <- shap_rank %>% slice_head(n = 20)
p_shap_bar <- ggplot(top20, aes(x = reorder(Lipid, MeanAbsSHAP), y = MeanAbsSHAP)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 Lipids by Mean |SHAP|",
       x = NULL, y = "Mean |SHAP| (days)") +
  theme_minimal()
p_shap_bar
ggsave("fig6_shap_bar.png", width = 6, height = 5, dpi = 300)

# Beeswarm
shap_long <- as_tibble(ts$shaps) %>%
  mutate(Sample = row_number()) %>%
  pivot_longer(-Sample, names_to = "Lipid", values_to = "SHAP") %>%
  filter(Lipid %in% top20$Lipid)
p_shap_bee <- ggplot(shap_long, aes(x = SHAP, y = reorder(Lipid, -SHAP))) +
  geom_jitter(height = 0.2, alpha = 0.4, size = 1) +
  labs(title = "Beeswarm of Top 20 Lipid SHAP Values",
       x = "SHAP value (days)", y = NULL) +
  theme_minimal()
p_shap_bee
ggsave("fig7_shap_beeswarm.png", width = 6, height = 6, dpi = 300)
