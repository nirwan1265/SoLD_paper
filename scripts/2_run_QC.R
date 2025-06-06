library(tidyverse)
library(reshape2)    # for melt()
library(FactoMineR)  # PCA
library(factoextra)  # nicer PCA plots
library(pheatmap)    # heat-map



# Control
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/")
A <- read.csv("data/SetA_lipid_FLO2019Control.csv")
B <- read.csv("data/SetB_lipid_FLO2022_lowP.csv")

# Retention time > 1 min
A <- A[which(A$row.retention.time >= 1), ]


colnames(A)


library(tidyverse)
library(stringr)

## ------------------------------------------------------------------ ##
## 1.  isolate the analytical columns (everything after col 13)       ##
## ------------------------------------------------------------------ ##
A_peaks <- A[ , 14:ncol(A)]

## ------------------------------------------------------------------ ##
## 2.  canonical column-name cleaner                                  ##
##      (works for QC, InjBL, ISTD & S1 injections)                   ##
## ------------------------------------------------------------------ ##
clean_name <- function(x){
  x %>%                                    # original long name
    str_replace("^.*?\\.\\.", "") %>%      # drop all up to last “..”
    str_replace("\\.mzML.*$", "")          # drop   .mzML.Peak.area suffix
}

new_names <- clean_name(names(A_peaks))

## optional: add PI to S1 columns  ----------------------------------- ##
new_names <- new_names %>% map_chr(function(nm){
  if(str_starts(nm, "S1_Run")){
    pi <- str_extract(nm, "PI\\d+")
    run<- str_extract(nm, "S1_Run\\d+")
    paste0(run, "_", pi)
  } else nm
})

colnames(A_peaks) <- new_names

## ------------------------------------------------------------------ ##
## 3.  tag each column & pull the *actual* Run number                 ##
## ------------------------------------------------------------------ ##
col_type <- function(x){
  case_when(
    str_starts(x,"QC_")    ~ "QC",
    str_starts(x,"InjBL")  ~ "Blank",
    str_starts(x,"ISTD")   ~ "ISTD",
    str_starts(x,"S1_Run") ~ "Sample",
    TRUE                   ~ "QC")
}

get_run_no <- function(x){
  as.numeric(str_extract(x, "(?<=Run)\\d+"))   # grabs the digits after “Run”
}

## ------------------------------------------------------------------ ##
## 4.  build the numeric matrix + metadata frame, then RE-ORDER       ##
## ------------------------------------------------------------------ ##
mat <- as.matrix(A_peaks);  mode(mat) <- "numeric"

qc_meta <- tibble(
  Injection = colnames(mat),
  Type      = col_type(Injection),
  Run       = get_run_no(Injection)           # ← numeric injection order
) |>
  arrange(Run)                                # sort by true queue order

# re-order the columns of the peak matrix the same way
mat <- mat[ , qc_meta$Injection]

# quick sanity check
table(qc_meta$Type)              # how many QCs, blanks, ISTDs, samples?

# total-ion-current dataframe
tic_df <- qc_meta |>
  mutate(TIC = colSums(mat, na.rm = TRUE))

## ------------------------------------------------------------------ ##
## 5.  TIC plot                                                       ##
## ------------------------------------------------------------------ ##
library(ggplot2)

quartz()
ggplot(tic_df, aes(Run, TIC, colour = Type)) +
  geom_line() +
  geom_point(size = 2) +
  scale_colour_manual(values = c(
    Sample = "steelblue",
    QC     = "tomato",
    Blank  = "grey40",
    ISTD   = "black")) +
  theme_bw(base_size = 14) +
  labs(x = "Injection order (Run #)",
       y = "Total ion current")





library(tidyverse)
library(FactoMineR)   # PCA core
library(factoextra)   # pretty ggplot helpers

## ------------------------------------------------------------------ ##
## 1.  prepare matrices                                               ##
## ------------------------------------------------------------------ ##
#  – mat  : numeric feature × injection matrix   (built earlier)
#  – qc_meta : Injection, Type, Order            (built earlier)

## optional: remove features that are zero everywhere (speeds PCA)
mat_filt <- mat[rowSums(mat) > 0, ]
str(mat)

## centre-scale & run PCA (keep, say, first 10 PCs)
pca_res <- PCA(t(mat_filt),            # transpose → injections as rows
               scale.unit = TRUE,
               graph      = FALSE,
               ncp        = 10)        # number of components

## ------------------------------------------------------------------ ##
## 2.  scatter of PC1 vs PC2 with QC colours                          ##
## ------------------------------------------------------------------ ##
p_df <- cbind(qc_meta, pca_res$ind$coord[, 1:2]) %>%        # PC coords
  rename(PC1 = Dim.1, PC2 = Dim.2)

quartz()
ggplot(p_df, aes(PC1, PC2, colour = Type))+
  geom_point(size = 3, alpha = 0.8)+
  scale_colour_manual(values = c(QC = "tomato",
                                 Sample = "steelblue",
                                 Blank = "grey40",
                                 ISTD = "black"))+
  theme_bw(14)+
  labs(title = "PCA of raw peak areas",
       subtitle = sprintf("Explained variance: PC1 %.1f %%  •  PC2 %.1f %%",
                          pca_res$eig[1,2], pca_res$eig[2,2]),
       x = "PC-1", y = "PC-2")









####  For lowinput

# 0. read & pre-filter Set B exactly like Set A …
B <- read.csv("data/SetB_lipid_FLO2022_lowP.csv")


# Retention time > 1 min
B <- B[which(B$row.retention.time >= 1), ]


B_peaks  <- B[ , 14:ncol(B)]

## ------------------------------------------------------------------ ##
## 2.  column-name cleaner *for Set B*                                 ##
## ------------------------------------------------------------------ ##
clean_B_name <- function(x){
  x %>% 
    str_replace("^.*?\\.FLO\\d+\\.", "") %>%   # keep everything *after* “.FLO22.”
    str_replace("\\.mzML.*$", "")              # drop “.mzML.Peak.area”
}

new_names <- clean_B_name(colnames(B_peaks))

## ---- attach PI tag to any S{digit}_Run ---- ##
new_names <- map_chr(new_names, function(nm){
  if( str_detect(nm, "^S\\d+_Run") ){
    pi_tag <- str_extract(nm, "PI\\d+")
    if(is.na(pi_tag)) pi_tag <- "PI0000"
    paste0(str_extract(nm, "S\\d+_Run\\d+"), "_", pi_tag)
  } else nm
})

colnames(B_peaks) <- new_names

## ------------------------------------------------------------------ ##
## 3.  helpers: injection *type*  &  run number                        ##
## ------------------------------------------------------------------ ##
col_type_B <- function(x){
  case_when(
    str_starts(x, "QC_")            ~ "QC",
    str_starts(x, "InjBL")          ~ "Blank",
    str_starts(x, "ISTD")           ~ "ISTD",
    str_detect(x, "^S\\d+_Run")     ~ "Sample",
    TRUE                            ~ "QC"
  )
}


# These files are no converted but they are QCs
#[1] "X20230621_JGI_RRA_508500_Sorgh_final.SetB_EXP120B_C18.LipidV7_USDAY63675_POS_MS2_0_QC_Pre_Rg132to1500.CE102040..QC_Run7"  
#[777] "X20230621_JGI_RRA_508500_Sorgh_final.SetB_EXP120B_C18.LipidV7_USDAY63675_POS_MS2_0_QC_Post_Rg132to1500.CE102040..QC_Run790"


get_run_no <- function(x){
  as.numeric(str_extract(x, "(?<=Run)\\d+"))
}

## ------------------------------------------------------------------ ##
## 4.  build matrix + metadata frame  (sorted by Run)                  ##
## ------------------------------------------------------------------ ##
mat_B          <- as.matrix(B_peaks);  mode(mat_B) <- "numeric"

# remove columns 778, 392, and 398 from mat_B
mat_B <- mat_B[ , -c(778, 392, 398)]  # remove columns with no data


meta_B <- tibble(
  Injection = colnames(mat_B),
  Type      = col_type_B(Injection),
  Run       = get_run_no(Injection)
) %>% arrange(Run)

mat_B <- mat_B[ , meta_B$Injection]           # reorder columns

message("Counts per category:")
print(table(meta_B$Type))

## ------------------------------------------------------------------ ##
## 5.  TIC vs. injection order plot                                    ##
## ------------------------------------------------------------------ ##
tic_df_B <- meta_B %>% 
  mutate(TIC = colSums(mat_B, na.rm = TRUE))

quartz()
ggplot(tic_df_B, aes(Run, TIC, colour = Type)) +
  geom_line() +
  geom_point(size = 2) +
  scale_colour_manual(values = c(
    Sample = "steelblue",
    QC     = "tomato",
    Blank  = "grey40",
    ISTD   = "black",
    Other  = "purple"
  )) +
  theme_bw(14) +
  labs(x = "Injection order (Run #)",
       y = "Total ion current",
       title = "Set-B TIC trace")

# ————————————————————————————————————————————————————————————————————— #
# 1) assume you have already built “mat_B” (numeric, 4418×778 with many NA’s)
# and “qc_meta” (a 778×3 tibble with columns Injection, Type, Run)
# ————————————————————————————————————————————————————————————————————— #

library(missMDA)       # for imputePCA()
library(FactoMineR)    # for PCA()
library(dplyr)
library(ggplot2)

# (a) drop any lipid that is zero everywhere (optional, just speeds up things)
mat_B_nonzero <- mat_B[rowSums(mat_B, na.rm=TRUE) > 0, ]


# # 1) Identify which columns are true “Sample” injections:
# is_sample_col <- grepl("^S1_Run", colnames(mat_B_nonzero))
# 
# # 2) Loop over each Sample‐column, find its minimum non‐NA, and fill NAs with min/3
# for (j in which(is_sample_col)) {
#   col_vals   <- mat_B_nonzero[, j]
#   min_non_na <- min(col_vals, na.rm = TRUE)
#   mat_B_nonzero[is.na(col_vals), j] <- min_non_na / 3
# }
# 
# is_sample_col <- grepl("^QC", colnames(mat_B_nonzero))
# for (j in which(is_sample_col)) {
#   col_vals   <- mat_B_nonzero[, j]
#   min_non_na <- min(col_vals, na.rm = TRUE)
#   mat_B_nonzero[is.na(col_vals), j] <- min_non_na / 3
# }
# is_sample_col <- grepl("^Inj", colnames(mat_B_nonzero))
# for (j in which(is_sample_col)) {
#   col_vals   <- mat_B_nonzero[, j]
#   min_non_na <- min(col_vals, na.rm = TRUE)
#   mat_B_nonzero[is.na(col_vals), j] <- min_non_na / 3
# }

# Now mat_B_nonzero has no NAs in any “S1_Run” column. You can proceed with PCA directly:
library(FactoMineR)
mat_for_pca <- mat_B_nonzero[rowSums(mat_B_nonzero) > 0, ]     # drop any lipids still all 0
pca_res <- PCA(t(mat_for_pca), scale.unit = TRUE, graph = FALSE, ncp = 10)

# And rebuild the PC1–PC2 plot as before:
p_df <- meta_B %>%
  arrange(Run) %>%
  bind_cols(as_tibble(pca_res$ind$coord[,1:2]) %>% rename(PC1 = Dim.1, PC2 = Dim.2))

library(ggplot2)
quartz()
ggplot(p_df, aes(x = PC1, y = PC2, colour = Type)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = c(
    Sample = "steelblue", QC = "tomato", Blank = "grey40", ISTD = "black", Other = "purple"
  )) +
  theme_bw(base_size = 14) +
  labs(
    title    = "Set B PCA (NAs in Sample columns → 1/3 min→ imputed)",
    subtitle = sprintf("Explained var: PC1 = %.1f%%   PC2 = %.1f%%",
                       pca_res$eig[1,2], pca_res$eig[2,2]),
    x = "PC 1", y = "PC 2"
  )
