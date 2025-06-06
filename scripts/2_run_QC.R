################################################################################
############ Load the libraries
################################################################################

library(tidyverse)
library(reshape2)    
library(FactoMineR)  
library(factoextra)  
library(pheatmap)
library(viridis)
library(stringr)
library(ggrepel)


################################################################################
############ Loading the Quantitative Peaks Intensities.
################################################################################


# Control
setwd("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SAP_lipids_GWAS/")
A <- read.csv("data/SetA_lipid_FLO2019Control.csv")

# Retention time > 1 min
A <- A[which(A$row.retention.time >= 1), ]


################################################################################
############ isolate the analytical columns (everything after col 13)
################################################################################

A_peaks <- A[ , 14:ncol(A)]


################################################################################
############ canonical column-name cleaner(works for QC, InjBL, ISTD & S1 injections)   
################################################################################

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


################################################################################
############ tag each column & pull the *actual* Run number             
################################################################################

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


################################################################################
############  Build the numeric matrix + metadata frame, then RE-ORDER      
################################################################################

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



################################################################################
############  TIC plot                                                      
################################################################################

quartz()
ggplot(tic_df, aes(Run, TIC, colour = Type, group = Type)) +
  # Geometries - hollow circles with colored outlines
  geom_line(linewidth = 0.6, alpha = 0.5) +
  geom_point(size = 3, shape = 21, fill = "white", stroke = 1.2) +
  
  # Color scheme (optimized for visibility)
  scale_colour_manual(
    values = c(
      Sample = "#4477AA",  # Distinct blue
      QC     = "#EE6677",  # Salmon for QC
      Blank  = "#777777",  # Neutral grey
      ISTD   = "#228833"   # Green for ISTD
    ),
    labels = c(
      Sample = "Sample",
      QC     = "Quality Control",
      Blank  = "Blank",
      ISTD   = "Internal Std."
    )
  ) +
  
  # Highlight max values
  ggrepel::geom_text_repel(
    data = tic_df %>% group_by(Type) %>% slice_max(TIC, n = 1),
    aes(label = paste0("Max: ", round(TIC, 0))),
    size = 3.2,
    nudge_y = 0.05 * max(tic_df$TIC),
    box.padding = 0.3,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  
  # Scales
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    labels = scales::scientific
  ) +
  
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.position = c(0.95, 0.95),  # Top right corner
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "grey70"),
    legend.key = element_blank(),
    legend.spacing = unit(0.2, "cm"),
    panel.grid.major.y = element_line(colour = "grey92")
  ) +
  
  # Labels
  labs(
    x = "Injection Order",
    y = "Total Ion Current (TIC)",
    #title = "Lipidomics Run Quality Control"
  ) +
  
  # Fix legend symbols to match plot points
  guides(
    colour = guide_legend(
      override.aes = list(
        shape = 21,       # Hollow circles
        fill = "white",   # White fill
        stroke = 1,       # Border thickness
        linetype = 0      # Remove line in legend
      )
    )
  )


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

# 1. Compute percent‐variance labels
var_pc1 <- round(pca_res$eig[1, 2], 1)
var_pc2 <- round(pca_res$eig[2, 2], 1)

# 2. Find “outlier” per Type (largest PC1^2 + PC2^2)
outliers <- p_df %>%
  group_by(Type) %>%
  slice_max(order_by = PC1^2 + PC2^2, n = 1) %>%
  ungroup() %>%
  mutate(label = Injection)  # use the Injection name as the label

# 3. Build the PCA scatter with hulls, ellipses, and outlier labels
hulls <- p_df %>%
  group_by(Type) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

quartz()
ggplot(p_df, aes(PC1, PC2, color = Type, fill = Type)) +
  # convex hull areas (transparent fill, no legend)
  geom_polygon(data = hulls, aes(PC1, PC2), alpha = 0.1, show.legend = FALSE) +
  # points (hollow circles with colored border and fill)
  geom_point(size = 4, shape = 21, stroke = 1, alpha = 0.9) +
  # 80% confidence ellipses
  stat_ellipse(level = 0.8, linewidth = 0.7, linetype = "dashed") +
  
  # manual colour + fill scale
  scale_color_manual(
    name = "Sample Type",
    values = c(
      Sample = "#4477AA",   # blue
      QC     = "#EE6677",   # salmon
      Blank  = "#777777",   # grey
      ISTD   = "#228833"    # green
    ),
    aesthetics = c("color", "fill")
  ) +
  
  # axes labels and title
  labs(
    x       = sprintf("PC1 (%.1f%% variance)", var_pc1),
    y       = sprintf("PC2 (%.1f%% variance)", var_pc2),
    #title   = "Lipidomics PCA: Sample Type Separation",
    caption = "80% confidence ellipses shown"
  ) +
  
  # minimal theme with customized text
  theme_minimal(base_size = 12) +
  theme(
    text            = element_text(family = "Arial"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    axis.line       = element_line(color = "black"),
    axis.ticks      = element_line(color = "black"),
    axis.title      = element_text(face = "bold", size = 12),
    axis.text       = element_text(color = "black"),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "grey70"),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.caption    = element_text(size = 9, color = "grey40")
  ) +
  
  # lock aspect ratio so one unit in PC1 equals one unit in PC2
  coord_fixed(ratio = 1) +
  
  # label each outlier with its Injection name
  geom_text_repel(
    data          = outliers,
    mapping       = aes(PC1, PC2, label = label),
    size          = 3,
    box.padding   = 0.5,
    show.legend   = FALSE
  )





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


qc_meta_B <- tibble(
  Injection = colnames(mat_B),
  Type      = col_type_B(Injection),
  Run       = get_run_no(Injection)
) %>% arrange(Run)

mat_B <- mat_B[ , qc_meta_B$Injection]           # reorder columns

message("Counts per category:")
print(table(meta_B$Type))

## ------------------------------------------------------------------ ##
## 5.  TIC vs. injection order plot                                    ##
## ------------------------------------------------------------------ ##
tic_df_B <- qc_meta_B %>% 
  mutate(TIC = colSums(mat_B, na.rm = TRUE))



quartz()
ggplot(tic_df_B, aes(Run, TIC, colour = Type, group = Type)) +
  # Geometries - hollow circles with colored outlines
  geom_line(linewidth = 0.6, alpha = 0.5) +
  geom_point(size = 3, shape = 21, fill = "white", stroke = 1.2) +
  
  # Color scheme (optimized for visibility)
  scale_colour_manual(
    values = c(
      Sample = "#4477AA",  # Distinct blue
      QC     = "#EE6677",  # Salmon for QC
      Blank  = "#777777",  # Neutral grey
      ISTD   = "#228833"   # Green for ISTD
    ),
    labels = c(
      Sample = "Sample",
      QC     = "Quality Control",
      Blank  = "Blank",
      ISTD   = "Internal Std."
    )
  ) +
  
  # Highlight max values
  ggrepel::geom_text_repel(
    data = tic_df_B %>% group_by(Type) %>% slice_max(TIC, n = 1),
    aes(label = paste0("Max: ", round(TIC, 0))),
    size = 3.2,
    nudge_y = 0.05 * max(tic_df_B$TIC),
    box.padding = 0.3,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  
  # Scales
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    labels = scales::scientific
  ) +
  
  # Theme
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.position = c(0.95, 0.95),  # Top right corner
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "grey70"),
    legend.key = element_blank(),
    legend.spacing = unit(0.2, "cm"),
    panel.grid.major.y = element_line(colour = "grey92")
  ) +
  
  # Labels
  labs(
    x = "Injection Order",
    y = "Total Ion Current (TIC)",
    #title = "Lipidomics Run Quality Control"
  ) +
  
  # Fix legend symbols to match plot points
  guides(
    colour = guide_legend(
      override.aes = list(
        shape = 21,       # Hollow circles
        fill = "white",   # White fill
        stroke = 1,       # Border thickness
        linetype = 0      # Remove line in legend
      )
    )
  )



# ————————————————————————————————————————————————————————————————————— #
# 1) assume you have already built “mat_B” (numeric, 4418×778 with many NA’s)
# and “qc_meta” (a 778×3 tibble with columns Injection, Type, Run)
# ————————————————————————————————————————————————————————————————————— #

# (a) drop any lipid that is zero everywhere (optional, just speeds up things)
mat_B_nonzero <- mat_B[rowSums(mat_B, na.rm=TRUE) > 0, ]


# Now mat_B_nonzero has no NAs in any “S1_Run” column. You can proceed with PCA directly:
mat_for_pca <- mat_B_nonzero[rowSums(mat_B_nonzero) > 0, ]     # drop any lipids still all 0
pca_res <- PCA(t(mat_for_pca), scale.unit = TRUE, graph = FALSE, ncp = 10)

# And rebuild the PC1–PC2 plot as before:
p_df <- meta_B %>%
  arrange(Run) %>%
  bind_cols(as_tibble(pca_res$ind$coord[,1:2]) %>% rename(PC1 = Dim.1, PC2 = Dim.2))

# 1. Compute percent‐variance labels
var_pc1 <- round(pca_res$eig[1, 2], 1)
var_pc2 <- round(pca_res$eig[2, 2], 1)

# 2. Find “outlier” per Type (largest PC1^2 + PC2^2)
outliers <- p_df %>%
  group_by(Type) %>%
  slice_max(order_by = PC1^2 + PC2^2, n = 1) %>%
  ungroup() %>%
  mutate(label = Injection)  # use the Injection name as the label

# 3. Build the PCA scatter with hulls, ellipses, and outlier labels
hulls <- p_df %>%
  group_by(Type) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

quartz()
ggplot(p_df, aes(PC1, PC2, color = Type, fill = Type)) +
  # convex hull areas (transparent fill, no legend)
  geom_polygon(data = hulls, aes(PC1, PC2), alpha = 0.1, show.legend = FALSE) +
  # points (hollow circles with colored border and fill)
  geom_point(size = 4, shape = 21, stroke = 1, alpha = 0.9) +
  # 80% confidence ellipses
  stat_ellipse(level = 0.8, linewidth = 0.7, linetype = "dashed") +
  
  # manual colour + fill scale
  scale_color_manual(
    name = "Sample Type",
    values = c(
      Sample = "#4477AA",   # blue
      QC     = "#EE6677",   # salmon
      Blank  = "#777777",   # grey
      ISTD   = "#228833"    # green
    ),
    aesthetics = c("color", "fill")
  ) +
  
  # axes labels and title
  labs(
    x       = sprintf("PC1 (%.1f%% variance)", var_pc1),
    y       = sprintf("PC2 (%.1f%% variance)", var_pc2),
    #title   = "Lipidomics PCA: Sample Type Separation",
    caption = "80% confidence ellipses shown"
  ) +
  
  # minimal theme with customized text
  theme_minimal(base_size = 12) +
  theme(
    text            = element_text(family = "Arial"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey92"),
    panel.grid.minor = element_blank(),
    axis.line       = element_line(color = "black"),
    axis.ticks      = element_line(color = "black"),
    axis.title      = element_text(face = "bold", size = 12),
    axis.text       = element_text(color = "black"),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "grey70"),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.caption    = element_text(size = 9, color = "grey40")
  ) +
  
  # lock aspect ratio so one unit in PC1 equals one unit in PC2
  coord_fixed(ratio = 1) +
  
  # label each outlier with its Injection name
  geom_text_repel(
    data          = outliers,
    mapping       = aes(PC1, PC2, label = label),
    size          = 3,
    box.padding   = 0.5,
    show.legend   = FALSE
  )
