#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
#### SOrghum Lipidomics Database (SOLD)
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


###############################################################################
## Supplementary Figures – Overlap of Control vs Low‑input lipid lists
## Three panels:
##   S1 – All detected lipids               (Venn)
##   S2 – Traditional classes only          (Venn)
##   S3 – Curated “Class” annotation set    (table + bar *optional*)
###############################################################################

# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------
library(vroom)      # fast CSV import
library(dplyr)      # data wrangling
library(tidyr)      # pivoting
library(ggplot2)    # plotting
library(stringr)    # string manipulation
library(forcats)    # factor reordering
library(tidyverse)
library(ggpubr)      # For stat_compare_means
library(scales) 
library(ggVennDiagram)   # nice, ggplot‑based 2‑set Venn
library(readr) 
library(sf)
library(ggVennDiagram) # For Venn diagrams)
library(viridis)

# ------------------------------------------------------------------------------
# 2. Read in the data
# ------------------------------------------------------------------------------
control  <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Control_all_lipids_final_non_normalized.csv")

lowinput <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv")


 # ------------------------------------------------------------------------------
# 3."Traditional” lipid classes and extractor
# ------------------------------------------------------------------------------
valid_classes <- c("TG","DG","MG",
                   "PC","PE","PG","PI",
                   "LPC","LPE",
                   "DGDG","MGDG",
                   "Cer","SM","FA","DGDG","SQDG","AEG")

class_pattern <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

extract_class <- function(x) str_extract(x, class_pattern)


# ------------------------------------------------------------------------------
# 4. Define aesthetics for plotting
# ------------------------------------------------------------------------------
nature_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, 
                              margin = margin(b = 10)),  # Title styling
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(color = "black"),  # Ensure axis labels are visible
    axis.text.y = element_text(face = "bold", hjust = 1),  # Bold class names
    axis.line = element_line(color = "black"),  # Add axis lines
    panel.grid = element_blank(),  # Remove all grid lines
    legend.position = "top",
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 10),
    plot.margin = margin(15, 15, 15, 15)  # Balanced white space
  )

# Color scheme (colorblind-friendly)
#condition_colors <- c(Control = "#4E79A7", LowInput = "#E15759") 
#condition_colors <- c(Control = "white", LowInput = "black") 
group_colors <- viridis(2)

condition_colors <- c(Control = "#440154FF", LowInput = "#FDE725FF") 


# -------------------------------------------------------------------------------
# 5. Supplementary Figure 1 – Lipid class counts (Working Version)
# -------------------------------------------------------------------------------

# Function to get class counts
get_class_counts <- function(df, condition_label) {
  df %>%
    dplyr::select(-Compound_Name) %>%    # drop sample IDs
    names() %>%                          # get lipid names
    extract_class() %>%                  # extract class prefix
    na.omit() %>%                        # drop non-matching
    as_tibble() %>%
    dplyr::count(value, name = "Count") %>% 
    dplyr::rename(Lipid_Class = value) %>%
    dplyr::mutate(Condition = condition_label)
}

# Get the class counts
counts_ctrl   <- get_class_counts(control,  "Control")
counts_low    <- get_class_counts(lowinput, "LowInput")
class_counts  <- bind_rows(counts_ctrl, counts_low)
print(class_counts, n = 30)

# How far the text sits 
offset <- max(class_counts$Count) * 0.04   # 3 % of the axis length

# Plot fir 1A
fig1a <- ggplot(class_counts,
                aes(x = fct_reorder(Lipid_Class, Count, .desc = TRUE),
                    y = Count,
                    fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.6,
           colour   = "black",
           linewidth = 0.25) +
  geom_text(aes(label = Count,
                y     = Count + max(Count) * 0.04),
            position = position_dodge(width = 0.8),
            hjust     = 0,
            size      = 3.3) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = condition_colors ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL,
       y = "Number of lipids detected") +
  nature_theme +                           # your base theme
  theme(
    # grid lines
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "#d0d0d0", linewidth = 0.2),
    # legend tweaks
    legend.position     = c(0.78, 0.93),   # inside top‑right
    legend.direction    = "horizontal",
    legend.background   = element_rect(fill = "white",
                                       colour = "black", linewidth = 0.2),
    legend.key.size     = unit(0.25, "cm"),
    legend.title        = element_blank(), # Nature rarely titles legends
    # axis text
    axis.text.y         = element_text(face = "bold")
  )


quartz()
print(fig1a)

# Save the plot
ggsave("SuppFig_1A_Lipid_Counts.png",
       plot = fig1a,
       width = 6, height = 6, dpi = 300,
       units = "in", bg = "white")  # white background for publication quality



# -------------------------------------------------------------------------------
# 6.Supplementary Figure 2 – Lipid‑class diversity using curated "Class" annotations
# -------------------------------------------------------------------------------

# Lipid class
lipid_class_info <- vroom::vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/data/lipid_class.csv", show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Class))
unique(lipid_class_info$Class)
unique(lipid_class_info$SubClass)

# Select which class you want to work with. Here I took class
lookup_tbl <- lipid_class_info %>%
  dplyr::select(Lipid = Lipids, Class) %>%
  dplyr::mutate(Lipid = str_trim(Lipid))

# Remove duplicates in the lookup table
lookup_tbl <- lookup_tbl %>%
  dplyr::select(Lipid, Class) %>%      # keep only mapping columns
  dplyr::distinct(Lipid, .keep_all = TRUE)   # drop duplicates


# Any lipid not in list will be tagged "Other"
lookup_fun <- function(nms) {
  tibble(Lipid = str_trim(nms)) %>%
    left_join(lookup_tbl, by = "Lipid") %>%
    mutate(Class = replace_na(Class, "Other"))
}

# Count lipids per class per condition
count_classes <- function(df, cond_label) {
  # Extract all lipid names (all column names except the first one)
  lipid_names <- names(df)[-1]
  
  # Map each lipid name to its Class, then group & summarize
  lookup_fun(lipid_names) %>%
    group_by(Class) %>%
    summarise(
      Count  = n(),
      Lipids = paste(Lipid, collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(Condition = cond_label) %>%
    select(Class, Count, Lipids, Condition)
}

counts_ctrl <- count_classes(control,  "Control")
counts_low  <- count_classes(lowinput, "LowInput")

class_counts <- bind_rows(counts_ctrl, counts_low)

# Order bars by total (#Control + #LowInput)
class_counts <- class_counts %>%
  group_by(Class) %>%
  mutate(total = sum(Count)) %>%
  ungroup() %>%
  arrange(desc(total))



# Remove Plasticizer, Herbicide, and Surfactant from  class_counts Class column
class_counts <- class_counts %>%
  filter(!Class %in% c("Plasticizer", "Herbicide", "Surfactant","Functional compound","Functional agent","Organic Compound","Lipid binding protein","Phenethylamine", "Alkaloid-like compound", "Psychoactive compound","Lipid-like compound","Prenol","Dicarboxylic acid","Phenylpropanoid","Flavonoid","Organic compound","Phenol derivative","Other"))
unique(class_counts$Class)
# 6. Plot Figure 1B 
# How far the numbers sit
offset <- max(class_counts$Count) * 0.04   # 3 % of the axis length

fig1b <- ggplot(class_counts,
                aes(x = fct_reorder(Class, total, .desc = TRUE),
                    y = Count,
                    fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.6,
           colour   = "black",
           linewidth = 0.25) +
  # ---------- COUNT LABELS ---------------------------------------------------
geom_text(aes(x = fct_reorder(Class, total, .desc = TRUE),
              y = Count + offset,                # little nudge outside bar
              label = Count),
          position = position_dodge(width = 0.8),
          hjust     = 0,
          size      = 3.3) +
  coord_flip(clip = "off") +                       # allow text past plotting area
  scale_fill_manual(values = condition_colors) +
  labs(x = NULL,
       y = "Number of lipids detected") +
  nature_theme +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "#d0d0d0", linewidth = 0.2),
    legend.position    = c(0.80, 0.93),
    legend.direction   = "horizontal",
    legend.background  = element_rect(fill = "white",
                                      colour = "black", linewidth = 0.2),
    legend.key.size    = unit(0.25, "cm"),
    legend.title       = element_blank(),
    axis.text.y        = element_text(face = "bold")
  )

quartz()
print(fig1b)

# Save the plot
ggsave("SuppFig_1B_Lipid_Class_Counts.png",
       plot = fig1b,
       width = 6, height = 6, dpi = 300,
       units = "in", bg = "white")  # white background for publication quality




# -------------------------------------------------------------------------------
# 7. Supplementary Figure 3 – Common Lipids - Venn Diagram
# -------------------------------------------------------------------------------



# 1) Filter for Control, extract the “Lipids” column, split each comma‐sep string,
#    and then unlist + unique() to get one vector of all unique lipid names
ctrl_names <- class_counts %>%
  filter(Condition == "Control") %>%    # keep only Control rows
  pull(Lipids) %>%                       # extract the “Lipids” column (each element is one comma‐sep string)
  str_split("\\s*,\\s*") %>%             # split each element on commas (allow optional spaces around “,”)
  unlist() %>%                           # flatten into one long character vector
  unique()                               # drop any duplicates

# 2) Do the same for LowInput:
low_names <- class_counts %>%
  filter(Condition == "LowInput") %>%
  pull(Lipids) %>%
  str_split("\\s*,\\s*") %>%
  unlist() %>%
  unique()

# 2) Prepare your two vectors of lipid‐feature names:
#    (just replace these with your actual character vectors)
ctrl_names  <- ctrl_names   # e.g. c("LipidA","LipidB","LipidC",…)
low_names   <- low_names    # e.g. c("LipidX","LipidY","LipidZ",…)


fill_color <- c("Control" = "#440154FF", "LowInput" = "#FDE725FF") 

# 3) Put them into a named list for ggvenn:
venn_list <- list(
  Control  = ctrl_names,
  LowInput = low_names
)

# 4) Call ggvenn with show_percentage = TRUE, plus custom blues:
# (3) Call ggvenn with an unnamed fill_color vector (first color → Control, second → LowInput):
library(ggvenn)
p <- ggvenn(
  venn_list,
  fill_color      = c("#440154FF", "#FDE725FF"),  # first = Control, second = LowInput
  stroke_size     = 0.8,    # thickness of the black outlines
  set_name_size   = 6,      # font size for “Control” / “LowInput”
  text_size       = 5,      # font size for the region labels (counts/percentages)
  show_percentage = TRUE    # automatically draws “count\n(XX%)” in each region
) +
  # 5) Add a title and increase margins so nothing gets clipped:
  #ggtitle("Lipid Feature Overlap Between Conditions") +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.margin     = margin(20, 20, 20, 20),  # extra space on all sides
    legend.position = "none"
  )

# 6) Print
quartz()
print(p)

# Save the plot
ggsave("SuppFig_1C_Lipid_Overlap_Venn.png",
       plot = p,
       width = 6, height = 6, dpi = 300,
       units = "in", bg = "white")  # white background for publication quality




