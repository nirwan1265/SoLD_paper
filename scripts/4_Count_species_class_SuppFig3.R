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


# Remove PC(17:0) from control and lowinput (internal standard) from the columns
control  <- control %>% dplyr::select(-`PC(17:0)`)
lowinput <- lowinput %>% dplyr::select(-`PC(17:0)`)

#### CHANGE THE NAME TO COMMONNAME
# Lipid class
lipid_class_info <- vroom::vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/data/lipid_class.csv", show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Class))


# Create a named vector of replacements: names = original, values = new names
name_map <- lipid_class_info %>%
  filter(!is.na(CommonName)) %>%
  select(Lipids, CommonName) %>%
  deframe()  # turns it into a named vector: Lipids -> CommonName

# Function to rename lipid columns in any dataset
rename_lipid_columns <- function(df, name_map) {
  original_names <- colnames(df)
  
  renamed <- original_names %>%
    # Replace only if a match exists
    sapply(function(x) if (x %in% names(name_map)) name_map[[x]] else x)
  
  # Apply renamed vector back
  colnames(df) <- renamed
  return(df)
}

# Apply to both control and lowinput
control  <- rename_lipid_columns(control,  name_map)
lowinput <- rename_lipid_columns(lowinput, name_map)

colnames(lowinput)
# ------------------------------------------------------------------------------
# 3."Traditional” lipid classes and extractor
# ------------------------------------------------------------------------------
valid_classes <- c("TG","DG","MG",
                   "PC","PE","PG","PI",
                   "LPC","LPE",
                   "DGDG","MGDG",
                   "Cer","SM","FA","DGDG","SQDG","AEG","GalCer","FA")

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
ggsave("SuppFig_3A_Lipid_Counts.png",
       plot = fig1a,
       width = 6, height = 6, dpi = 300,
       units = "in", bg = "white")  # white background for publication quality









# -------------------------------------------------------------------------------
# 6.Supplementary Figure 2 – Lipid‑class diversity using curated "Class" annotations
# -------------------------------------------------------------------------------

# Lipid class
lipid_class_info <- vroom::vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/data/lipid_class.csv", show_col_types = FALSE) %>%
  dplyr::filter(!is.na(Class))

# Lipid class Traditional lipids
traditional_lipid_classes <- c("Glycerolipid", "Glycerophospholipid", "Glycoglycerolipid",
                         "Sphingolipid", "Sterol", "Betaine lipid", "Fatty acid and derivative")
lipid_class_traditional <- lipid_class_info %>% 
  dplyr::filter(Class %in% traditional_lipid_classes)

# Lipid class Traditional lipids
nontraditional_lipid_class <- c("N-acylethanolamine","Terpenoid","Prenol","Tetrapyrrole","Vitamin")

lipid_class_nontraditional <- lipid_class_info %>% 
  dplyr::filter(Class %in% nontraditional_lipid_class)


# Use this updated lipid_class_info
lipid_class_info <- lipid_class_info %>%
  mutate(final_name = ifelse(is.na(CommonName), Lipids, CommonName))

# Redefine traditional and nontraditional sets
lipid_class_traditional <- lipid_class_info %>%
  filter(Class %in% traditional_lipid_classes)

lipid_class_nontraditional <- lipid_class_info %>%
  filter(Class %in% nontraditional_lipid_class)

# Use final_name column now
traditional_lipids     <- lipid_class_traditional$final_name
nontraditional_lipids  <- lipid_class_nontraditional$final_name

# Assuming you renamed control/lowinput already
trad_cols_control     <- intersect(traditional_lipids, colnames(control))
nontrad_cols_control  <- intersect(nontraditional_lipids, colnames(control))

trad_cols_lowinput    <- intersect(traditional_lipids, colnames(lowinput))
nontrad_cols_lowinput <- intersect(nontraditional_lipids, colnames(lowinput))

# Subset
control_traditional     <- control  %>% select(Compound_Name, all_of(trad_cols_control))
control_nontraditional  <- control  %>% select(Compound_Name, all_of(nontrad_cols_control))
lowinput_traditional    <- lowinput %>% select(Compound_Name, all_of(trad_cols_lowinput))
lowinput_nontraditional <- lowinput %>% select(Compound_Name, all_of(nontrad_cols_lowinput))

colnames(control_traditional)

# — 1. Build a lookup table from your final_name → Class
lookup <- lipid_class_info %>%
  mutate(final_name = ifelse(is.na(CommonName), Lipids, CommonName)) %>%
  select(final_name, Class)

# — 2. Function to count lipids by Class for any subset df
count_classes <- function(df, cond_label, set_label) {
  tibble(final_name = setdiff(names(df), "Compound_Name")) %>%
    left_join(lookup, by = "final_name") %>%
    mutate(
      Class     = replace_na(Class, "Other"),
      Condition = cond_label,
      Set       = set_label
    ) %>%
    group_by(Set, Class, Condition) %>%
    summarise(Count = n(), .groups = "drop")
}

# — 3. Get counts for all four combinations
counts_trad     <- count_classes(control_traditional,    "Control",  "Traditional")
counts_nontrad  <- count_classes(control_nontraditional, "Control",  "Nontraditional")
counts_trad_L   <- count_classes(lowinput_traditional,  "LowInput", "Traditional")
counts_nontrad_L<- count_classes(lowinput_nontraditional,"LowInput", "Nontraditional")

class_counts <- bind_rows(counts_trad, counts_trad_L,
                          counts_nontrad, counts_nontrad_L) %>%
  group_by(Set, Class) %>%
  mutate(total = sum(Count)) %>%
  ungroup() %>%
  arrange(Set, desc(total)) %>%
  mutate(Class = fct_inorder(Class))

# — 4. Draw the figure, exactly in your Figure 1B style, but with one facet per Set
offset <- max(class_counts$Count) * 0.04

fig1b_updated <- ggplot(class_counts,
                        aes(x = fct_reorder(Class, total, .desc = TRUE),
                            y = Count,
                            fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.6,
           colour   = "black",
           linewidth = 0.25) +
  geom_text(aes(y = Count + offset, label = Count),
            position = position_dodge(width = 0.8),
            hjust     = 0,
            size      = 3.3) +
  coord_flip(clip = "off") +
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
  ) +
  facet_wrap(~ Set, scales = "free_y", ncol = 1)

# — 5. Display / save
quartz()
print(fig1b_updated)

ggsave("SuppFig3B_trad_nontrad_counts.png",
       plot = fig1b_updated,
       width = 6, height = 8, dpi = 300, bg = "white")





#### VENN DIAGRAM



library(ggvenn)
library(patchwork)

# 1) Build the two named lists for ggvenn()
venn_trad <- list(
  Control  = setdiff(names(control_traditional),    "Compound_Name"),
  LowInput = setdiff(names(lowinput_traditional),   "Compound_Name")
)

venn_nontrad <- list(
  Control  = setdiff(names(control_nontraditional), "Compound_Name"),
  LowInput = setdiff(names(lowinput_nontraditional),"Compound_Name")
)

# 2) Traditional lipids Venn
p_trad <- ggvenn(
  venn_trad,
  fill_color      = c("#440154FF", "#FDE725FF"),
  stroke_size     = 0.8,
  set_name_size   = 6,
  text_size       = 5,
  show_percentage = TRUE
) +
  ggtitle("Traditional Lipids") +
  theme(
    # push circles down by adding bottom space *under* the title
    plot.title      = element_text(
      face   = "bold",
      hjust  = 0.5,
      size   = 14,
      margin = margin(b = 20)    # <— 20 pts of space *below* the title
    ),
    plot.margin     = margin(t = 10, r = 20, b = 5, l = 20),
    legend.position = "none"
  )

p_nontrad <- ggvenn(
  venn_nontrad,
  fill_color      = c("#440154FF", "#FDE725FF"),
  stroke_size     = 0.8,
  set_name_size   = 6,
  text_size       = 5,
  show_percentage = TRUE
) +
  ggtitle("Non-Traditional Lipids") +
  theme(
    plot.title      = element_text(
      face   = "bold",
      hjust  = 0.5,
      size   = 14,
      margin = margin(b = 20)    # <— same trick here
    ),
    plot.margin     = margin(t = 10, r = 20, b = 5, l = 20),
    legend.position = "none"
  )
quartz(); print(p_trad)

# 6) (Optional) Save
# ggsave("Fig_venn_trad_nontrad.png", venn_figure, width = 10, height = 5, dpi = 300)

# Save the plot
ggsave("SuppFig_3C_Lipid_Overlap_Venn_traditional.png",
       plot = p_trad,
       width = 10, height = 10, dpi = 300,
       units = "in", bg = "white")  # white background for publication quality
# Save the plot
ggsave("SuppFig_3D_Lipid_Overlap_Venn_nontraditional.png",
       plot = p_nontrad,
       width = 10, height = 10, dpi = 300,
       units = "in", bg = "white")  # white background for publication quality



