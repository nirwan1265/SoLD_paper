###############################################################################
###############################################################################
## Figures
##   Fig1A – Lipid species count for Glycerolipid and Glycerophospholipid
##   Fig1B – Lipid class counts (Total)
##   Fig1C – Venn Diagram
###############################################################################
###############################################################################

###############################################################################
### Load the libraries
###############################################################################

# Load
library(vroom)
library(dplyr)
library(tidyr)
library(ggplot2)    
library(stringr)    
library(forcats)    
library(tidyverse)
library(ggpubr)     
library(scales) 
library(ggVennDiagram)  
library(readr) 
library(sf)
library(viridis)
library(ggvenn)
library(patchwork)

###############################################################################
### Read in the data
###############################################################################
# Control
control  <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
colnames(control)[1] <- "Compound_Name"  

# Lowinput
lowinput  <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
colnames(lowinput)[1] <- "Compound_Name"  


###############################################################################
### Load the Lipid Class
###############################################################################
lipid_class_info <- vroom::vroom("data/lipid_class/Final_lipid_classes.csv", show_col_types = FALSE) %>% dplyr::distinct(Lipids, .keep_all = TRUE)  

# Check the classes
table(lipid_class_info$SubClass)
table(lipid_class_info$Class)


###############################################################################
### Define aesthetics for plotting
###############################################################################
plot_theme <- theme_minimal(base_size = 16) +
  theme(
    plot.title     = element_text(
      size   = 16,
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

# Color scheme (colorblind-friendly)
group_colors <- viridis(2)
condition_colors <- c(Control = "#440154FF", LowInput = "#FDE725FF") 


###############################################################################
### Extract Glycerolipid and Glycerophosholipid species
###############################################################################

# Define the species
valid_classes <- c("TG","DG","MG",
                   "PC","PE","PG","PS",
                   "LPC","LPE",
                   "DGDG","MGDG","SQDG"
)
# Match Patter
class_pattern <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")
# Extract the lipids
extract_class <- function(x) str_extract(x, class_pattern)


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

# How far the text sits above the bars
offset <- max(class_counts$Count) * 0.04   # 3 % of the axis length


###############################################################################
### Figure 1A
###############################################################################

# Plot
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
  plot_theme +                           # your base theme
  theme(
    # grid lines
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "#d0d0d0", linewidth = 0.2),
    # legend tweaks
    legend.position     = c(0.78, 0.93),   
    legend.direction    = "horizontal",
    legend.background   = element_rect(fill = "white",
                                       colour = "black", linewidth = 0.2),
    legend.key.size     = unit(0.25, "cm"),
    legend.title        = element_blank(), 
    # axis text
    axis.text.y         = element_text(face = "bold")
  )

# Check
quartz()
print(fig1a)

# Save the plot
ggsave("fig/supp/SuppFig_4A_Lipid_Counts.png",
       plot = fig1a,
       width = 6, height = 6, dpi = 300,
       units = "in", bg = "white")  # white background for publication quality




###############################################################################
### Extract the lipid classes
###############################################################################

# Define which Classes
super_wanted <- c("Glycerolipid","Glycerophospholipid",
                  "Sphingolipid","Sterol","Betaine lipid",
                  "Fatty acid","Ether lipid","Terpenoid","Prenol")

# Build lookup of lipid classes
lipid_lookup <- lipid_class_info %>%
  filter(Class %in% super_wanted) %>%
  select(Lipids, Class)

# Function to count species per Class
get_class_counts <- function(df, label){
  df %>%
    select(-Compound_Name) %>%          # drop sample‐ID col
    names() %>%                          # get vector of lipid names
    enframe(name = NULL, value = "Lipids") %>%
    inner_join(lipid_lookup, by = "Lipids") %>%
    distinct(Lipids, Class) %>%
    count(Class, name = "Count") %>%
    mutate(Condition = label)
}

# Apply the function
counts_ctrl <- get_class_counts(control,  "Control")
counts_low  <- get_class_counts(lowinput, "LowInput")
class_counts <- bind_rows(counts_ctrl, counts_low)

# Fix order of Class on the axis
class_counts <- class_counts %>%
  mutate(Class = factor(Class, levels = super_wanted))



###############################################################################
### Figure 1B
###############################################################################

# Plot
offset <- max(class_counts$Count) * 0.04

fig1b <- ggplot(class_counts,
                aes(x = fct_reorder(Class, Count, .desc = TRUE),
                    y = Count,
                    fill = Condition)) +
  geom_col(position   = position_dodge(width = 0.8),
           width      = 0.6,
           colour     = "black",
           linewidth  = 0.25) +
  geom_text(aes(label = Count,
                y     = Count + offset),
            position = position_dodge(width = 0.8),
            hjust    = 0,
            size     = 3.3) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = condition_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL,
       y = "Number of lipid species") +
  theme_classic(base_size = 14) +
  plot_theme +
  theme(
    panel.grid.major.x     = element_line(color = "#d0d0d0", size = 0.2),
    panel.grid.major.y     = element_blank(),
    panel.grid.minor       = element_blank(),
    axis.text.y            = element_text(face = "bold"),
    legend.position        = c(0.78, 0.93),
    legend.direction       = "horizontal",
    legend.key.size        = unit(0.25, "cm"),
    legend.title           = element_blank(),
    legend.background      = element_rect(fill = "white",
                                          colour = "black",
                                          linewidth = 0.2),
    legend.box.background  = element_rect(fill = "white",
                                          colour = "black",
                                          linewidth = 0.2)
  )

# Check
quartz()
print(fig1b)

# Save the plot
ggsave("fig/supp/SuppFig_4B_Lipid_Class_Counts.png",
       plot = fig1b,
       width = 8, height = 5, dpi = 300,
       units = "in", bg = "white")  # white background for publication quality



###############################################################################
### Venn diagram
###############################################################################

# Build the two named lists for lipids 
venn_lipid <- list(
  Control = setdiff(names(control), "Compound_Name"),
  LowInput = setdiff(names(lowinput), "Compound_Name")
)


###############################################################################
### Figure 1C
###############################################################################

# Plot
p_lipid <- ggvenn(
  venn_lipid,
  fill_color      = c("#440154FF", "#FDE725FF"),
  stroke_size     = 0.8,
  set_name_size   = 6,
  text_size       = 5,
  show_percentage = TRUE
) +
  #ggtitle("All Detected Lipids") +
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
  ) +
  plot_theme +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())


# Check
quartz()
p_lipid

# Save the plot
ggsave("fig/supp/SuppFig_4C_Lipid_Overlap_Venn.png",
       plot = p_lipid,
       width = 6, height = 6, dpi = 300,
       units = "in", bg = "white")  
