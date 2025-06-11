#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
################################################################################
#### SOrghum Lipidomics Database (SOLD)
################################################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

################################################################################
################################################################################
#### FIGURE 1A
################################################################################
################################################################################

###############################################################################
# 1.  LIBRARIES
###############################################################################

suppressPackageStartupMessages({
  library(vroom);   library(dplyr);   library(tidyr);  library(stringr)
  library(ggplot2); library(viridis); library(gridExtra); library(grid)
})

################################################################################
# 1.  LOAD DATA  (Control & LowInput)
################################################################################
control  <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Control_all_lipids_final_non_normalized.csv")
lowinput <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv")


################################################################################
# 2.  SELECT THE CLASSES 
################################################################################

valid_classes <- c("TG","DG","MG","PC","PE","PG","PI",
                   "LPC","LPE","DGDG","MGDG","Cer","SM","FA","SQDG","AEG","PA","PS")
                   
valid_classes <- c("PC","PE","PA","PG","PS")
#valid_classes <- c("TG","DG","MG","PC","PE",
#                   "DGDG","MGDG","SQDG")

class_pat <- paste0("\\b(", paste(valid_classes, collapse = "|"), ")\\b")

reshape_plate <- function(df, label){
  df %>%
    pivot_longer(-Compound_Name, names_to = "Lipid", values_to = "Intensity") %>%
    mutate(Class = str_extract(Lipid, class_pat)) %>%
    filter(!is.na(Class)) %>%
    group_by(Sample = Compound_Name, Class) %>%
    summarise(sum_int = sum(Intensity, na.rm = TRUE), .groups = "drop") %>%
    mutate(Condition = label)
}

ctrl_long <- reshape_plate(control,  "Control")
low_long  <- reshape_plate(lowinput, "LowInput")


###############################################################################
# 3.  % OF TIC  → mean %  per condition
###############################################################################

pct_tbl <- bind_rows(ctrl_long, low_long) %>%
  group_by(Sample, Condition) %>%
  mutate(pct = 100 * sum_int / sum(sum_int)) %>%
  ungroup()

pct_mean <- pct_tbl %>%
  group_by(Condition, Class) %>%
  summarise(pct_mean = mean(pct), .groups = "drop")

###############################################################################
# 4.  FIXED PALETTE AND THEME
###############################################################################

classes_all <- sort(unique(pct_mean$Class))
pal          <- viridis(length(classes_all), option = "D")
names(pal)   <- classes_all

nature_theme <- theme_minimal(base_size = 14) +
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

###############################################################################
# 5.  STACKED BAR PLOT (labels ≥ 3 %)
###############################################################################
threshold <- 3

p_bar <- ggplot(pct_mean,
                aes(x = pct_mean, y = Condition,
                    fill = factor(Class, levels = classes_all))) +
  geom_col(width = .8, colour = "white") +
  geom_text(aes(label = ifelse(pct_mean >= threshold,
                               sprintf("%.1f%%", pct_mean), "")),
            position = position_stack(vjust = 0.5),
            colour   = "white", size = 3) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,100),
                     breaks = seq(0,100,25),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "% of TIC", y = NULL) +
  theme_classic(base_size = 14) +
  theme(axis.line.x = element_line(size = .4),
        axis.line.y = element_line(size = .4),
        axis.text.x    = element_text(face = "bold", size = 14),  # enlarge x‐axis tick labels
        axis.text.y    = element_text(face = "bold", size = 14),  # enlarge y‐axis tick labels
        axis.title.x   = element_text(face = "bold", size = 16)   # enlarge x‐axis label
  ) 
###############################################################################
# 6.  NUMERIC TABLE  +  COLOUR SWATCH COLUMN FOR LEGENDS
###############################################################################

tbl <- pct_mean %>%
  # Keep the raw numeric pct_mean so we can test < 0.05:
  mutate(pct_numeric = pct_mean) %>%
  
  # Convert numeric → string, with "<0.1 %" whenever pct_numeric < 0.05
  mutate(pct = case_when(
    pct_numeric < 0.05 ~ "<0.1 %",
    TRUE               ~ sprintf("%.1f %%", pct_numeric)
  )) %>%
  
  # Now drop the helper column and rename Class → Species
  select(Species = Class, Condition, pct) %>%
  
  # Pivot so that each row is one Species, with separate Control/LowInput columns
  pivot_wider(
    names_from  = Condition,
    values_from = pct,
    values_fill = ""    # (if you prefer “” instead of “NA”)
  ) %>%
  
  # Arrange rows in the same order as classes_all (your existing palette order)
  arrange(match(Species, classes_all))

# ── 2. Prepend an empty “Swatch” column so we can draw coloured squares later
tbl <- tibble(Swatch = "", !!!tbl)

# ── 3. Create the tableGrob with larger text sizes
tg <- tableGrob(
  tbl,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 14,
    core = list(fg_params = list(fontsize = 14, hjust = 0.5)),
    colhead = list(fg_params = list(fontsize = 14, fontface = "bold"))
  )
)

# ── 4. Overwrite each “Swatch” cell in column 1 with a coloured rectangle
core_rows <- which(tg$layout$name == "core-bg" & tg$layout$l == 1)
for (i in seq_along(core_rows)) {
  tg$grobs[[ core_rows[i] ]] <- rectGrob(
    width = unit(4, "mm"), height = unit(4, "mm"),
    gp = gpar(
      fill = pal[ tbl$Species[i] ],  # use tbl$Species here
      col  = pal[ tbl$Species[i] ]
    )
  )
}


###############################################################################
# 7.  LAYOUT:  STACKED BAR | TABLE
###############################################################################

fig1a <- arrangeGrob(                       
  p_bar, tg,
  ncol   = 2,
  widths = unit.c(unit(1, "null"), unit(0.30, "null"))
)

# view in R
quartz()   
grid::grid.draw(fig1a)


###############################################################################
# 8.  SAVE THE PLOT
###############################################################################

#ggsave("Fig1a_lipid_species.png",fig1a, width = 21, height = 8, units = "in", bg = "white")
ggsave("SuppFig_TIC_PC_PA_PE_PG_PS.png",fig1a, width = 21, height = 8, units = "in", bg = "white")


################################################################################
################################################################################
#### FIGURE 1B
################################################################################
################################################################################

###############################################################################
# 1.  libraries
###############################################################################

suppressPackageStartupMessages({
  library(vroom);    library(dplyr);    library(tidyr);   library(stringr)
  library(ggplot2);  library(viridis);  library(gridExtra); library(grid)
})

###############################################################################
# 2.  LOAD THE DATA
###############################################################################

control   <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Control_all_lipids_final_non_normalized.csv")
lowinput  <- vroom("/Users/nirwantandukar/Documents/Research/data/SAP/non_normalized_intensities/Lowinput_all_lipids_final_non_normalized.csv")


lipid_class_info <- vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/data/lipid_class.csv",
                          show_col_types = FALSE) %>% 
  filter(!is.na(SubClass)) %>% 
  transmute(Lipid = Lipids, SuperClass = Class)


# keep only the five super-classes you want on the pie
wanted_super <- c("Terpenoid",
                  "Glycerophospholipid",
                  "Glycerolipid",
                  "Glycoglycerolipid",
                  "Fatty acid and derivative","Sphingolipid",
                  "Steroid", "Tetrapyrrole", "Sterol",
                  "N-acylethanolamine", "Betaine lipid","Vitamin" )

lipid_class_info <- lipid_class_info %>% 
  filter(SuperClass %in% wanted_super)

###############################################################################
# 3.  HELPER FUNCTION TO SUM THE PER CLASS
###############################################################################

reshape_plate_SC <- function(df, label){
  df %>%
    pivot_longer(-Compound_Name, names_to = "Lipid", values_to = "Intensity") %>%
    left_join(lipid_class_info, by = "Lipid") %>%
    filter(!is.na(SuperClass)) %>%
    group_by(Sample = Compound_Name, SuperClass) %>%
    summarise(sum_int = sum(Intensity, na.rm = TRUE), .groups = "drop") %>%
    mutate(Condition = label)
}

ctrl_SC  <- reshape_plate_SC(control,  "Control")
low_SC   <- reshape_plate_SC(lowinput, "LowInput")


###############################################################################
# 4.  PERC OF TOTAL INTENSITIES PER SAMPLE → mean % per condition
###############################################################################

pct_SC <- bind_rows(ctrl_SC, low_SC) %>%
  group_by(Sample, Condition) %>%
  mutate(pct = 100 * sum_int / sum(sum_int)) %>%
  ungroup()

pct_SC_mean <- pct_SC %>% 
  group_by(Condition, SuperClass) %>% 
  summarise(pct_mean = mean(pct), .groups = "drop")

###############################################################################
# 5.  FIXED PALETTE
###############################################################################

classes_all <- wanted_super
pal          <- viridis(length(classes_all), option = "D")
names(pal)   <- classes_all

###############################################################################
#  6.  STACKED BAR PLOT (labels ≥ 3 %)
###############################################################################
threshold <- 3         # print label only if slice ≥ 3 %

p_sc <- ggplot(pct_SC_mean,
               aes(x = pct_mean,
                   y  = Condition,
                   fill = factor(SuperClass, levels = classes_all))) +
  geom_col(width = .8, colour = "white") +
  geom_text(aes(label = ifelse(pct_mean >= threshold,
                               sprintf("%.1f%%", pct_mean), "")),
            position = position_stack(vjust = 0.5),
            colour = "white", size = 3.2) +
  scale_fill_manual(values = pal, guide = "none") +   # ← remove legend
  scale_x_continuous(expand = c(0,0),
                     limits  = c(0,100),
                     breaks  = seq(0,100,25),
                     labels  = function(x) paste0(x, "%")) +
  labs(x = "% of TIC", y = NULL) +
  theme_classic(base_size = 14) +
  theme(axis.line.x = element_line(size = .4),
        axis.line.y = element_line(size = .4),
        axis.text.x    = element_text(face = "bold", size = 14),  # enlarge x‐axis tick labels
        axis.text.y    = element_text(face = "bold", size = 14),  # enlarge y‐axis tick labels
        axis.title.x   = element_text(face = "bold", size = 16))   # enlarge x‐axis label)

###############################################################################
# 7.  NUMERIC TABLE  +  COLOUR SWATCH COLUMN FOR LEGENDS
###############################################################################

tbl_sc <- pct_SC_mean %>%
  group_by(SuperClass, Condition) %>%
  summarise(pct = mean(pct_mean), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = Condition,
    values_from = pct,
    values_fill = 0
  ) %>%
  # Format each numeric percentage, replacing anything < 0.05 (~0.0%) with "<0.1 %"
  mutate(across(
    Control:LowInput,
    ~ case_when(
      . < 0.05 ~ "<0.1 %",
      TRUE     ~ sprintf("%.1f %%", .)
    )
  )) %>%
  # Arrange in your desired order, then rename SuperClass → Class
  arrange(match(SuperClass, classes_all)) %>%
  rename(Class = SuperClass) %>%
  # Add blank Swatch column and reorder
  mutate(Swatch = "") %>%
  select(Swatch, Class, Control, LowInput)

tg_sc <- tableGrob(
  tbl_sc, rows = NULL,
  theme = ttheme_minimal(
    base_size = 14,
    core    = list(fg_params = list(fontsize = 14, hjust = 0.5)),
    colhead = list(fg_params = list(fontsize = 14, fontface = "bold"))
  )
)

core_rows <- which(tg_sc$layout$name == "core-bg" & tg_sc$layout$l == 1)
for (i in seq_along(core_rows)) {
  tg_sc$grobs[[core_rows[i]]] <-
    rectGrob(
      width = unit(4, "mm"),
      height = unit(4, "mm"),
      gp = gpar(
        fill = pal[tbl_sc$Class[i]],    # ← now use the renamed “Class” column
        col  = pal[tbl_sc$Class[i]]
      )
    )
}

###############################################################################
# 8.  LAYOUT:  STACKED BAR | TABLE
###############################################################################

fig1b <- gridExtra::arrangeGrob(
  p_sc, tg_sc,
  ncol   = 2,
  widths = unit.c(unit(1, "null"),  # bar stretches
                  unit(0.30, "null")))  # table column

# view in R
quartz()   
grid::grid.draw(fig1b)


###############################################################################
# 9.  SAVE THE PLOT
###############################################################################

ggplot2::ggsave("Fig1b_lipid_class.png",
                fig1b,
                width = 21, height = 8, units = "in", bg = "white")








library(vroom); library(dplyr); library(tidyr)
library(ggplot2); library(viridis); library(patchwork); library(stringr)

# ── 1.  Load the class dictionary with Super- & SubClass info ────────────────
lipid_info <- vroom("/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Github/SoLD/data/lipid_class.csv",
                    show_col_types = FALSE) %>%
  filter(!is.na(SubClass)) %>% 
  filter(Class %in% c("Terpenoid",
                      "Glycerophospholipid",
                      "Glycerolipid",
                      "Glycoglycerolipid","Tetrapyrrole")) %>%
  
  transmute(Lipid = Lipids,
            SuperClass = Class,
            SubClass)

# ── 2.  Helper: sample × SubClass sums for one condition ─────────────────────
reshape_plate_sub <- function(df, label){
  df %>%
    pivot_longer(-Compound_Name,
                 names_to  = "Lipid",
                 values_to = "Intensity") %>%
    left_join(lipid_info, by = "Lipid") %>%        # add Super/Sub info
    filter(!is.na(SuperClass)) %>%                 # keep mapped lipids
    group_by(Sample = Compound_Name,
             Condition = label,
             SuperClass,
             SubClass) %>%                         # <- aggregation level
    summarise(sum_int = sum(Intensity, na.rm = TRUE),
              .groups = "drop")
}

control_sub  <- reshape_plate_sub(control,  "Control")
lowinput_sub <- reshape_plate_sub(lowinput, "LowInput")

# ── 3.  Convert to % of SuperClass total per sample, then average ─────────────
pct_tbl <- bind_rows(control_sub, lowinput_sub) %>%
  group_by(Sample, Condition, SuperClass) %>%
  mutate(pct = 100 * sum_int / sum(sum_int)) %>%
  ungroup()

pct_mean <- pct_tbl %>%
  group_by(Condition, SuperClass, SubClass) %>%
  summarise(pct_mean = mean(pct), .groups = "drop")


###############################################################################
# 0. libraries (already loaded above)  – no changes
###############################################################################

###############################################################################
# 1–3  load + reshape  (exactly the code you already have)
###############################################################################
# ── 1. class dictionary kept as `lipid_info`
# ── 2. reshape_plate_sub() + control_sub / lowinput_sub
# ── 3. pct_tbl  →  pct_mean  (mean % of superclass total)
#     pct_mean columns:  Condition | SuperClass | SubClass | pct_mean
###############################################################################

###############################################################################
# 4.  single, fixed palette for every sub-class
###############################################################################
subclasses_all <- sort(unique(pct_mean$SubClass))
pal_sub        <- viridis::viridis(length(subclasses_all), option = "C")
names(pal_sub) <- subclasses_all

###############################################################################
# 1. helper: build one “bar + table” grob   for a given super-class
###############################################################################
make_panel <- function(df_long, super_name,
                       thr = 3,            # label threshold
                       pal_fun = viridis::viridis) {
  
  ## -------- data for this super-class ----------
  dat <- df_long %>% filter(SuperClass == super_name)
  
  ## 1-A  local palette for its sub-classes
  sub_levels <- sort(unique(dat$SubClass))
  pal_local  <- pal_fun(length(sub_levels), option = "D")
  names(pal_local) <- sub_levels
  
  ## 1-B  bar plot (no legend)
  p <- ggplot(dat,
              aes(x = pct_mean, y = Condition,
                  fill = factor(SubClass, levels = sub_levels))) +
    geom_col(width = .8, colour = "white") +
    geom_text(aes(label = ifelse(pct_mean >= thr,
                                 sprintf("%.1f%%", pct_mean), "")),
              position = position_stack(vjust = 0.5),
              colour = "white", size = 2.8) +
    scale_x_continuous(expand = c(0,0), limits = c(0,100),
                       breaks = c(0,25,50,75,100),
                       labels = function(x) paste0(x,"%")) +
    scale_fill_manual(values = pal_local, guide = "none") +
    labs(title = super_name, x = NULL, y = NULL) +
    theme_classic(base_size = 14) +
    theme(plot.title  = element_text(face = "bold", hjust = .5),
          axis.line.x = element_line(size=.35),
          axis.line.y = element_line(size=.35))
  
  ## 1-C  numeric table
  tbl <- dat %>%
    select(SubClass, Condition, pct_mean) %>%
    mutate(pct = sprintf("%.1f %%", pct_mean)) %>%
    select(-pct_mean) %>%
    pivot_wider(names_from = Condition, values_from = pct) %>%
    arrange(match(SubClass, sub_levels)) %>%
    mutate(Swatch = "") %>%
    select(Swatch, SubClass, Control, LowInput)
  
  tg <- tableGrob(tbl, rows = NULL,
                  theme = ttheme_minimal(base_size = 14,
                                         core    = list(fg_params=list(hjust=.5)),
                                         colhead = list(fg_params=list(fontface="bold"))))
  
  # place a coloured square in every core row, column 1
  core_rows <- which(tg$layout$name == "core-bg" & tg$layout$l == 1)
  for (i in seq_along(core_rows))
    tg$grobs[[core_rows[i]]] <-
    rectGrob(width = unit(3.5,"mm"), height = unit(3.5,"mm"),
             gp=gpar(fill = pal_local[tbl$SubClass[i]],
                     col  = pal_local[tbl$SubClass[i]]))
  
  # header cell blank
  tg$grobs[[ which(tg$layout$name=="colhead-fg" & tg$layout$l==1) ]] <-
    textGrob(" ")
  
  ## 1-D  combine bar + table
  arrangeGrob(p, tg, ncol=2,
              widths = unit.c(unit(1,"null"), unit(0.55,"null")))
}

###############################################################################
# 2. build panels for every super-class
###############################################################################
supers <- c("Terpenoid","Glycerophospholipid",
            "Glycerolipid","Glycoglycerolipid","Tetrapyrrole")

panel_list <- lapply(supers, make_panel, df_long = pct_mean)

###############################################################################
# 3. arrange everything   (2 columns × 5 rows)
###############################################################################
library(gridExtra)

combo <- arrangeGrob(grobs = panel_list,
                     ncol = 2,
                     top = textGrob("Sub-class composition of each super-class",
                                    gp=gpar(fontface="bold", fontsize=10)))

quartz()                      # view
grid::grid.draw(combo)

ggsave("subclass_panels_with_tables.png",
       combo, width = 16, height = 6, units = "in", bg = "white")


# find version of R
R.version.string



