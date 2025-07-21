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

control <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  select(-c(2,3,4)) %>%
  rename(Compound_Name = 1)

lowinput <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>%
  select(-c(2,3,4)) %>%
  rename(Compound_Name = 1) 


################################################################################
# 2.  SELECT THE CLASSES 
################################################################################

valid_classes <- c("TG","DG","MG",
                   "PC","PE","PI",
                   "DGDG","MGDG",
                   "SQDG","SM","AEG",
                   "LPC","LPE","PG","PA","Cer","GalCer","FA")

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

# Remove rows with string less than 5 characters from Sample
pct_tbl <- pct_tbl %>%
  filter(nchar(Sample) >= 5)


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
  theme_classic(base_size = 16) +
  theme(axis.line.x = element_line(size = .4),
        axis.line.y = element_line(size = .4),
        axis.text.x    = element_text(face = "bold", size = 16),  # enlarge x‐axis tick labels
        axis.text.y    = element_text(face = "bold", size = 16),  # enlarge y‐axis tick labels
        axis.title.x   = element_text(face = "bold", size = 16)   # enlarge x‐axis label
  ) + 
  nature_theme

###############################################################################
# 6.  NUMERIC TABLE  +  COLOUR SWATCH COLUMN FOR LEGENDS
###############################################################################

tbl <- pct_mean %>%
  # Keep the raw numeric pct_mean so we can test < 0.05:
  dplyr::mutate(pct_numeric = pct_mean) %>%
  
  # Convert numeric → string, with "<0.1 %" whenever pct_numeric < 0.05
  dplyr::mutate(pct = case_when(
    pct_numeric < 0.05 ~ "<0.1 %",
    TRUE               ~ sprintf("%.1f %%", pct_numeric)
  )) %>%
  
  # Now drop the helper column and rename Class → Species
  dplyr::select(Species = Class, Condition, pct) %>%
  
  # Pivot so that each row is one Species, with separate Control/LowInput columns
  pivot_wider(
    names_from  = Condition,
    values_from = pct,
    values_fill = ""    # (if you prefer “” instead of “NA”)
  ) %>%
  
  # Arrange rows in the same order as classes_all (your existing palette order)
  dplyr::arrange(match(Species, classes_all))

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

ggsave("fig/main/Fig1a_lipid_species.png",fig1a, width = 21, height = 8, units = "in", bg = "white")
#ggsave("SuppFig_TIC_PC_PA_PE_PG_PS.png",fig1a, width = 21, height = 8, units = "in", bg = "white")








################################################################################
################################################################################
#### Supplementary Figure 4
################################################################################
################################################################################


### Lipid Class
lipid_class_info <- vroom("data/lipid_class/final_lipid_classes.csv",
                          show_col_types = FALSE) %>%
  #filter(!is.na(SubClass)) %>%
  transmute(Lipid = Lipids, SuperClass = Class)

# All the super class
wanted_super <- c("Glycerolipid", "Glycerophospholipid",
                   "Sphingolipid", "Sterol", "Betaine lipid", "Fatty acid","Ether lipid","Terpenoid")

lipid_class_info <- lipid_class_info %>% 
  filter(SuperClass %in% wanted_super)


### Sum per class
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


### Percentage of total intensities per sample → mean % per condition
pct_SC <- bind_rows(ctrl_SC, low_SC) %>%
  group_by(Sample, Condition) %>%
  mutate(pct = 100 * sum_int / sum(sum_int)) %>%
  ungroup()

pct_SC_mean <- pct_SC %>% 
  group_by(Condition, SuperClass) %>% 
  summarise(pct_mean = mean(pct), .groups = "drop")


### Get the classes and define the colors
classes_all <- wanted_super
pal          <- viridis(length(classes_all), option = "D")
names(pal)   <- classes_all


### Plot (labels ≥ 3 %)
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
        axis.title.x   = element_text(face = "bold", size = 16)) +  # enlarge x‐axis label)
  nature_theme


### Numeric table + color swatch for legends
tbl_sc <- pct_SC_mean %>%
  dplyr::group_by(SuperClass, Condition) %>%
  dplyr::summarise(pct = mean(pct_mean), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = Condition,
    values_from = pct,
    values_fill = 0
  ) %>%
  # Format each numeric percentage, replacing anything < 0.05 (~0.0%) with "<0.1 %"
  dplyr::mutate(across(
    Control:LowInput,
    ~ case_when(
      . < 0.05 ~ "<0.1 %",
      TRUE     ~ sprintf("%.1f %%", .)
    )
  )) %>%
  # Arrange in your desired order, then rename SuperClass → Class
  dplyr::arrange(match(SuperClass, classes_all)) %>%
  dplyr::rename(Class = SuperClass) %>%
  # Add blank Swatch column and reorder
  dplyr::mutate(Swatch = "") %>%
  dplyr::select(Swatch, Class, Control, LowInput)

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


### Layout of the plot
fig1b <- gridExtra::arrangeGrob(
  p_sc, tg_sc,
  ncol   = 2,
  widths = unit.c(unit(1, "null"),  # bar stretches
                  unit(0.30, "null")))  # table column

# view in R
quartz()   
grid::grid.draw(fig1b)


### Save the plot
ggplot2::ggsave("fig/supp/SuppFig4A_lipid_class.png",
                fig1b,
                width = 21, height = 8, units = "in", bg = "white")






###############################################################################
### Supplementary Figure 5
###############################################################################

### Load the class dictionary with Super- & SubClass info
lipid_info <- vroom("data/lipid_class/final_lipid_classes.csv",
                    show_col_types = FALSE) %>%
  filter(!is.na(SubClass)) %>% 
  transmute(Lipid = Lipids,
            SuperClass = Class,
            SubClass)

# SUM intensities per Sample × Class × SubClass ─────────────────────────
agg_plate <- function(df, condition_label) {
  df %>%
    pivot_longer(
      -Compound_Name,
      names_to  = "Lipid",
      values_to = "Intensity"
    ) %>%
    left_join(lipid_info, by = "Lipid") %>%
    filter(!is.na(SuperClass)) %>%
    group_by(
      Sample    = Compound_Name,
      Condition = condition_label,
      SuperClass,
      SubClass
    ) %>%
    summarise(sum_int = sum(Intensity, na.rm = TRUE), .groups = "drop")
}

control_sub  <- agg_plate(control,  "Control")
lowinput_sub <- agg_plate(lowinput, "LowInput")

# ── 4. COMPUTE percent of each Class’s total per sample ─────────────────────────
pct_tbl <- bind_rows(control_sub, lowinput_sub) %>%
  group_by(Sample, Condition, SuperClass) %>%
  mutate(
    total_cls = sum(sum_int, na.rm = TRUE),
    pct       = if_else(total_cls > 0, 100 * sum_int / total_cls, 0)
  ) %>%
  ungroup()

# ── 5. AVERAGE across samples to get mean % per Condition × Class × SubClass ──
pct_mean <- pct_tbl %>%
  group_by(Condition, SuperClass, SubClass) %>%
  summarise(pct_mean = mean(pct, na.rm = TRUE), .groups = "drop") %>%
  # replace any NaN (0/0) with 0
  mutate(pct_mean = if_else(is.nan(pct_mean), 0, pct_mean))

# ── 6. SET UP color palettes ───────────────────────────────────────────────────
classes_all <- unique(pct_mean$SuperClass)
pal_class   <- viridis(length(classes_all), option = "D")
names(pal_class) <- classes_all

# ── 7. FUNCTION to build one Class‑panel ──────────────────────────────────────
make_class_panel <- function(df, cls, label_thresh = 3) {
  # 1) subset & ensure both Conditions
  dat <- df %>%
    filter(SuperClass == cls) %>%
    complete(
      SubClass,
      Condition = c("Control","LowInput"),
      fill      = list(pct_mean = NA_real_)
    )
  
  # 2) build a local palette
  subs    <- sort(unique(dat$SubClass))
  pal_sub <- viridis(length(subs), option = "C")
  names(pal_sub) <- subs
  
  # 3) bar plot with <0.1% labels
  p <- ggplot(dat, aes(x = pct_mean, y = Condition, fill = SubClass)) +
    geom_col(colour = "white", width = .8, na.rm = TRUE) +
    geom_text(aes(
      label = case_when(
        pct_mean >= label_thresh                ~ sprintf("%.1f%%", pct_mean),
        pct_mean > 0 & pct_mean < label_thresh  ~ "<0.1%",
        TRUE                                    ~ ""
      )
    ),
    position = position_stack(vjust = .5),
    colour = "white", size = 3, na.rm = TRUE) +
    scale_fill_manual(values = pal_sub) +
    scale_x_continuous(
      expand = c(0,0),
      limits = c(0,100),
      breaks = seq(0,100,25),
      labels = function(x) paste0(x,"%")
    ) +
    labs(title = cls, x = NULL, y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      plot.title   = element_text(face="bold", hjust=.5),
      axis.text    = element_text(size=12)
    )
  
  # 4) numeric table with <0.1% formatting
  tbl <- dat %>%
    select(SubClass, Condition, pct_mean) %>%
    mutate(
      pct = case_when(
        is.na(pct_mean)          ~ NA_character_,
        pct_mean < 0.1 & pct_mean > 0 ~ "<0.1 %",
        TRUE                     ~ sprintf("%.1f %%", pct_mean)
      )
    ) %>%
    select(-pct_mean) %>%
    pivot_wider(
      names_from  = Condition,
      values_from = pct,
      values_fill = list(pct = NA_character_)
    ) %>%
    arrange(match(SubClass, subs)) %>%
    mutate(Swatch = "") %>%
    select(Swatch, SubClass, Control, LowInput)
  
  tg <- tableGrob(
    tbl, rows = NULL,
    theme = ttheme_minimal(
      core    = list(fg_params=list(hjust=.5, fontsize=12)),
      colhead = list(fg_params=list(fontface="bold", fontsize=13))
    )
  )
  
  # fill the swatch squares
  core_rows <- which(tg$layout$name == "core-bg" & tg$layout$l == 1)
  for (i in seq_along(core_rows)) {
    fill_col <- pal_sub[tbl$SubClass[i]]
    tg$grobs[[core_rows[i]]] <-
      rectGrob(
        gp    = gpar(fill = fill_col, col = fill_col),
        width = unit(4, "mm"),
        height= unit(4, "mm")
      )
  }
  
  # blank the Swatch header
  head_idx <- which(tg$layout$name=="colhead-fg" & tg$layout$l==1)
  tg$grobs[[head_idx]] <- textGrob(" ")
  
  # 5) combine bar + table
  arrangeGrob(
    p, tg,
    ncol   = 2,
    widths = unit.c(unit(1, "null"), unit(0.4, "null"))
  )
}


# ── 8. BUILD ALL PANELS ───────────────────────────────────────────────────────
panel_list <- lapply(classes_all, make_class_panel, df = pct_mean)

# ── 9. DISPLAY (e.g. 2 columns) ───────────────────────────────────────────────
quartz()
grid.arrange(grobs = panel_list, ncol = 2)

quartz()

#save
ggsave("fig/supp/SuppFig5_lipid_class_subclass.png",
       arrangeGrob(grobs = panel_list, ncol = 2),
       width = 40, height = 16, units = "in", bg = "white")


