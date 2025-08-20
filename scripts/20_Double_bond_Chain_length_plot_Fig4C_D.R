library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
# --- If needed: install.packages(c("dplyr","tidyr","stringr","ggplot2"))


####### DOUBLE BOND



# Read the raw files 
control  <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_control_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
colnames(control)[1] <- "Compound_Name"  


lowinput  <- vroom("data/SPATS_fitted/non_normalized_intensities/Final_subset_lowinput_all_lipids_fitted_phenotype_non_normalized.csv") %>% dplyr::select(-c(2,3,4))
colnames(lowinput)[1] <- "Compound_Name"  


# Theme
plot_theme <- theme_minimal(base_size = 24) +
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
    
    legend.position      = c(0.95, 0.95),
    legend.justification = c("right","top"),
    legend.background    = element_rect(fill="white", color="grey70", size=0.4),
    legend.direction     = "vertical",
    legend.spacing.y     = unit(0.2,"cm"),
    legend.title         = element_blank(),
    legend.text          = element_text(size=16),
    
    plot.margin    = margin(15, 15, 15, 15)
  )



# 0) Bind Control + LowInput (you already have these tibbles)
combined <- bind_rows(control  %>% mutate(Condition = "Control"),
                      lowinput %>% mutate(Condition = "LowInput"))

# (optional) limit to lipid classes you care about
valid_classes <- c("TG","DG","MG","PC","PE","PI","PG","PA","PS","MGDG","DGDG","SQDG","LPC","LPE")
class_pat <- paste0("\\b(", paste(valid_classes, collapse="|"), ")\\b")

# 1) Long format at the **species** level (do NOT aggregate by class here)
long_species <- combined %>%
  pivot_longer(-c(Compound_Name, Condition),
               names_to = "Lipid", values_to = "Intensity") %>%
  dplyr::filter(str_detect(Lipid, class_pat)) %>%     # keep only lipid species
  dplyr::rename(Sample = Compound_Name)

# 2) TIC normalize within each sample, then log10
long_species <- long_species %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(TIC = sum(Intensity, na.rm = TRUE),
         rel_abund = Intensity / TIC) %>%
  ungroup() %>%
  group_by(Sample) %>%
  dplyr::mutate(eps = {mp <- min(rel_abund[rel_abund > 0], na.rm = TRUE); ifelse(is.finite(mp), mp*0.5, 0)},
         log_rel = log10(rel_abund + eps)) %>%
  ungroup()

# 3) Z-score **per lipid species** across all samples (Control + LowInput)
long_species <- long_species %>%
  dplyr::group_by(Lipid) %>%
  dplyr::mutate(z = (log_rel - mean(log_rel, na.rm = TRUE)) / sd(log_rel, na.rm = TRUE)) %>%
  ungroup()

# 4) Compute delta-Z per species: LowInput mean(Z) - Control mean(Z)
effect_by_species <- long_species %>%
  dplyr::group_by(Lipid, Condition) %>%
  dplyr::summarise(mean_z = mean(z, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = mean_z) %>%
  dplyr::mutate(delta_z = LowInput - Control)

# 5) Parse total number of double bonds from lipid names
#    Works for MG(18:3), DG(16:0/18:2), TG(16:0/18:1/18:2), MGDG(18:3/18:3), etc.
count_unsat <- function(name) {
  m <- str_match_all(name, "(\\d+):(\\d+)")[[1]]
  if (nrow(m) == 0) return(NA_real_)
  sum(as.numeric(m[,3]))
}
effect_by_species <- effect_by_species %>%
  mutate(total_unsat = vapply(Lipid, count_unsat, numeric(1))) %>%
  filter(!is.na(total_unsat) & is.finite(delta_z))

# (optional) also keep the class for faceting later if you want
effect_by_species <- effect_by_species %>%
  mutate(Class = str_extract(Lipid, class_pat))

# Remove LPC, LPE, P, PG, PS
effect_by_species <- effect_by_species %>%
  filter(!Class %in% c("LPC", "LPE", "PG", "PA", "PS"))

# 6A) Global curve (all classes pooled), LOESS smooth of delta-Z vs total unsaturation
quartz()
ggplot(effect_by_species, aes(x = total_unsat, y = delta_z)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "loess", span = 0.8, se = TRUE, color = "blue") +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Total unsaturated bonds",
       y = "ΔZ (LowInput – Control)",
       title = "Condition effect vs total unsaturation (species-level)") +
  plot_theme

# Save the plot
ggsave("fig/main/Unsaturation_condition.png", width = 10, height = 6)


# 6B) If you’d rather stratify by class (nice biological readout):
quartz()
ggplot(effect_by_species, aes(total_unsat, delta_z)) +
  geom_point(alpha = 0.15) +
  geom_smooth(method = "loess", span = 0.8, se = TRUE, color = "blue") +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~ Class, scales = "free_y") +
  labs(x = "Total unsaturated bonds", y = "ΔZ (LowInput – Control)") +
  plot_theme

# Save the faceted plot
ggsave("fig/main/Unsaturation_individual.png", width = 10, height = 6)









####### CHAIN LENGTH









library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# classes to keep (edit if you want to restrict)
valid_classes <- c("TG","DG","MG","PC","PE","PI","DGDG","MGDG","SQDG","LPC","LPE","PG","PA","PS")
class_pat <- paste0("\\b(", paste(valid_classes, collapse="|"), ")\\b")

# 1) total chain length parser: sum the "carbons" from every N:N in the name
sum_carbons <- function(x) {
  m <- stringr::str_match_all(x, "(\\d+):(\\d+)")[[1]]
  if (nrow(m) == 0) return(NA_real_)             # no acyl info
  sum(as.numeric(m[,2]))                          # column 2 = first capture group (carbons)
}

# 2) add Class and chain length
species_annot <- long_species %>%
  mutate(Class      = str_extract(Lipid, class_pat),
         chain_len  = vapply(Lipid, sum_carbons, numeric(1))) %>%
  filter(!is.na(Class), !is.na(chain_len))

# 3) ΔZ per lipid species
delta_by_species <- species_annot %>%
  group_by(Lipid, Class, chain_len, Condition) %>%
  summarise(z_mean = mean(z, na.rm = TRUE), .groups="drop") %>%
  tidyr::pivot_wider(names_from = Condition, values_from = z_mean) %>%
  mutate(delta_z = LowInput - Control) %>%
  filter(!is.na(delta_z), chain_len >= 10, chain_len <= 60)  # tidy range; tweak as needed


# Remove LPC, LPE, P, PG, PS
delta_by_species <- delta_by_species %>%
  filter(!Class %in% c("LPC", "LPE", "PG", "PA", "PS"))

# 4) Plot (overall)
quartz()
ggplot(delta_by_species, aes(x = chain_len, y = delta_z)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(alpha = 0.25, size = 1) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 0.9, colour = "blue") +
  labs(x = "Total chain length (sum of carbons)", y = "ΔZ (LowInput − Control)",
       title = "Condition effect vs total chain length (species-level)") +
  plot_theme

# Save the plot
ggsave("fig/main/Chain_length_condition.png", width = 10, height = 6)


# Optional: facet by class to see which families drive the trend
quartz()
ggplot(delta_by_species, aes(chain_len, delta_z)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 0.9, colour = "blue") +
  facet_wrap(~ Class, scales = "free_x") +
  labs(x = "Total chain length", y = "ΔZ (LowInput − Control)") +
  plot_theme

# Save the faceted plot
ggsave("fig/main/Chain_length_individual.png", width = 10, height = 6)
