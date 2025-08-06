library(ggplot2)
library(dplyr)
library(forcats)
library(tidyverse)


# your data frame
df <- read.csv("table/Linex2/LION-enrichment.csv")

# Filter p.value for <0.05
df <- df %>%
  filter(p.value < 0.05)


# add the –log10 metrics:
df2 <- df %>%
  mutate(
    logP = -log10(p.value),
    logQ = -log10(FDR.q.value),
    Description = fct_reorder(Discription, logP)     # so bigger dots are lower on the axis
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

lion_lipid <- ggplot(df2, aes(x=logQ, y=Description)) +
  # dashed line at FDR = 0.05 → –log10(0.05)
  geom_vline(xintercept = -log10(0.05), color="red", linetype="dashed") +
  # split into two facets by up/down
  facet_grid(Regulated ~ ., scales="free_y", space="free_y") +
  # dots sized by number of annotated lipids, colored by –log10(p)
  geom_point(aes(size=Annotated, color=logP)) +
  scale_color_viridis_c(
    name = expression(-log[10](p~value)),
    option = "C"
  ) +
  scale_size_continuous(name = "Annotated hits") +
  labs(
    x = expression(-log[10](FDR~q~value)),
    y = NULL,
    title = "LION enrichment analysis (ranking mode)\nlowinput vs. control"
  ) +
  nature_theme +
  theme(
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor   = element_blank()
  )

# Save the plot
ggsave("fig/main/Fig7_LION_enrichment.png",plot=last_plot(), width = 12, height = 14, dpi = 300, bg = "white")
