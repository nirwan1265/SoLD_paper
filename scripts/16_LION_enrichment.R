library(tidyverse)
# your data frame
df <- read.csv("table/Linex2/LION-enrichment.csv")
colnames(df) <- c("TermID", "Description", "Annotated", "p.value", "q.value")



df <- df %>%
  # sort by ascending q‐value (smallest = most significant)
  arrange(q.value) %>%
  # turn Description into a factor in that order, then reverse so the top row is most significant
  mutate(
    Description = fct_rev(fct_inorder(Description)),
    logQ = -log10(q.value),
    logP = -log10(p.value)
  )

# plot
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
quartz()
enrichment <- ggplot(df2, aes(x = logQ, y = Description)) +
  # vertical line at FDR q = 0.05
  geom_vline(xintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_point(aes(size = Annotated, color = logP)) +
  scale_color_viridis_c(name = expression(-log[10](p~value)), option = "C") +
  scale_size_continuous(name = "Annotated hits") +
  labs(
    x = expression(-log[10](FDR~q~value)),
    y = NULL,
    title = "LION Enrichment"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.key.width    = unit(1.2, "cm")
  ) + nature_theme

# Save the plot
ggsave("fig/main/Fig7_LION_enrichment.png", enrichment, width = 8, height = 12, dpi = 300, bg = "white")
