library(dplyr)
library(purrr)
library(magrittr)
library(ggplot2)
library(aplot)
library(ggpubr)
library(scales)

load("Figure3a.RData")
head(trans2)


trans2 <- trans2 %>%
  mutate(Diagnosis.x = if_else(Diagnosis.x == "RA untreated", "RA", Diagnosis.x))
unique(trans2$Diagnosis.x)

diagnosis_levels <- c("CRC", "GC", "EC", "RA", "HC")


diagnosis_colors <- c(
  "CRC" = "#66c2a5",
  "GC" = "#fc8d62",
  "EC" = "#e78ac3",
  "RA" = "#8da0cb",
  "HC" = "#a6d854"
)

my_comparisons <- list(c("HC", "GC"), c("HC", "EC"), c("HC", "CRC"), c("RA", "HC"))


non_pma_data <- trans2 %>%
  filter(
    Diagnosis.x %in% diagnosis_levels,
    treatment.y == "non-PMA",
    !is.na(centered.min)
  ) %>%
  mutate(Diagnosis.x = factor(Diagnosis.x, levels = diagnosis_levels, ordered = T))
unique(non_pma_data$Diagnosis.x)
nGD_non_pma_plot <- ggplot(non_pma_data, aes(x = Diagnosis.x, y = centered.min)) +
  geom_boxplot(aes(color = Diagnosis.x), show.legend = FALSE, outlier.shape = NA) +
  geom_point(aes(color = Diagnosis.x), position = position_jitter(width = 0.1), size = 1, alpha = 0.7, show.legend = FALSE) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = label_scientific()
  ) +
  theme_bw() +
  scale_color_manual(values = diagnosis_colors) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", vjust = 0.5, hide.ns = TRUE) +
  labs(
    y = "centered nGD",
    x = NULL
  )



print(nGD_non_pma_plot)

ggsave("Figure3a.png", plot = nGD_non_pma_plot, height = 2.5, width = 2.8)
ggsave("Figure3a.pdf", plot = nGD_non_pma_plot, height = 2.5, width = 3.2)
