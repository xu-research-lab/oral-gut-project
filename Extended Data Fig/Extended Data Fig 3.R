library(dplyr)
library(purrr)
library(ggplot2)
library(aplot)
library(ggpubr)
library(reshape2)
library(patchwork)
library(stringr)

load("Extended Data Fig 3.RData")

color_mapping <- c(
  "between individuals" = "#F8766D",
  "oral-gut within individual" = "#00BA38",
  "with-without PMA treatments" = "#619CFF"
)


p_density <- ggplot(dis_merged_filtered, aes(centered.min, fill = group.x)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = color_mapping) +
  geom_vline(aes(xintercept = thresholds)) +
  facet_wrap(~SBGS, scale = "free") +
  xlim(c(0, 2.5)) +
  theme_bw() +
  xlab(label = "Centered nGD") +
  ylab(label = "Density") +
  labs(fill = "Group") +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 6.8)
  )


ggsave("Extended Data Fig 3.png", plot = p_density, height = 15, width = 15)
ggsave("Extended Data Fig 3.pdf", plot = p_density, height = 15, width = 15)
