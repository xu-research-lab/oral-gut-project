library(dplyr)
library(purrr)
library(ggplot2)
library(aplot)
library(ggpubr)
library(reshape2)
library(patchwork)
library(scatterpie)

load("Extended Data Fig 7c.RData")

ggplot(subset(sub_otu, value > 0), aes(paste(V2, "\n", Var1), value * 100, color = Response_6)) +
    geom_boxplot() +
    stat_compare_means(label = "p.format", show.legend = F) +
    ylim(c(0, 1.6)) +
    scale_y_log10() +
    geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
    geom_text(aes(x = paste(V2, "\n", Var1), y = 1, label = sig), color = "black") +
    theme_bw() +
    scale_color_manual(values = c("#377EB8", "#E41A1C"), name = "Response") +
    theme(axis.text.y = element_text(hjust = 0, vjust = 0.5), strip.background = element_blank(), axis.text.x = element_text(size = 7), panel.grid = element_blank()) +
    xlab(NULL) +
    coord_flip() +
    ylab("relative abundance (%)")

ggsave("Extended Data Fig 7c.pdf", height = 2.5, width = 5)
