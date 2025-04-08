library(tidyverse)
library(gghalves)
library(ggplot2)
library(ggpubr)
load("SGB_counts_abundance.RData")
total <- ggplot(SGB_counts_abundance,aes(x = Diagnosis, y = `total observed`, color = Response_6)) +
  # 左半边小提琴图（NR组）
  geom_half_violin(
    data = ~ subset(.x, Response_6 == levels(factor(.x$Response_6))[1]),  # 修正筛选方式
    fill = NA, side = "l", show.legend = FALSE, linewidth = 0.2
  ) +
  # 右半边小提琴图（R组）
  geom_half_violin(
    data = ~ subset(.x, Response_6 == levels(factor(.x$Response_6))[2]),
    fill = NA, side = "r", show.legend = FALSE, linewidth = 0.2
  ) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.3, seed = 123),size=1,alpha=0.7,show.legend = F)+
  geom_half_boxplot(data = filter(SGB_counts_abundance,Response_6 == "NR"),aes(x = Diagnosis, y = `total observed`), width = 0.25, side = "l",fill = "#E41A1C",coef = 0,alpha = 0.5,outlier.shape = NA, linewidth = 0.2) +
  geom_half_boxplot(data = filter(SGB_counts_abundance, Response_6 == "R"),aes(x = Diagnosis, y = `total observed`), width = 0.25, side = "r",fill = "#377EB8",coef = 0,alpha = 0.5,outlier.shape = NA, linewidth = 0.2) +
  theme_bw()+
  stat_compare_means(
    comparisons = list(c("CRC", "GC"), c("CRC", "EC")),
    label = "p.signif", vjust = 0.5, hide.ns = TRUE
  ) +
  stat_compare_means(
    method = "wilcox.test", label = "p.signif",
    label.y.npc = "bottom", show.legend = FALSE,label.y = c(-10, -10, -10)
  ) +
  theme(
    strip.background = element_blank(),plot.title = element_text(size = 12, hjust = 0.5),axis.text.x = element_text(size = 12)
  ) +
  ylab("SGB counts") +
  xlab(NULL) +
  scale_color_manual(values = c("#E41A1C", "#377EB8"), name = "Response") +
  labs(title = "Total") + theme(plot.title = element_text(hjust = 0.5),legend.position = "none")+
  scale_y_continuous(breaks = c(200, 400, 600))


putative_oral <- ggplot(SGB_counts_abundance,aes(x = Diagnosis, y = `oral-origin`, color = Response_6)) +
  # 左半边小提琴图（NR组）
  geom_half_violin(
    data = ~ subset(.x, Response_6 == levels(factor(.x$Response_6))[1]),  # 修正筛选方式
    fill = NA, side = "l", show.legend = FALSE, linewidth = 0.2
  ) +
  # 右半边小提琴图（R组）
  geom_half_violin(
    data = ~ subset(.x, Response_6 == levels(factor(.x$Response_6))[2]),
    fill = NA, side = "r", show.legend = FALSE, linewidth = 0.2
  ) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.3, seed = 123),size=1,alpha=0.7,show.legend = F)+
  geom_half_boxplot(data = filter(SGB_counts_abundance, Response_6 == "NR"),aes(x = Diagnosis, y = `oral-origin`), width = 0.25, side = "l",fill = "#E41A1C",coef = 0,alpha = 0.5,outlier.shape = NA, linewidth = 0.2) +
  geom_half_boxplot(data = filter(SGB_counts_abundance, Response_6 == "R"),aes(x = Diagnosis, y = `oral-origin`), width = 0.25, side = "r",fill = "#377EB8",coef = 0,alpha = 0.5,outlier.shape = NA, linewidth = 0.2) +
  theme_bw()+
  stat_compare_means(
    comparisons = list(c("CRC", "GC"), c("GC", "EC"), c("CRC", "EC")),
    label = "p.signif", vjust = 0.5, hide.ns = TRUE
  ) +
  stat_compare_means(
    method = "wilcox.test", label = "p.signif",
    label.y.npc = "bottom", show.legend = FALSE,label.y = c(-10, -10, -10)
  ) +
  theme(
    strip.background = element_blank(),plot.title = element_text(size = 12, hjust = 0.5),axis.text.x = element_text(size = 12)
  ) +
  ylab("SGB counts") +
  xlab(NULL) +
  scale_color_manual(values = c("#E41A1C", "#377EB8"), name = "Response") +
  labs(title = "Putative oral") + theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
  
  
putative_oral                                          
oral_origin <- ggplot(SGB_counts_abundance,aes(x = Diagnosis, y = `cumlative abundance`*100, color = Response_6)) +
  # 左半边小提琴图（NR组）
  geom_half_violin(
    data = ~ subset(.x, Response_6 == levels(factor(.x$Response_6))[1]),  # 修正筛选方式
    fill = NA, side = "l", show.legend = FALSE, linewidth = 0.2
  ) +
  # 右半边小提琴图（R组）
  geom_half_violin(
    data = ~ subset(.x, Response_6 == levels(factor(.x$Response_6))[2]),
    fill = NA, side = "r", show.legend = FALSE, linewidth = 0.2
  ) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.3, seed = 123),size=1,alpha=0.7,show.legend = F)+
  geom_half_boxplot(data = filter(subset(SGB_counts_abundance, Response_6 == "NR")),aes(x = Diagnosis, y = `cumlative abundance`*100), width = 0.25, side = "l",fill = "#E41A1C",coef = 0,alpha = 0.5,outlier.shape = NA, linewidth = 0.2) +
  geom_half_boxplot(data = filter(subset(SGB_counts_abundance, Response_6 == "R")),aes(x = Diagnosis, y = `cumlative abundance`*100), width = 0.25, side = "r",fill = "#377EB8",coef = 0,alpha = 0.5,outlier.shape = NA, linewidth = 0.2) +
  theme_bw()+
  stat_compare_means(
    comparisons = list(c("CRC", "GC"), c("GC", "EC")),
    label = "p.signif", vjust = 0.5, hide.ns = TRUE
  ) +
  stat_compare_means(
    method = "wilcox.test", label = "p.signif",
    label.y.npc = "bottom", show.legend = FALSE,label.y = c(-3, -3, -3)
  ) +
  theme(
    strip.background = element_blank(),plot.title = element_text(size = 12, hjust = 0.5),axis.text.x = element_text(size = 12)
  ) +
  ylab("Cumulative abundance") +
  xlab(NULL) +
  scale_color_manual(values = c("#E41A1C", "#377EB8"), name = "Response") + scale_y_log10(limits = c(1e-03, 1e+03)) +
  labs(title = "Putative oral") + theme(plot.title = element_text(hjust = 0.5)) 
oral_origin
library(patchwork)
all <- total+putative_oral+oral_origin
all
ggsave("p_observed_shared_big_cohort_boxplot_violin.pdf",plot = all,height = 2.5,width = 8)



