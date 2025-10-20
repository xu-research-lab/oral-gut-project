library(dplyr)
library(purrr)
library(magrittr)
library(ggplot2)
library(microeco)

load("../data/Figure1/data_for_GI_gut_cancer.RData")
rownames(sam1)<-sam1$X
sam1<-sam1[,-1]
otu_table2<-otu_table1[,rownames(sam1)]
otu_table2<-otu_table2*100
rownames(otu_table2)<-gsub("t__","",rownames(otu_table2))
tax_id1<-tax_id[,-9]
dataset<-microtable$new(sample_table=sam1,otu_table=otu_table2,tax_table = tax_id1)

lefse<-trans_diff$new(dataset=dataset,method="lefse",group="Diagnosis",alpha=0.05,taxa_level = "SGB",lefse_min_subsam = 5,lefse_norm = 10e+6)
res_diff<-lefse$res_diff
res_diff$SGB<-sub(".*?\\|SGB", "SGB",res_diff$Taxa)
res_diff$Taxa<-sub(".*?s__", "",res_diff$Taxa)
res_diff$Taxa<-sub("[|]", " ",res_diff$Taxa)

res_diff<-merge(res_diff,tax_id[,c("SGB","group")],by.x = 9,by.y = 1,all.x = T)

plot_data <- res_diff%>%
  filter(!grepl("GGB", Taxa) & LDA >= 2.5 & Group %in% c("CRC", "GC", "EC")) %>%
  mutate(label_color = ifelse(!is.na(group), "red", "black")) %>%
  mutate(Group = factor(Group, levels = c("CRC", "GC", "EC"))) %>%
  arrange(Group, LDA) %>%
  mutate(Taxa = factor(Taxa, levels = unique(Taxa)))

y_axis_colors <- plot_data$label_color

lefse_p_final <- ggplot(plot_data, aes(x = LDA, y = Taxa, fill = Group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    name = "Diagnosis",
    values = c("CRC" = "#66c2a5", "GC" = "#fc8d62", "EC" = "#e78ac3"),
    limits = c("CRC", "GC", "EC")
  ) +
  labs(x = "LDA Score", y = NULL) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8, colour = y_axis_colors),
    axis.ticks.x = element_line(colour = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave("../results/lefse_Diagnosis_big_cohort.pdf",height = 7,width = 5)