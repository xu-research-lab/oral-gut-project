library(Maaslin2)

load("data_for_GI_gut_cancer.RData")

otu_table2<-data.frame(t(otu_table1))
otu_table2<-otu_table2[rownames(sam1),]

fit_big_cohort_all = Maaslin2(
  input_data = otu_table2,
  input_metadata = sam1,
  output = "Maaslin2_big_cohort_all",
  fixed_effects = "Response_6",random_effects = c("Diagnosis","Age","BMI_category"),normalization = "NONE")

sam_CRC<-subset(sam,Diagnosis=="CRC"&grepl("ICI",Treatment_class)&Response_6 %in% c("R","NR"))

fit_big_cohort_CRC = Maaslin2(
  input_data = otu_table2,
  input_metadata = sam_CRC,
  output = "Maaslin2_big_cohort_CRC",
  fixed_effects = "Response_6",random_effects = c("Diagnosis","Age","BMI_category"),normalization = "NONE")
  
sam_EC<-subset(sam,Diagnosis=="EC"&grepl("ICI",Treatment_class)&Response_6 %in% c("R","NR"))

fit_big_cohort_EC = Maaslin2(
  input_data = otu_table2,
  input_metadata = sam_EC,
  output = "Maaslin2_big_cohort_EC",
  fixed_effects = "Response_6",random_effects = c("Diagnosis","Age","BMI_category"),normalization = "NONE")
  
fit_big_cohort_GC = Maaslin2(
  input_data = otu_table2,
  input_metadata = sam_GC,
  output = "Maaslin2_big_cohort_GC",
  fixed_effects = "Response_6",random_effects = c("Diagnosis","Age","BMI_category"),normalization = "NONE")
  
  
###manually combine output files of "Maaslin2_big_cohort_all","Maaslin2_big_cohort_CRC","Maaslin2_big_cohort_GC"and "Maaslin2_big_cohort_EC"

library(dplyr)
library(microeco)
library(magrittr)
library(ggplot2)
library(aplot)
load("Figure1e.RData")

load("masslin2_result.RData")
plot_data_p1_big <- subset(result2, !grepl("GGB", id))


plot_data_p1_big$group.x <- factor(plot_data_p1_big$group.x, levels = c("All", "CRC", "GC", "EC"), ordered = TRUE)



id_to_group_map <- plot_data_p1_big %>%
  select(id, group.y) %>%
  distinct() %>%
  arrange(id)

y_colors_p1_big <- id_to_group_map %>%
  mutate(label_color = ifelse(group.y == "shared", "red", "black")) %>%
  pull(label_color)



p1_big_final <- ggplot(plot_data_p1_big, aes(group.x, id, fill = coef)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  geom_text(aes(label = sig), vjust = 0.7, size = 5) +
  theme_bw() +
  labs(x = NULL, y = NULL) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8, colour = y_colors_p1_big),
  )


ggsave("Figure1e.pdf", plot = p1_big_final, height = 7, width = 5)