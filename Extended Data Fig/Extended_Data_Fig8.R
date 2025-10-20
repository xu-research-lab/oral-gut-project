library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
gene<-read.table("../data/Extended Data Fig/SCFA_genes.txt",header=T,row.names = NULL)
kinase<-subset(gene,enzyme=="kinase")
kinase_melt<-melt(kinase)
kinase_melt$Species<- gsub("s__", "", kinase_melt$Species)
kinase_melt_nonpma<-subset(kinase_melt,grepl("Feces_",variable)&!grepl("PMA_",variable))

kinase_melt_nonpma<-kinase_melt_nonpma %>%
  group_by(Species,variable) %>%
  summarise(value_sum=sum(value))

kinase_mean <- kinase_melt_nonpma %>%
  group_by(Species)%>%
  summarise(average = mean(value_sum))%>%
  arrange(desc(average)) %>%  # 先按平均值降序排序
  mutate(Species= factor(Species, levels = Species))  # 设置因子水平为排序后的顺序
kinase_mean<-subset(kinase_mean,Species!="unclassified")

unselected_mean<-kinase_mean[-c(1:10),]
unselected_mean$group<-"Average of other Species"
colnames(unselected_mean)[2]<-"value"

selected<-subset(kinase_melt_nonpma,Species%in% kinase_mean$Species[1:10])

selected<-selected[,c("Species","value_sum")]
selected$group<-selected$Species
colnames(selected)<-c("Species","value","group")

plot_data<-rbind(selected,unselected_mean)

med<-plot_data%>%
  filter(value>0)%>%
  group_by(group)%>%
  summarise(med=median(value))%>%
  arrange(desc(med)) %>%  # 先按平均值降序排序
  mutate(group = factor(group, levels = group)) 

plot_data$group<-factor(plot_data$group,levels = med$group)

p_butyrate<-ggplot(subset(plot_data,value>0), aes(value, group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(width = 0.0005, height = 0.1), shape = 1) +
    theme_bw() +
    scale_y_discrete(limits = rev(levels(plot_data$group))) +
    scale_x_log10() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "abundance", y = NULL, title = "Butyrate kinase contributor")
print(p_butyrate)
ggsave("../results/butyrate kinase contributor.pdf",height = 5,width = 5)
