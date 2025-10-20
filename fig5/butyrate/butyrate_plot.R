library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

load("../data/Figure5/butyrate/butyrate kinase contributor.RData")# dataframe `gene` contains the identified key gene abundance involoved in SCFA production 
load("../data/Figure5/butyrate/scfa_gene_features.RData") ## dataframe `scfa_gene_features` alreadly merged the gene profiles, transmitted abundnce and metadata together 
kinase<-subset(gene,enzyme=="kinase")
kinase_melt<-melt(kinase)
kinase_melt$Species_cleaned <- gsub("s__", "", kinase_melt$Species_cleaned )
kinase_melt_nonpma<-subset(kinase_melt,grepl("Feces_",variable)&!grepl("PMA_",variable))

kinase_melt_nonpma<-kinase_melt_nonpma %>%
  group_by(Species_cleaned,variable) %>%
  summarise(value_sum=sum(value))

kinase_mean <- kinase_melt_nonpma %>%
  group_by(Species_cleaned)%>%
  summarise(average = mean(value_sum))%>%
  arrange(desc(average)) %>%  
  mutate(Species_cleaned = factor(Species_cleaned, levels = Species_cleaned))  
kinase_mean<-subset(kinase_mean,Species_cleaned!="unclassified")

unselected_mean<-kinase_mean[-c(1:10),]
unselected_mean$group<-"Average of other Species"
colnames(unselected_mean)[2]<-"value"

selected<-subset(kinase_melt_nonpma,Species_cleaned %in% kinase_mean$Species_cleaned[1:10])

selected<-selected[,c("Species_cleaned","value_sum")]
selected$group<-selected$Species_cleaned
colnames(selected)<-c("Species_cleaned","value","group")

plot_data<-rbind(selected,unselected_mean)
table(plot_data$group)

med<-plot_data%>%
  filter(value>0)%>%
  group_by(group)%>%
  summarise(med=median(value))%>%
  arrange(desc(med)) %>%  
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
