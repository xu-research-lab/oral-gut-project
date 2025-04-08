library(dplyr)
library(tidyr)

load("butyrate kinase contributor.RData")# dataframe `gene` contains the identified key gene abundance involoved in SCFA production 
load("scfa_gene_features.RData") ## dataframe `scfa_gene_features` alreadly merged the gene profiles, transmitted abundnce and metadata together 
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
ggsave("butyrate kinase contributor.pdf",height = 5,width = 5)

save(plot_data,p_butyrate,gene,kinase,kinase_mean,kinase_melt,kinase_melt_nonpma,file = "butyrate kinase contributor.RData")


###plot butyrate kinase association with transmitted abundance
ggplot(subset(enzyme_features,!is.na(treatment)&group=="Gut"&enzyme=="Butyrate kinase"),aes(Feces,transmitted.abundance.in.Gut))+geom_point(shape=1)+
       stat_smooth(method = "lm",se = F)+facet_wrap(~treatment)+
       stat_cor(method = "spearman",color="blue")+theme_bw()+
       theme(strip.background = element_blank(),panel.grid = element_blank())+
       scale_y_continuous(expand = c(0,0))+xlab("Abundance of butyrate kinase")+ylab("Transmitted abundance")
	   
	    ggsave("butyrate_kinase.pdf",height = 2.5,width = 5)