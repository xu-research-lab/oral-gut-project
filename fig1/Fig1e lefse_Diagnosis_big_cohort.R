library(dplyr)
library(purrr)
library(magrittr)
library(microeco)

load("data_for_GI_gut_cancer.RData")
dataset<-microtable$new(sample_table=sam1,otu_table=otu_table1,tax_table = tax_id1)


lefse<-trans_diff$new(dataset=dataset,method="lefse",group="Diagnosis",alpha=0.05,taxa_level = "Strain",lefse_min_subsam = 5,lefse_norm = 10e+6)

group_p<-ggplot(subset(res_diff,!grepl("GGB",Taxa)&LDA>=2.5),aes("group",Taxa,fill=group))+geom_tile()+theme_bw()+
scale_x_discrete(expand = c(0,0))+xlab(NULL)+ylab(NULL)+
theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.text.y=element_text(size=8))

lefse_p<-ggplot(subset(res_diff,!grepl("GGB",Taxa)&LDA>=2.5),aes(LDA,Taxa,fill=Group))+geom_col(width = 0.7)+theme_bw()+
scale_x_continuous(expand = c(0,0))+ scale_fill_manual(values = c("#66c2a5", "#e78ac3", "#fc8d62"),name="Diagnosis")+
theme(axis.text.y=element_blank(),axis.ticks.y = element_blank(),panel.grid.major.x = element_blank())+
ylab(NULL)+xlab("LDA score")


lefse_p%>%insert_left(group_p,0.05)
ggsave("lefse_Diagnosis_big_cohort.pdf",height = 7,width = 5)

save(res_diff,dataset,otu_table1,sam1,lefse,lefse_p,group_p,file="lefse_Diagnosis_big_cohort.RData")