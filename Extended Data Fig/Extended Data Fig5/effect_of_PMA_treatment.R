library(reshape2)
library(tidyverse)
library(rstatix)
library(ggplot2)
SGB<-read.table("../data/Extended Data Fig/Extended Data Fig5/SGB_profile.tsv",header = T,sep="\t")
SGB<-melt(SGB)
SGB$variable<-gsub("d_metagenome","",SGB$variable)
SGB<-SGB%>%mutate(site=case_when(grepl("Feces",variable)~"Feces",
                                 grepl("Saliva",variable)~"Saliva"),
                  treatment=case_when(grepl("PMA_",variable)~"PMA",
                                 !grepl("PMA_",variable)~"non-PMA"))
##remove samples from public cohort
SGB<-subset(x = SGB,!grepl("Control|reated",variable))

SGB$variable<-gsub("PMA_","",SGB$variable)

treatment<-SGB%>%pivot_wider(names_from = treatment,values_from = value)

colnames(treatment)<-c("sgb","id","site","nonpma","pma")
treatment$nonpma[is.na(treatment$nonpma)]<-0
treatment$pma[is.na(treatment$pma)]<-0
treatment<-subset(treatment,nonpma+pma>0)


###draw species abundance correlation analysis for paired non-PMA and PMA group

species_number<-treatment%>%
  group_by(sgb,site)%>%
  summarise(number=n())%>%
  subset(number>=20)

filtered<-treatment %>%
  left_join(species_number, by = c("sgb", "site")) %>%
  filter(!is.na(number))

####correlation analysis for each pair of data
stat<-filtered%>%
  group_by(sgb,site)%>%
  cor_test(nonpma,pma,method = "pearson",use = "pairwise.complete.obs")%>%
  adjust_pvalue(method = "BH")

stat1<-subset(stat,p.adj<=0.05)
## plot distribution of pearson correlation coefficients
p<-ggplot(stat,aes(cor,color=site))+geom_density()+theme_bw()+xlab(label="pearson correlation coefficients")+geom_vline(xintercept = 0.6,linetype="dashed")
#plot_grid(p1,p2,ncol = 2,labels = c("A","B"),rel_widths = c(1.2,1))
p
ggsave("../results/effect of PMA treatment.pdf",width = 8,height = 4)