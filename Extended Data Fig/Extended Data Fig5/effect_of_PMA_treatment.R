####plot PMA 
library(reshape2)
setwd("/beegfs/huangxiaochang/Oral_gut_cancer/mpa4/PMA_treatment_effect")
SGB<-read.table("SGB_profile.tsv",header = T,sep="\t")
rownames(SGB)<-SGB[,1]
SGB<-SGB[,-1]
feces<-SGB[,grepl("Feces",colnames(SGB))]
saliva<-SGB[,grepl("Saliva",colnames(SGB))]

feces<-feces[,gsub("Saliva","Feces",colnames(saliva))]
saliva<-melt(as.matrix(saliva))
nonPMA<-nonPMA[,gsub("PMA_","",colnames(nonPMA))]
nonPMA<-nonPMA[,gsub("PMA_","",colnames(PMA))]
nonPMA<-melt(as.matrix(nonPMA))
PMA<-melt(as.matrix(PMA))

PMA$Var2<-gsub("PMA_","",PMA$Var2)
###get SGB abundance table with PMA and nonPMA treatment sample paired
treatment<-merge(nonPMA,PMA,by=c(1,2),all = T)
colnames(treatment)<-c("sgb","id","nonpma","pma")
treatment<-subset(treatment,nonpma+pma>0)

#removed<-subset(treatment,nonpma>0&pma==0)

#added<-subset(treatment,nonpma==0&pma>0)

#removed$site[grepl("Saliva",removed$id)]<-"Saliva"
#removed$site[grepl("Feces",removed$id)]<-"Feces"
#added$site[grepl("Feces",added$id)]<-"Feces"
#added$site[grepl("Saliva",added$id)]<-"Saliva"
#removed$site<-as.factor(removed$site)

#added$group<-"PMA added"
#removed$group<-"PMA removed"
#added$nonpma<-added$pma

####combined the PMA-removed and PMA-add species
#df<-rbind(removed,added)

#df$site<-gsub("Feces","Gut",df$site)
#df$site<-gsub("Saliva","Oral",df$site)

#p1<-ggplot(df,aes(nonpma,color=site,linetype=group))+geom_density()+xlim(c(0,0.075))+theme_bw()+xlab(label="SGB relative abundance (%)")



###draw species abundance correlation analysis for paired nonPMA and PMA group
library(rstatix)
library(cowplot)
treatment$site[grepl("Saliva",treatment$id)]<-"Saliva"
treatment$site[grepl("Feces",treatment$id)]<-"Feces"

treatment$site<-as.factor(treatment$site)
species_number<-aggregate(treatment$sgb,by=list(treatment$site,treatment$sgb),length)
species_number<-subset(species_number,x>=20)
filtered<-merge(species_number[,1:2],treatment,by.x=c(1,2),by.y=c(5,1),all.x = T)

####correlation analysis for each pair of data
stat<-filtered%>%
  group_by(Group.1,Group.2)%>%
  cor_test(nonpma,pma,method = "pearson",use = "pairwise.complete.obs")%>%
  adjust_pvalue(method = "BH")


stat1<-subset(stat,p.adj<=0.05)
table(stat1$Group.1)
## plot distribution of pearson correlation coefficients
p2<-ggplot(stat,aes(cor,color=site))+geom_density()+theme_bw()+xlab(label="pearson correlation coefficients")+geom_vline(xintercept = 0.6,linetype="dashed")
#plot_grid(p1,p2,ncol = 2,labels = c("A","B"),rel_widths = c(1.2,1))

ggsave("effect of PMA treatment.pdf",width = 8,height = 4)