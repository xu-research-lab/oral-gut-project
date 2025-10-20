
library(ggplot2)
library(dplyr)
library(tidyr)
library(scatterpie)
library(reshape2)
load("../data/Extended Data Fig/shared_trans_species_stat.RData")
df<-subset(shared_new7,!is.na(YES))
df<-df%>%mutate(shared.YES=oral.gut-(NO+YES))
df1<-df[,c(1:4,6,28,15:16)]
df1<-subset(df1,shared.YES>=0)
df1<-df1%>%mutate(oral_total=oral+shared.YES+NO+YES)

df1<-df1%>%mutate(trans_rate=YES/oral_total)

df1$X2<-gsub("_group","",df1$X2)
df1$X1<-gsub("_unclassified_SGB19850","",df1$X1)
df1$id<-paste(df1$X1,df1$X2)
df1_long<-melt(df1,id.vars = c("group","X2", "X1", "Diagnosis","oral_total","trans_rate"))
df1_long <- df1 %>%
  pivot_longer(
    cols = c(oral, shared.YES, NO, YES),
    names_to = "variable",
    values_to = "value"
  )

df1_long<-df1_long%>%mutate(percent=value/oral_total*100)


df1_long$variable<-gsub("oral","shared.NO",df1_long$variable)
df1_long<-df1_long%>%mutate(variable=case_when(variable=="YES"~"oral-gut transmission",
                                               variable=="NO"~"oral-gut non-transmission",
                                               variable=="shared.NO"~"oral-only",
                                               variable=="shared.YES"~"oral-gut shared (unresolved strain)"))
											   
df2_long<-subset(df1_long,id %in% counts$id)

# 首先创建包含GC、EC和CRC的"all"组数据
df_all <- df2_long %>%
  filter(Diagnosis %in% c("GC", "EC", "CRC")) %>%
  group_by(group, X2, X1, variable) %>%
  summarise(
    oral_total = sum(oral_total),
    value = sum(value),
    .groups = "drop"
  ) %>%
  mutate(
    Diagnosis = "all",
    percent = value / oral_total * 100,
    trans_rate = NA,  # 由于合并后trans_rate可能不适用，设为NA
    id = paste(X1, X2)  # 重新创建id列
  )

# 选择需要的列，确保列顺序与原始数据一致
df_all <- df_all %>% select(names(df2_long))

# 将新数据与原始数据合并
df_combined <- bind_rows(df2_long, df_all)

# 查看结果
df_combined$id<-factor(df_combined$id,levels = counts$id)

ggplot()+
     geom_col(data = subset(df2_long),aes(percent,id,fill=variable))+
     facet_grid(~Diagnosis+group)+xlab("proportion (%)")+ylab(NULL)+theme_bw()+
theme(axis.text.x = element_text(hjust = 1,angle = 45),axis.text = element_text(size = 7,color=counts$color),legend.position = "top")+
labs(fill="Oral strain status")
ggsave("../results/strain_shared_trans_proportion.pdf",height = 7,width = 12)

#save(df,df1,shared_new7,df1_long,file="non_pma_pie_matrix.RData")


###plot only non-PMA

df3_long<-merge(df2_long,stat1,by.x = c(2,1,4),by.y = 1:3,all = T)
df3_long$percent[is.na(df3_long$percent)&df3_long$oral_counts>0]<-100
df3_long$variable[is.na(df3_long$variable)]<-"oral-only"

df3_long$variable<-gsub("transmitted","transmission",df3_long$variable)

df3_long$variable<-factor(df3_long$variable,levels=c("oral-only","oral-gut shared (unresolved strain)","oral-gut transmission","oral-gut non-transmission"))

df3_long$Diagnosis<-factor(df3_long$Diagnosis,levels=c("CRC","GC","EC","RA","HC"))
tmp<-subset(df3_long,group!="PMA"&!is.na(X2))%>%group_by(X2,variable)%>%summarise(value=sum(percent))
tmp1<-subset(tmp,value<600)

ggplot()+
    geom_col(data = subset(df3_long,group!="PMA"&X2!="SGB8002"&Diagnosis!="Treated"&X2%in% unique(tmp1$X2)),aes(percent,id.y,fill=variable))+
    facet_grid(~Diagnosis)+
	xlab("individual proportion (%)")+
	ylab(NULL)+theme_bw()+
    theme(axis.text.x = element_text(hjust = 1,angle = 45),
	axis.text = element_text(size = 7),legend.position = "right")+
    labs(fill="Oral strain status")

ggsave("../results/strain_shared_trans_proportion_nonpma.pdf",height = 7,width = 10)


###plot PMA

tmp<-subset(df3_long,group=="PMA"&!is.na(X2))%>%group_by(X2,variable)%>%summarise(value=sum(percent))
tmp1<-subset(tmp,value<600)

ggplot()+
  geom_col(data = subset(df3_long,group=="PMA"&X2!="SGB8002"&Diagnosis!="Treated"&X2%in% unique(tmp1$X2)),aes(percent,id.y,fill=variable))+
  facet_grid(~Diagnosis)+
  xlab("individual proportion (%)")+
  ylab(NULL)+theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45),
        axis.text = element_text(size = 7),legend.position = "right")+
  labs(fill="Oral strain status")

ggsave("../results/strain_shared_trans_proportion_pma.pdf",height = 7,width = 10)