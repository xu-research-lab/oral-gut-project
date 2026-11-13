load("../data/Extended Data Fig/Extended_Data_Fig1/shared_new9.RData")
library(ggplot2)
library(ggpubr)
shared_new9$Diagnosis <- factor(shared_new9$Diagnosis, levels = c("CRC", "GC", "EC","RA","HC"))
p1 <- ggplot(subset(shared_new9,group=="non-PMA"),aes(x=rate,y=oral_median,color=Diagnosis))+
  geom_point()+theme_bw()+scale_y_log10()+
  stat_smooth(method = "lm",se=F,linewidth=0.5)+stat_cor(method="spearman",size=1.8,label.x = 0.6,label.y.npc = "bottom")+
  labs(x="Shared rate",y="Oral median abundance")+
  scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))
p1
ggsave("../results/rate_oral_median.pdf",plot=p1,height = 2.6,width = 3.75)