library(ggplot2)
library(tidyr)
library(dplyr)
library(aplot)
library(ggpubr)
gut_boxplot_shared<-ggplot(subset(result,treatment=="non-PMA"&Diagnosis!="Treated"),aes(Diagnosis,gut_sum))+
  geom_boxplot(aes(color=Diagnosis),show.legend = F,outlier.shape = NA)+
  #geom_point(aes(color=Diagnosis),width = 0.1,show.legend = F)+
  geom_jitter(aes(color=Diagnosis),width = 0.05,height = 0,show.legend = F)+
  theme_bw()+
  #theme(panel.grid = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size = 6))+
  ylab(label="Cumulative abundance in gut")+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8DA0CB", "HC"="#A6D854"))+
  stat_compare_means(comparisons = list(c("CRC", "HC"),c("EC", "HC"),c("GC", "HC"),c("CRC", "GC"),c("EC", "RA"),c("RA","GC")),label = "p.signif",tip.length = 0.01,label.y = c(42.8,45.8,48.8,51.8,54.8,57.8),size = 3,vjust = 0.5)
gut_boxplot_shared
ggsave(filename = "gut_boxplot_shared_new_20250819.pdf",plot = gut_boxplot_shared,width = 3,height = 2.2)

gut<- ggplot(subset(result,treatment=="non-PMA"&Diagnosis!="Treated"),aes(oral_sum,gut_sum,color=Diagnosis))+
  geom_point()+theme_bw()+
  stat_cor(show.legend = F,size=2.5,label.x = 0)+
  stat_smooth(method = "lm",se = F,linewidth=0.5)+
  #theme(axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 6))+
  xlab(label="Cumulative abundance in oral")+
  ylab(label="Cumulative abundance in gut")+
  #facet_grid(~treatment)+
  scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8DA0CB", "HC"="#A6D854"))
gut
ggsave(filename = "gut_new_20250819.pdf",plot = gut,width = 3.75,height = 2.6)
oral_gut_shared_SGB_counts<-ggplot(subset(result,treatment=="non-PMA"&Diagnosis!="Treated"),aes(Diagnosis,shared_counts))+
  geom_boxplot(aes(color=Diagnosis),show.legend = F,outlier.shape = NA)+
  #geom_point(aes(color=Diagnosis),width = 0.1,show.legend = F)+
  geom_jitter(aes(color=Diagnosis),width = 0.05,height = 0,show.legend = F)+
  theme_bw()+
  #theme(panel.grid = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size = 6))+
  ylab(label="Shared SGB counts")+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8DA0CB", "HC"="#A6D854"))+
  stat_compare_means(comparisons = list(c("HC", "RA"),c("EC", "HC"),c("GC", "HC"),c("CRC", "GC")),label = "p.signif",tip.length = 0.01,label.y = c(67.8,71.8,75.8,79.8),size = 3,vjust = 0.5)
oral_gut_shared_SGB_counts
ggsave(filename = "oral_gut_shared_SGB_counts_new_20250821.pdf",plot = oral_gut_shared_SGB_counts,width = 2.8,height = 2.2)
