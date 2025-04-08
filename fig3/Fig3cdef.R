library(ggplot2)
library(ggpubr)
library(aplot)
library(ggExtra)
gut_boxplot_trans<-ggplot(subset(features,Diagnosis!="Treated"&treatment=="non-PMA"),aes(Diagnosis,transmitted.SGB.counts))+geom_boxplot(aes(color=Diagnosis),show.legend = F,outlier.shape = NA)+
  geom_jitter(aes(color=Diagnosis),width = 0.1,show.legend = F)+theme_bw()+
  #theme(panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x  = element_blank())+
  ylab(label="Transmitted SGB counts")+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  stat_compare_means(comparisons = list(c("EC", "HC"),c("GC", "HC"),c("CRC", "GC"),c("RA","GC")),label = "p.signif",tip.length = 0.01,label.y = c(11.8,12.8,13.8,14.8),size = 3.5,vjust = 0.5)
gut_boxplot_trans
ggsave(filename = "gut_boxplot_trans_new_20250820.pdf",plot = gut_boxplot_trans,width = 3,height = 2.5)

gut_boxplot_shared<-ggplot(subset(features,Diagnosis!="Treated"&treatment=="non-PMA"),aes(Diagnosis,Oral.gut.shared.species.counts))+
  geom_boxplot(aes(,color=Diagnosis),show.legend = F,outlier.shape = NA)+geom_jitter(aes(color=Diagnosis),width = 0.2,show.legend = F)+
  theme_classic()+theme(panel.grid = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())+
  ylab(NULL)+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  stat_compare_means(comparisons = list(c("HC", "RA"),c("EC", "HC"),c("GC", "HC"),c("CRC", "GC"),c("CRC", "RA")),label = "p.signif",tip.length = 0.01,label.y = c(55.8,59.8,63.8,67.8,71.8),size = 2)
gut_boxplot_shared
gut<- ggplot(subset(features,Diagnosis %in% c("CRC", "GC", "EC","RA","HC")&treatment=="non-PMA"),aes(transmitted.SGB.counts,Oral.gut.shared.species.counts,color=Diagnosis))+
  geom_point(aes(color=Diagnosis))+theme_bw()+
  stat_cor(show.legend = F,size=2.5)+
  stat_smooth(method = "lm",se = F,linewidth=0.5)+
  xlab(label="Transmitted SGB counts")+
  ylab(label="Shared SGB counts")+
  #facet_grid(~treatment)+
  scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "#D3D3D3")
gut
ggsave(filename = "gut_counts_new_20250819.pdf",plot = gut,width = 3.75,height = 2.6)
gut_trans_shared_counts<-gut%>%insert_right(gut_boxplot_shared,0.5)%>%insert_top(gut_boxplot_trans,0.5)
gut_trans_shared_counts
ggsave(filename = "final_color_five_withoutns_non_PMA_transmitted_abundance_vs_shared_abundance.pdf",plot = gut_trans_shared_counts,width = 5,height = 4)


gut_boxplot_trans_abundance<-ggplot(subset(features,Diagnosis!="Treated"&treatment=="non-PMA"),aes(Diagnosis,transmitted.abundance.in.Gut))+geom_boxplot(aes(color=Diagnosis),show.legend = F,outlier.shape = NA)+
  #geom_point(aes(color=Diagnosis),height = 0.2,show.legend = F)+
  geom_jitter(aes(color=Diagnosis),width = 0.1,show.legend = F)+
  theme_bw()+
  #theme(panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x  = element_blank(),axis.text.y = element_text(size = 6))+
  ylab(label="Transmitted abundance (gut)")+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  stat_compare_means(comparisons = list(c("CRC", "HC"),c("EC", "HC"),c("GC", "HC"),c("CRC", "GC"),c("EC", "RA"),c("RA","GC")),label = "p.signif",tip.length = 0.01,label.y = c(34.8,36.8,38.8,40.8,42.8,44.8),size = 3.5,vjust = 0.5)
gut_boxplot_trans_abundance
ggsave(filename = "gut_boxplot_trans_abundance_new_20250820.pdf",plot = gut_boxplot_trans_abundance,width = 3,height = 2.5)
gut_boxplot_shared_abundance<-ggplot(subset(features,Diagnosis!="Treated"&treatment=="non-PMA"),aes(Diagnosis,shared_abundance.Feces))+
  geom_boxplot(aes(,color=Diagnosis),show.legend = F,outlier.shape = NA)+
  #geom_point(aes(color=Diagnosis),width = 0.2,show.legend = F)+
  geom_jitter(aes(color=Diagnosis),width = 0.03,show.legend = F)+
  theme_classic()+theme(panel.grid = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size = 6))+
  ylab(NULL)+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  stat_compare_means(comparisons = list(c("HC", "RA"),c("EC", "HC"),c("GC", "HC"),c("CRC", "GC"),c("EC", "RA"),c("RA","GC")),label = "p.signif",tip.length = 0.01,label.y = c(38.8,41.8,44.8,47.8,50.8,53.8),size = 2,vjust = 0.5)
gut_boxplot_shared_abundance
gut<- ggplot(subset(features,Diagnosis %in% c("CRC", "GC", "EC","RA","HC")&treatment=="non-PMA"),aes(transmitted.abundance.in.Gut,shared_abundance.Feces,color=Diagnosis))+
  geom_point(aes(color=Diagnosis))+theme_bw()+
  stat_cor(show.legend = F,size=2.5)+
  stat_smooth(method = "lm",se = F,linewidth=0.5)+
  #theme(axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 6))+
  xlab(label="Transmitted abundance (gut)")+
  ylab(label="Shared abundance (gut)")+
  #facet_grid(~treatment)+
  scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "#D3D3D3",size=0.5)
gut
ggsave(filename = "gut_abundance_new_20250819.pdf",plot = gut,width = 3.75,height = 2.6)
gut_trans_shared_counts_abundance<-gut%>%insert_right(gut_boxplot_shared_abundance,0.5)%>%insert_top(gut_boxplot_trans_abundance,0.5)
gut_trans_shared_counts_abundance
ggsave(filename = "abundance_non-PMA_20250808.pdf",plot = gut_trans_shared_counts_abundance,width = 5,height = 4)

save(list=c("features"),file = "Fig3cdef.RData")

