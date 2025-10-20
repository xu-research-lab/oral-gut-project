load("../data/Extended Data Fig/Extended Data Fig5/Extended_Data_Fig5e.RData")
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(aplot)
data51 <- data11[,c("id","Oral.gut.shared.species.counts","treatment","Diagnosis")]
data52 <- data51 %>%
  group_by(id, Diagnosis) %>%
  spread(key = treatment, value = Oral.gut.shared.species.counts) %>%
  ungroup()
gut_nonpma3<-ggplot(subset(data52),aes(`non-PMA`,Diagnosis))+geom_boxplot(aes(color=Diagnosis),show.legend = F,outlier.shape = NA)+
  #geom_point(aes(color=Diagnosis),show.legend = F)+
  geom_jitter(aes(color=Diagnosis),height = 0.1,width = 0,show.legend = F)+
  theme_classic()+
  theme(panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x  = element_blank(),axis.text.y = element_text(size = 6))+
  ylab(NULL)+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  stat_compare_means(comparisons = list(c("CRC","GC")),label = "p.signif",tip.length = 0.01,label.y = c(72.8,75.8,78.8),size = 2,vjust = 0.5)
gut_nonpma3
gut_pma3<-ggplot(subset(data52),aes(Diagnosis,PMA))+
  geom_boxplot(aes(,color=Diagnosis),show.legend = F,outlier.shape = NA)+
  #geom_point(aes(color=Diagnosis),show.legend = F)+
  geom_jitter(aes(color=Diagnosis),height = 0,width = 0.1,show.legend = F)+
  theme_classic()+theme(panel.grid = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size = 6))+
  ylab(NULL)+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  stat_compare_means(comparisons = list(c("GC", "CRC"),c("EC","GC")),label = "p.signif",tip.length = 0.01,label.y = c(72.8,75.8),size = 2,vjust = 0.5)
gut_pma3
gut3<- ggplot(subset(data52),aes(`non-PMA`,PMA,color=Diagnosis))+
  geom_point(aes(color=Diagnosis))+theme_bw()+
  stat_cor(show.legend = F,size=2.5)+
  stat_smooth(method = "lm",se = F,linewidth=0.5)+
  #theme(axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 6))+
  xlab(label="Shared SGB counts\n(non-PMA)")+
  ylab(label="Shared SGB counts\n(PMA)")+
  #facet_grid(~treatment)+
  scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "#D3D3D3",size=0.5)+xlim(0, 85)+ylim(0, 85)
gut3
gut_all3<-gut3%>%insert_right(gut_pma3,0.5)%>%insert_top(gut_nonpma3,0.5)
gut_all3
ggsave(filename = "../results/Extended_Data_Fig5e.pdf",plot = gut_all3,width = 5,height = 4)