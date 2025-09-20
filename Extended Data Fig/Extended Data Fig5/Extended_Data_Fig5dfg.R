load("Extended_Data_Fig5defg.RData")
library(ggplot2)
library(ggpubr)
library(aplot)
library(ggExtra)
library(dplyr)
library(stringr)
library(tidyr)
features_merge1$Diagnosis <- factor(features_merge1$Diagnosis,levels=c("CRC","GC","EC"))
#gut
data11 <- features_merge1[,c("id","gut","treatment","Diagnosis")]
data12 <- data11 %>%
  group_by(id, Diagnosis) %>%
  spread(key = treatment, value = gut) %>%
  ungroup()
select <- subset(data11,id=="Feces_1345909")
gut_nonpma<-ggplot(subset(data12),aes(`non-PMA`,Diagnosis))+geom_boxplot(aes(color=Diagnosis),show.legend = F,outlier.shape = NA)+
  #geom_point(aes(color=Diagnosis),show.legend = F)+
  geom_jitter(aes(color=Diagnosis),height = 0.1,width = 0,show.legend = F)+
  theme_classic()+
  theme(panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x  = element_blank(),axis.text.y = element_text(size = 6))+
  ylab(NULL)+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  stat_compare_means(comparisons = list(),label = "p.signif",tip.length = 0.01,label.y = c(32.8,34.8,36.8),size = 2,vjust = 0.5)
gut_nonpma
gut_pma<-ggplot(subset(data12),aes(Diagnosis,PMA))+
  geom_boxplot(aes(,color=Diagnosis),show.legend = F,outlier.shape = NA)+
  #geom_point(aes(color=Diagnosis),show.legend = F)+
  geom_jitter(aes(color=Diagnosis),width = 0.1,height = 0,show.legend = F)+
  theme_classic()+theme(panel.grid = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size = 6))+
  ylab(NULL)+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  stat_compare_means(comparisons = list(),label = "p.signif",tip.length = 0.01,label.y = c(32.8,34.8,36.8),size = 2,vjust = 0.5)
gut_pma
gut<- ggplot(subset(data12),aes(`non-PMA`,PMA,color=Diagnosis))+
  geom_point(aes(color=Diagnosis))+theme_bw()+
  stat_cor(show.legend = F,size=2.5)+
  stat_smooth(method = "lm",se = F,linewidth=0.5)+
  #theme(axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 6))+
  xlab(label="Transmitted abundance\n(non-PMA)")+
  ylab(label="Transmitted abundance\n(PMA)")+
  #facet_grid(~treatment)+
  scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "#D3D3D3",size=0.5)+xlim(0, 55)+ylim(0, 55)
gut
gut_all<-gut%>%insert_right(gut_pma,0.5)%>%insert_top(gut_nonpma,0.5)
gut_all
ggsave(filename = "Extended_Data_Fig5f.pdf",plot = gut_all,width = 5,height = 4)
#counts
data31 <- features_merge1[,c("id","counts","treatment","Diagnosis")]
data32 <- data31 %>%
  group_by(id, Diagnosis) %>%
  spread(key = treatment, value = counts) %>%
  ungroup()
select <- subset(data11,id=="Feces_1345909")
gut_nonpma2<-ggplot(subset(data32),aes(`non-PMA`,Diagnosis))+geom_boxplot(aes(color=Diagnosis),show.legend = F,outlier.shape = NA)+
  #geom_point(aes(color=Diagnosis),show.legend = F)+
  geom_jitter(aes(color=Diagnosis),height = 0.1,width = 0,show.legend = F)+
  theme_classic()+
  theme(panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x  = element_blank(),axis.text.y = element_text(size = 6))+
  ylab(NULL)+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  stat_compare_means(comparisons = list(),label = "p.signif",tip.length = 0.01,label.y = c(32.8,34.8,36.8),size = 2,vjust = 0.5)
gut_nonpma2
gut_pma2<-ggplot(subset(data32),aes(Diagnosis,PMA))+
  geom_boxplot(aes(,color=Diagnosis),show.legend = F,outlier.shape = NA)+
  #geom_point(aes(color=Diagnosis),show.legend = F)+
  geom_jitter(aes(color=Diagnosis),width = 0.1,height = 0,show.legend = F)+
  theme_classic()+theme(panel.grid = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size = 6))+
  ylab(NULL)+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  stat_compare_means(comparisons = list(),label = "p.signif",tip.length = 0.01,label.y = c(32.8,34.8,36.8),size = 2,vjust = 0.5)
gut_pma2
gut2<- ggplot(subset(data32),aes(`non-PMA`,PMA,color=Diagnosis))+
  geom_point(aes(color=Diagnosis))+theme_bw()+
  stat_cor(show.legend = F,size=2.5)+
  stat_smooth(method = "lm",se = F,linewidth=0.5)+
  #theme(axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 6))+
  xlab(label="Transmitted strain counts\n(non-PMA)")+
  ylab(label="Transmitted strain counts\n(PMA)")+
  #facet_grid(~treatment)+
  scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "#D3D3D3",size=0.5)+xlim(0, 16)+ylim(0, 16)
gut2
gut_all2<-gut2%>%insert_right(gut_pma2,0.5)%>%insert_top(gut_nonpma2,0.5)
gut_all2
ggsave(filename = "Extended_Data_Fig5g.pdf",plot = gut_all2,width = 5,height = 4)
#shared
data41 <- features_merge1[,c("id","shared_abundance.Feces","treatment","Diagnosis")]
data42 <- data41 %>%
  group_by(id, Diagnosis) %>%
  spread(key = treatment, value = shared_abundance.Feces) %>%
  ungroup()
select <- subset(data11,id=="Feces_1345909")
gut_nonpma3<-ggplot(subset(data42),aes(`non-PMA`,Diagnosis))+geom_boxplot(aes(color=Diagnosis),show.legend = F,outlier.shape = NA)+
  #geom_point(aes(color=Diagnosis),show.legend = F)+
  geom_jitter(aes(color=Diagnosis),height = 0.1,width = 0,show.legend = F)+
  theme_classic()+
  theme(panel.grid = element_blank(),axis.text.x = element_blank(),axis.ticks.x  = element_blank(),axis.text.y = element_text(size = 6))+
  ylab(NULL)+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  stat_compare_means(comparisons = list(),label = "p.signif",tip.length = 0.01,label.y = c(32.8,34.8,36.8),size = 2,vjust = 0.5)
gut_nonpma3
gut_pma3<-ggplot(subset(data42),aes(Diagnosis,PMA))+
  geom_boxplot(aes(,color=Diagnosis),show.legend = F,outlier.shape = NA)+
  #geom_point(aes(color=Diagnosis),show.legend = F)+
  geom_jitter(aes(color=Diagnosis),width = 0.1,height = 0,show.legend = F)+
  theme_classic()+theme(panel.grid = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.text.x = element_text(size = 6))+
  ylab(NULL)+xlab(NULL)+scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  stat_compare_means(comparisons = list(),label = "p.signif",tip.length = 0.01,label.y = c(32.8,34.8,36.8),size = 2,vjust = 0.5)
gut_pma3
gut3<- ggplot(subset(data42),aes(`non-PMA`,PMA,color=Diagnosis))+
  geom_point(aes(color=Diagnosis))+theme_bw()+
  stat_cor(show.legend = F,size=2.5)+
  stat_smooth(method = "lm",se = F,linewidth=0.5)+
  #theme(axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 6))+
  xlab(label="cumulative abundance\n(non-PMA)")+
  ylab(label="cumulative abundance\n(PMA)")+
  #facet_grid(~treatment)+
  scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cd", "HC"="#a6d854"))+
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "#D3D3D3",size=0.5)+xlim(0, 60)+ylim(0, 60)
gut3
gut_all3<-gut3%>%insert_right(gut_pma3,0.5)%>%insert_top(gut_nonpma3,0.5)
gut_all3
ggsave(filename = "Extended_Data_Fig5d.pdf",plot = gut_all3,width = 5,height = 4)
