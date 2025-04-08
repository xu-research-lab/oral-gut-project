library(ggplot2)
library(ggpubr)
load("bugbase.RData")

bugbase$Diagnosis<-factor(bugbase$Diagnosis,levels = c("CRC","GC","EC"))

bugbase_plot<-ggplot(subset(bugbase,site=="gut"),aes(ratio,transmitted.abundance.in.Gut,color=Diagnosis))+
  geom_point(size=1)+stat_smooth(method = "lm",se=F,size=1)+stat_cor(method = "spearman",size=2)+theme_bw()+
  facet_wrap(~treatment,scale="free")+scale_y_log10()+scale_x_log10()+
  scale_color_discrete(type = c("#66c2a5","#fc8d62","#e78ac3"))+ylab(label="transmssion abundance of gut")+
  xlab(label="relative abundance")+
  theme(axis.text = element_text(size=7),strip.background = element_blank(),legend.position="top")

ggsave(plot = bugbase_plot,filename = "bugbase_plot.pdf",height = 3,width = 5)
