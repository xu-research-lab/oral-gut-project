library(ggplot2)
library(dplyr)
library(ggpubr)

load("alpha_beta_meta.RData")

p1 <- ggplot(subset(alpha_beta_meta,Diagnosis %in% c("CRC", "GC", "EC","RA","HC")&treatment=="non-PMA"),aes(observed.oral,observed.gut,color=Diagnosis))+
  geom_point(shape=19)+
  stat_cor(show.legend = F,size = 2.5)+
  stat_smooth(method = "lm",se = F,linewidth=0.5)+
  theme_bw()+
  #theme(axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 6))+
  xlab(label="Observed SGBs (oral) ")+
  ylab(label="Observed SGBs (gut) ")+
  #facet_grid(~treatment)+
  scale_color_manual(values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3","RA"="#8da0cb", "HC"="#a6d854"))#+
  theme(strip.background = element_blank(),strip.text = element_text(size = 11))
p1
ggsave("Figure2A_change_color_non-PMA.pdf",plot=p1,height = 2.6,width = 3.75)

