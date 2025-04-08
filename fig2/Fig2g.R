library(ggplot2)
library(dplyr)
library(aplot)
p_shared1<-ggplot(subset(df_sig2,rate>0&!grepl("Candidatus",id)&group=="non-PMA"),aes(Diagnosis,id))+
  geom_point(color ="#BC3C29",aes(size=rate))+
  geom_text(aes(label=sig),size=3,nudge_y = -0.25)+
  theme_bw()+
  theme(axis.text = element_text(size = 8),axis.text.y = element_text(color=counts1$color))+
  xlab(NULL)+ylab(NULL)+scale_size_continuous(name="shared rate",range = c(0.1,3),limits = c(0, 1),breaks = seq(0.2, 0.8, by = 0.3))+ggtitle("Shared rate")+theme(plot.title = element_text(size = 8,hjust=0.5))
p_shared1
p_cor1<-ggplot(subset(stat,!(Diagnosis=="Treated")&id>=6&group=="non-PMA"),aes(Diagnosis,paste(X1,X2),fill=cor))+geom_tile()+geom_text(aes(label=sig),size=3,vjust=0.8)+
  theme_bw()+
  theme(axis.text = element_text(size = 7),axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.grid = element_blank(),strip.background = element_blank(),strip.text.y.right = element_text(angle = 0,hjust=0),panel.spacing = unit(0.1,"mm"))+
  ylab(NULL)+xlab(NULL)+scale_x_discrete(expand = c(0,0))+
  scale_fill_gradient2(name="spearman R",low="#377EB8",mid = 'white', high = '#E41A1C',na.value = "white")+ggtitle("Oral-gut correlection")+theme(plot.title = element_text(size = 8,hjust=0.5))
p_cor1
p_all <- p_shared1%>%insert_right(p_cor1)
p_all
ggsave("p_shared_new20250819.pdf",p_all,height = 7,width = 7)

