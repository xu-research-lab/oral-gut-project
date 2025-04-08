p1<-ggplot(subset(all,species=="Streptococcus salivarius"), aes(both_YES+gut_YES, both_YES+oral_YES)) + 
  geom_point(data = subset(subset(all,species=="Streptococcus salivarius"), !is.na(GO)), aes(color = GO), shape = 16) +
  geom_point(data = subset(subset(all,species=="Streptococcus salivarius"), is.na(GO)), color = "grey50", shape = 1,size=0.5) +
  facet_wrap(~species, scales = "free") +
  theme_bw()+xlab(label="gene present in gut")+
  ylab(label="gene present in oral")+
  geom_abline(slope=1,intercept = (-2.160913+3*8.197742),linetype="dashed")+
  geom_abline(slope =1,intercept = (-2.160913-3*8.197742),linetype="dashed")+
  theme(legend.position = "none")+geom_abline(slope=1)

p2<-ggplot(subset(all,species=="Streptococcus parasanguinis"), aes(both_YES+gut_YES, both_YES+oral_YES)) + 
  geom_point(data = subset(subset(all,species=="Streptococcus parasanguinis"), !is.na(GO)), aes(color = GO), shape = 16) +
  geom_point(data = subset(subset(all,species=="Streptococcus parasanguinis"), is.na(GO)), color = "grey50", shape = 1,size=0.5) +
  facet_wrap(~species, scales = "free") +
  theme_bw()+xlab(label="gene present in gut")+
geom_abline(slope=1,intercept = (-1.478430-3*4.430210),linetype="dashed")+
  geom_abline(slope =1,intercept = -1.478430+3*4.430210,linetype="dashed")+geom_abline(slope=1)+
  ylab(label=NULL)+xlim(c(0,60))+theme(legend.key.size = unit(0.1,"cm"))

Fig3.D<-p1%>%insert_right(p2)

ggsave("go richement.pdf",height =3 ,width = 7)