p_nGD_response_nonpma<-ggplot(subset(trans3,!is.na(Response_6)&treatment.x=="non-PMA"),aes(Diagnosis.x,centered.min,text=Response_6))+
    geom_boxplot(aes(color=Response_6),outlier.shape = NA)+
    geom_point(aes(color=Response_6),position = position_jitterdodge(jitter.width = 0.2),show.legend = T)+
    stat_compare_means(method = "wilcox.test",label = "p.signif",hide.ns = T,vjust = 1)+
    ylim(c(0,2))+theme_bw()+scale_color_discrete(name="Response",type = c( "#377EB8","#E41A1C"))+
    geom_text(aes(label=NA),show.legend = F)+xlab(NULL)+ylab("nGD")
	
	ggsave(p_nGD_response_nonpma,"nGD_response_nonPMA_diff.pdf",height = 2.5,width = 3.5)
	