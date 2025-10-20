library(ggplot2)
library(dplyr)
library(tidyr)
load("../data/Figure4/fig4d.RData")
p_data <- Species_shared_trans_diff %>% 
  pivot_longer(
    cols = shared_rate:tran_rate,     
    names_to = "rate_type",
    values_to = "rate"
  )
p_data<-p_data%>%mutate(sig=ifelse(adjusted_pvalue_shared<=0.25|adjusted_pvalue_trans<=0.25,"*",NA),sig=ifelse(is.na(adjusted_pvalue_trans)&rate_type=="tran_rate",NA,sig))
ggplot(subset(p_data, group=="non-PMA")) + 
  geom_point(aes(x = rate, y = id,color=Response_6),size = 3) +
  geom_line(aes(group = interaction(SGB, rate_type), x = rate, y = id), color = "gray") +
  geom_text(data=subset(p_data, rate_type == "shared_rate"),aes(x=0.45,y=id,label = sig,vjust=0.7), size = 5,color = "black") +
  geom_text(data=subset(p_data, rate_type == "tran_rate"),aes(x=0.45,y=id,label = sig,vjust=0.7), size = 5,color = "black") +
  scale_color_manual(name="Response",values = c("R" = "#377EB8", "NR" = "#E41A1C")) +
  #facet_wrap(~rate, nrow = 1) +
  facet_grid(cols = vars(rate_type))+
  theme_bw() +
  xlab(NULL)+ylab(NULL) +
  expand_limits(x = 0)+
  theme(#axis.text = element_text(size = 7,color = counts$color),
    axis.ticks.y = element_blank(),axis.text.y = element_text(size = 6))

ggsave("../results/trans_R_NR_shared_rate_new_segment.pdf",height = 3,width = 6.5)
