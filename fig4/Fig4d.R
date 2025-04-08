library(ggplot2)
library(dplyr)
library(tidyr)
p6 <- ggplot(subset(tmp4, id!="NA")) + 
  geom_point(aes(x = rate, y = id,color=Response_6),size = 3) +
  geom_line(aes(group = interaction(SGB, rate_type), x = rate, y = id), color = "gray") +
  geom_text(data=subset(tmp4, Response_6 == "R"),aes(x=-0.05,y=id,label = sigall,vjust=0.7), size = 5,color = "black") +
  geom_text(data=subset(tmp4, Response_6 == "NR"),aes(x=0.45,y=id,label = sigall,vjust=0.7), size = 5,color = "black") +
  scale_color_manual(name="Response",values = c("R" = "#377EB8", "NR" = "#E41A1C")) +
  #facet_wrap(~rate, nrow = 1) +
  facet_grid(rows = vars(group), cols = vars(rate_type))+
  theme_bw() +
  xlab(NULL)+ylab(NULL) +
  expand_limits(x = 0)+
  theme(#axis.text = element_text(size = 7,color = counts$color),
    axis.ticks.y = element_blank(),axis.text.y = element_text(size = 6))
p6
ggsave("trans_R_NR_shared_rate_new_segment.pdf",p6,height = 3,width = 6.5)

