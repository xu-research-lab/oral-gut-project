library(ggplot2)
library(dplyr)
p5 <- ggplot(subset(remergeshare_new9, group=="non-PMA"&id!="NA")) + 
  geom_segment(aes(x = `trans. rate`, xend = rate, y = id, yend = id),color="grey70" ) +
  geom_point(aes(x = `trans. rate`, y = id,color=oral.x,shape="trans.rate"),size = 3) +
  geom_point(aes(x = rate, y = id,shape="shared.rate"),size = 2.5,color="grey70") +
  geom_text(aes(x=-0.1,y=id,label = sigtrans,vjust=0.78), size = 5) +
  geom_vline(xintercept = 0, color = "red", linetype = "solid", size = 0.5) +
  scale_shape_manual(values = c("trans.rate" = 16, "shared.rate" = 17), 
                     name = "Type") +
  facet_wrap(~Diagnosis, nrow = 1) +
  theme_bw() +
  xlab(NULL)+ylab(NULL) +
  expand_limits(x = 0)+
  theme(axis.text = element_text(size = 7,color = counts$color),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 9),  # 调整 x 轴刻度标签尺寸
        axis.text.y = element_text(size = 10))+
  scale_color_gradient(name="Present(oral)")
p5
ggsave("trans_non-PMA_shared_rate_new_segment_20250821.pdf",p5,height = 8,width = 12.5)

