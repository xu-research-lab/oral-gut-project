library(ggplot2)
library(dplyr)
load("../data/Extended Data Fig/Extended Data Fig5/Extended_Data_Fig5ab host_DNA_stat.RData")
df <- reads %>%
  mutate(subjectID = factor(subjectID, levels = unique(subjectID[order(site,treatment.1,-host_removed)])))

scientific_labels <- function(x) {
  parse(text = gsub("e\\+?", " %*% 10^", scales::scientific_format()(x)))
}


# 创建柱状图
saliva<-subset(df,site=="Saliva")
saliva <- saliva %>%
  mutate(subjectID = factor(subjectID, levels = unique(subjectID[order(site,treatment.1,-host_removed)])))


p_host_removed_saliva <- ggplot(saliva, aes(x = subjectID, y = host_removed, fill = subjectID)) +
  geom_col(stat = "identity", color = "black", linewidth = 0.01) +
  facet_grid(treatment.1 ~ site, scales = "fixed") +
  scale_y_log10() +  # 只设置对数变换，不设limits
  coord_cartesian(ylim = c(1e+6, 7.5e+07)) +  # 用这个来设置显示范围
  labs(title = "Host Removed Reads by Subject",
       x = "Subject ID",
       y = "Host Removed Reads Count") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

library(ggpubr)
p_host_rate<-ggplot(df,aes(site,rate*100,color=treatment.1))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2),size=0.5)+
  scale_y_log10()+
  stat_compare_means(label = "p.signif",show.legend = F,vjust = 0.3)+
  coord_cartesian(ylim = c(1e-5, 1.5e+2))+
  theme_bw()+
  labs(title = "propotion of host reads",x=NULL,y="propotion (%)",color="treatment")+
  theme(plot.title = element_text(hjust=0.5))

ggsave(plot = p_host_rate,filename = "boxplot_host_rate.pdf",height = 2.5,width = 3.5)

p_reads_stat<-ggplot(df,aes(site,host_removed,color=treatment.1))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2),size=0.5)+
  scale_y_log10()+
  stat_compare_means(label = "p.signif",show.legend = F,vjust = 0.3)+
  coord_cartesian(ylim = c(1e+6, 1e+8))+
  theme_bw()+
  labs(title = "clean data size",x=NULL,y="paired reads counts",color="treatment")+
  theme(plot.title = element_text(hjust=0.5))

ggsave(plot = p_reads_stat,filename = "../results/boxplot_data_reads.pdf",height = 2.5,width = 3.5)



