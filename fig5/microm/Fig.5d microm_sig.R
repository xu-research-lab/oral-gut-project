
load("microm_sig.RData")
## file microm - association between transmitted abundance and all predicted metabolites
## sig_metabolites - significant associations
microm<-read.delim("clipboard",header = T)
sig_metabolites <- microm[abs(microm$statistic) >= 0.2&microm$n>=5&microm$q<=0.25, ]

# 检查是否有满足条件的行
if(nrow(sig_metabolites) > 0) {
  # 绘制柱状图
  ggplot(sig_metabolites, aes(x = reorder(name, statistic), y = statistic)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +  # 使长名称更易阅读
    labs(x = "Metabolite", y = "Spearman r") +
    theme_bw()
} else {
  message("没有满足条件的metabolites（所有statistic绝对值都小于2）")
}
ggsave("microm_sig.pdf",height = 2,width = 2.5)
save(microm,sig_metabolites,file="microm_sig.RData")
