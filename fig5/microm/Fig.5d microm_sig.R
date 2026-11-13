
library(ggplot2)
load("../data/Figure5/microm/microm_sig.RData")
## file microm - association between transmitted abundance and all predicted metabolites
## sig_metabolites - significant associations
#microm<-read.delim("clipboard",header = T)
sig_metabolites <- microm[abs(microm$statistic) >= 0.2&microm$n>=5&microm$q<=0.25, ]

# Check whether any rows meet the condition
if(nrow(sig_metabolites) > 0) {
  # Draw bar plot
  ggplot(sig_metabolites, aes(x = reorder(name, statistic), y = statistic)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +  # Flip coordinates to make long names easier to read
    labs(x = "Metabolite", y = "Spearman r") +
    theme_bw()
} else {
  message("No metabolites meet the condition (all |statistic| < 2)")
}
ggsave("../results/microm_sig.pdf",height = 2,width = 2.5)
