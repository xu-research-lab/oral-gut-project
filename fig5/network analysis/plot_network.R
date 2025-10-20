library(ggraph)
library(ggplot2)
library(igraph)
load("../data/Figure5/network analysis/feces_nonPMA.RData")
pdf("../results/feces_nonPMA.pdf",height = 6,width=6)
plot(g.feces_nonPMA,vertex.color=V(g.feces_nonPMA)$vertex.color,
     vertex.size=log2(V(g.feces_nonPMA)$vertex.size),
     vertex.label=NA,
     edge.lty=ifelse(E(g.feces_nonPMA)$sparcc<0,2,1),
     edge.width=E(g.feces_nonPMA)$weight*2,
     vertex.frame.color = NA)
dev.off()


