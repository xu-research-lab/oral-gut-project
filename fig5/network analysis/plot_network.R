setwd("/home/huangxc/software/SparCC3-master/feces_HC_pvals")
g.feces_nonPMA<-read.graph("network.gml",format = "gml")
shared_trans$SGB<-gsub("s__","",shared_trans$SGB)
temp<-V(g.feces_nonPMA)$name
temp<-data.frame(order=1:length(temp),sgb=temp)
temp$sgb<-gsub('[.]t__',"|",temp$sgb)
temp<-merge(temp,shared_trans[shared_trans$treatment=="HC",c("SGB","shared_rate","transmission.rate")],by.x = 2,by.y = 1,all.x = T)
numbers<-data.frame(sgb=unique(trans3$SGB),code=1:37)
temp<-merge(temp,numbers,by=1,all.x = T)
temp<-merge(temp,feces_HC_count,by=1,all.x = T)
temp$color[temp$shared_rate==0|is.na(temp$shared_rate)]<-"blue"
temp$color[temp$shared_rate>0]<-"orange"
temp$color[temp$`transmission.rate`>0]<-"red"
temp<-temp[order(temp$order),]

V(g.feces_nonPMA)$vertex.color<-temp$color
V(g.feces_nonPMA)$vertex.size<-temp$counts
V(g.feces_nonPMA)$vertex.label<-temp$code

g.feces_nonPMA<-delete.vertices(g.feces_nonPMA,which(degree(g.feces_nonPMA)==0))
#V(g.feces_nonPMA)$vertex.size[is.na(V(g.feces_nonPMA)$vertex.size)]<-1
pdf("feces_PMA.pdf",height = 6,width=6)
plot(g.feces_nonPMA,vertex.color=V(g.feces_nonPMA)$vertex.color,vertex.size=log2(V(g.feces_nonPMA)$vertex.size),vertex.label=V(g.feces_nonPMA)$vertex.label,vertex.label.cex=0.5,edge.lty=ifelse(E(g.feces_nonPMA)$sparcc<0,2,1),edge.width=E(g.feces_nonPMA)$weight*2)
dev.off()
node_attribute<-temp
save(g.feces_nonPMA,node_attribute,numbers,shared_trans,feces_nonPMA_count,file="feces_HC.RData")
