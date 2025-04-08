
setwd("/home/huangxc/software/SparCC3-master")
sgb<-read.table("SGB_profile.tsv",header = T,sep="\t",row.names = 1)
saliva<-sgb[,grepl("Saliva",colnames(sgb))]
feces<-sgb[,grepl("Feces",colnames(sgb))]
saliva_nonPMA<-saliva[,!grepl("PMA_|Control|reated",colnames(saliva))]
feces_nonPMA<-feces[,!grepl("PMA_|Control|reated",colnames(feces))]
saliva_PMA<-saliva[,grepl("PMA_",colnames(saliva))]
feces_PMA<-feces[,grepl("PMA_",colnames(feces))]
saliva_HC<-saliva[,grepl("_Control|Treated",colnames(saliva))]
feces_HC<-feces[,grepl("_Control|Treated",colnames(feces))]


saliva_PMA<-saliva_PMA[rowSums(saliva_PMA)>0,]
saliva_nonPMA<-saliva_nonPMA[rowSums(saliva_nonPMA)>0,]
feces_PMA<-feces_PMA[rowSums(feces_PMA)>0,]
feces_nonPMA<-feces_nonPMA[rowSums(feces_nonPMA)>0,]

feces_HC<-feces_HC[rowSums(feces_HC)>0,]
saliva_HC<-saliva_HC[rowSums(saliva_HC)>0,]

write.table(saliva_nonPMA*100000,file = "saliva_nonPMA.tsv",quote = F,row.names = T,sep="\t")
write.table(saliva_PMA*100000,file = "saliva_PMA.tsv",quote = F,row.names = T,sep="\t")
write.table(feces_nonPMA*100000,file = "feces_nonPMA.tsv",quote = F,row.names = T,sep="\t")
write.table(feces_PMA*100000,file = "feces_PMA.tsv",quote = F,row.names = T,sep="\t")
write.table(feces_HC*100000,file = "feces_HC.tsv",quote = F,row.names = T,sep="\t")
write.table(saliva_HC*100000,file = "saliva_HC.tsv",quote = F,row.names = T,sep="\t")


saliva_nonPMA<-melt(as.matrix(saliva_nonPMA))

saliva_nonPMA<-saliva_nonPMA[saliva_nonPMA$value>0,]
saliva_nonPMA_count<-aggregate(saliva_nonPMA$Var1,by=list(saliva_nonPMA$Var1),length)

saliva_nonPMA_count$Group.1<-gsub(".*s__","",saliva_nonPMA_count$Group.1)
saliva_nonPMA_count$Group.1<-gsub("t__","",saliva_nonPMA_count$Group.1)

colnames(saliva_nonPMA_count)<-c("sgb","counts")



saliva_PMA<-melt(as.matrix(saliva_PMA))

saliva_PMA<-saliva_PMA[saliva_PMA$value>0,]
saliva_PMA_count<-aggregate(saliva_PMA$Var1,by=list(saliva_PMA$Var1),length)

saliva_PMA_count$Group.1<-gsub(".*s__","",saliva_PMA_count$Group.1)
saliva_PMA_count$Group.1<-gsub("t__","",saliva_PMA_count$Group.1)

colnames(saliva_PMA_count)<-c("sgb","counts")

saliva_HC<-melt(as.matrix(saliva_HC))

saliva_HC<-saliva_HC[saliva_HC$value>0,]
saliva_HC_count<-aggregate(saliva_HC$Var1,by=list(saliva_HC$Var1),length)

saliva_HC_count$Group.1<-gsub(".*s__","",saliva_HC_count$Group.1)
saliva_HC_count$Group.1<-gsub("t__","",saliva_HC_count$Group.1)

colnames(saliva_HC_count)<-c("sgb","counts")


feces_nonPMA<-melt(as.matrix(feces_nonPMA))

feces_nonPMA<-feces_nonPMA[feces_nonPMA$value>0,]
feces_nonPMA_count<-aggregate(feces_nonPMA$Var1,by=list(feces_nonPMA$Var1),length)

feces_nonPMA_count$Group.1<-gsub(".*s__","",feces_nonPMA_count$Group.1)
feces_nonPMA_count$Group.1<-gsub("t__","",feces_nonPMA_count$Group.1)

colnames(feces_nonPMA_count)<-c("sgb","counts")



feces_PMA<-melt(as.matrix(feces_PMA))

feces_PMA<-feces_PMA[feces_PMA$value>0,]
feces_PMA_count<-aggregate(feces_PMA$Var1,by=list(feces_PMA$Var1),length)

feces_PMA_count$Group.1<-gsub(".*s__","",feces_PMA_count$Group.1)
feces_PMA_count$Group.1<-gsub("t__","",feces_PMA_count$Group.1)

colnames(feces_PMA_count)<-c("sgb","counts")

feces_HC<-melt(as.matrix(feces_HC))

feces_HC<-feces_HC[feces_HC$value>0,]
feces_HC_count<-aggregate(feces_HC$Var1,by=list(feces_HC$Var1),length)

feces_HC_count$Group.1<-gsub(".*s__","",feces_HC_count$Group.1)
feces_HC_count$Group.1<-gsub("t__","",feces_HC_count$Group.1)

colnames(feces_HC_count)<-c("sgb","counts")
