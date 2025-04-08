library(optparse)

args <- OptionParser(usage = "%sheet name", add_help_option = TRUE, prog='randomforest.R')
args_list <- parse_args(args, positional_arguments=TRUE)

library(randomForest)
library(pROC)
library(caret)
library(ggplot2)
library(readxl)

all<-read.csv(args_list$args[2],header=T,row.names=1)
#all<-read.csv("oral_gut.csv",header=T,row.names=1)
meta<-read.table("metadata.txt",header=T,row.names=1,sep="\t")

data<-merge(meta[,c("Gender","Treatment","MSI_Status","Response_3","Response_6","Operation","Tumor_Staging","Diagnosis")],all,by=0,all=T)
rownames(data)<-data[,1]
data<-data[,-1]
data<-data[!grepl("_B",rownames(data)),]
#data<-data[!grepl("_C",rownames(data)),]
data<-subset(data,!grepl("PMA",rownames(data)))
#data<-subset(data,grepl("PMA",rownames(data))|grepl("_Control",rownames(data)))
data[is.na(data)]<-0

randomfor <- function(dat = data, test) {
  work_table <- data.frame(dat[, test], dat[, c(9:ncol(dat) - 1)])
  work_table[, 1] <- as.factor(work_table[, 1])
  
  # creat 5 fold cross-validate test data 
  n <- createFolds(work_table[, 1], k = 5, returnTrain = FALSE)
  obs_p_ran <- data.frame()
  
  for(i in 1:length(n)) {
    test_randomforest <- randomForest(work_table[-n[[i]], 1] ~ ., 
                                      data = work_table[-n[[i]], -1], 
                                      ntree = 1000, importance = TRUE)
    pre_ran_i <- predict(test_randomforest, newdata = work_table[n[[i]], ], type = "prob")
    obs_p_ran_i <- data.frame(id = rownames(work_table)[n[[i]]], 
                              prob = pre_ran_i[, 1])
    obs_p_ran <- rbind(obs_p_ran, obs_p_ran_i)
  }
  
  # merge results
  all_mean <- data.frame(id = rownames(work_table),
                         obs = work_table[, 1],
                         NR = 1 - obs_p_ran$prob[match(rownames(work_table), obs_p_ran$id)],
                         R = obs_p_ran$prob[match(rownames(work_table), obs_p_ran$id)])
  
  ran_roc <- roc(all_mean$obs, all_mean$R)
  return(ran_roc)
}

# Create new diagnosis groups for the combined comparisons
data$Diagnosis_Group <- data$Diagnosis

# CRC vs GC+EC
print("CRC vs GC+EC")
data1 <- data[which(data$Diagnosis %in% c("CRC", "GC", "EC")),]
data1$Diagnosis_Group <- ifelse(data1$Diagnosis == "CRC", "CRC", "GC_EC")
CRC_vs_GCEC <- randomfor(dat=data1, test="Diagnosis_Group")

# GC vs CRC+EC
print("GC vs CRC+EC")
data2 <- data[which(data$Diagnosis %in% c("CRC", "GC", "EC")),]
data2$Diagnosis_Group <- ifelse(data2$Diagnosis == "GC", "GC", "CRC_EC")
GC_vs_CRCEC <- randomfor(dat=data2, test="Diagnosis_Group")

# EC vs CRC+GC
print("EC vs CRC+GC")
data3 <- data[which(data$Diagnosis %in% c("CRC", "GC", "EC")),]
data3$Diagnosis_Group <- ifelse(data3$Diagnosis == "EC", "EC", "CRC_GC")
EC_vs_CRCGC <- randomfor(dat=data3, test="Diagnosis_Group")

print("first_step")
data4<-data[which(data$Diagnosis=="CRC"|data$Diagnosis=="GC"),]
CRC_GC<-randomfor(dat=data4,test="Diagnosis")

print("second_step")
data5<-data[which(data$Diagnosis=="CRC"|data$Diagnosis=="EC"),]
CRC_EC<-randomfor(dat=data5,test="Diagnosis")

data6<-data[which(data$Diagnosis=="GC"|data$Diagnosis=="EC"),]
EC_GC<-randomfor(dat=data6,test="Diagnosis")


## Each cancer vs HC
#print("CRC vs HC")
#data4 <- data[which(data$Diagnosis %in% c("CRC", "HC")),]
#CRC_vs_HC <- randomfor(dat=data4, test="Diagnosis")

#print("GC vs HC")
#data5 <- data[which(data$Diagnosis %in% c("GC", "HC")),]
#GC_vs_HC <- randomfor(dat=data5, test="Diagnosis")

#print("RC vs HC")
#data6 <- data[which(data$Diagnosis %in% c("RA", "HC")),]
#EC_vs_HC <- randomfor(dat=data6, test="Diagnosis")

# All cancers (GC+EC+CRC) vs HC
#print("All cancers vs HC")
#data7 <- data[which(data$Diagnosis %in% c("CRC", "GC", "EC", "HC")),]
#data7$Diagnosis_Group <- ifelse(data7$Diagnosis == "HC", "HC", "All_Cancers")
#AllCancers_vs_HC <- randomfor(dat=data7, test="Diagnosis_Group")

# Collect all AUC values
auc <- c(
  CRC_vs_GCEC = CRC_vs_GCEC$auc,
  GC_vs_CRCEC = GC_vs_CRCEC$auc,
  EC_vs_CRCGC = EC_vs_CRCGC$auc,
  CRC_GC=CRC_GC$auc,
  EC_GC=EC_GC$auc,
  CRC_EC=CRC_EC$auc
 # RA_vs_HC = EC_vs_HC$auc,
 # AllCancers_vs_HC = AllCancers_vs_HC$auc
)

# Print AUC values
print("OK")
write.table(auc,paste0(args_list$args[1],"_Diagnosis_AUC_nonPMA.txt"),quote=F,sep="\t")
