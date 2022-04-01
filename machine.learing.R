###################### ARGS input
args=commandArgs(T)
outdir=args[1]   ## out dir
trainvar=args[2] ## training data variations ;row.name=sample; col.names=genes name;
traingro=args[3] ## training data group; case=1; control=0; only one column= 1 1 1 0 0 0 ...
predvar=args[4]  ## predict data variations ;row.name=sample; col.names=genes name;
predgro=args[5]  ## predict data group; case=1; control=0; only one column= 1 1 1 0 0 0 ...
casename=args[6] ## cases group name, such as "PTB"
controlname=args[7] ## control group name, such as "FTB"
#scaleYes=args[8]  ## whether scale, if Yes,value=1; if no, value=0;
#topCut=args[9]   ## default=30

#######################
library(ggplot2)
library(pROC)
library(ROCR)
library(ggsci)
require(caret)
library(dplyr)
library(glmnet)
library(randomForest)
library(gbm)
library(kernlab)
library(pheatmap)
source("train2pred.R")

######################### files input #########################################
gen1<-read.table(trainvar,header = T,sep="\t",row.names = 1) #
gen_raw=gen1  ##for heatmap
gen1<-log2(gen1+1)
gen2=scale(gen1, center = TRUE, scale = TRUE);
gen1<-as.data.frame(gen2)
varLen<-length(names(gen1))
varAll<-names(gen1)
gro<-read.table(traingro,sep="\t")
colnames(gro)<-c("gr")
gen1$PE<-gro$gr
gen1$PE<-factor(gen1$PE,levels=c(0,1),labels=c("No","Yes"))
gen_raw$PE<-gen1$PE

kk<-read.table(predvar,header = T,sep="\t",row.names = 1)
kk_raw<-kk
kk<-log2(kk+1)
kk <- as.data.frame(scale(kk,center = attr(gen2, "scaled:center"), scale = attr(gen2, "scaled:scale")))
grop<-read.table(predgro,sep="\t")
colnames(grop)<-c("gr")
kk$PE<-grop$gr
kk$PE<-factor(kk$PE,levels=c(0,1),labels=c("No","Yes"))
kk_raw$PE<-kk$PE


##################### features selection based on the importance of LR and RF ROC #################
features<-c() 
seeds<-c(1:20)*100
tryn<-length(seeds)
for (seed1 in seeds)
{
  set.seed(seed1) #2
  gen=gen1
  trControl <- trainControl(method = 'repeatedcv', number = 7,repeats = 1,classProbs = TRUE, summaryFunction = twoClassSummary,savePredictions = T)
###### glmnet for top 100 selection
  myParamGrid=expand.grid(lambda=c(0,0.001, 0.5, 0.01, 0.03, 0.05, 0.1,1),alpha=c(0.099,1,0.1))
  M1 <- train(PE~.,data=gen,method='glmnet',trControl=trControl,metric = "ROC",tuneGrid=myParamGrid)
  
  for (top1 in c(100,50,20,10,5))
  {
    if(varLen<5){var1=varAll;break;}  ## PPROM+GDM < 30
    if(varLen<top1){next;}
    importance = varImp(M1,scale = FALSE)
    imsort<-sort(importance$importance$Overall,decreasing=T) 
    cutoff<-imsort[top1]
    cutid<-which(importance$importance$Overall>=cutoff)
    var1<-(row.names(importance$importance))[cutid]
    var<-c(var1,"PE")
    gen<- gen[,names(gen) %in% var]
    
    myParamGrid <- expand.grid(mtry=c(1,2,3,4)) #c(1,2,3,4)
    M1 <- train(PE~.,data=gen,method="rf",trControl=trControl,metric = "ROC",tuneGrid = myParamGrid)
  }
  features<-c(features,var1)
}

featuresFrq<-as.data.frame(table(features))
var2<-as.character(featuresFrq[which(featuresFrq$Freq>=tryn*0.5),1])
#var2<-unique(c("CLCN3","DAPP1","MOB1B","RGS18","RAB27B","PPBP","MAP3K7CL"))
var<-c(var2,"PE")
gen<- gen1[,names(gen1) %in% var]
gen_raw1<-gen_raw[,names(gen_raw) %in% var]
write.table(var2,file=paste(outdir,"/features.txt",sep=""),quote = F)


##########################################print importance for features
imp<-function(M1,model,outdir)
{
	importance = varImp(M1,scale = TRUE)
	impScoreOrder<-importance$importance[order(-importance$importance[,1]),]
	impGneOrder<-row.names(importance$importance)[order(-importance$importance[,1])]
	impOder<-data.frame(impGneOrder,impScoreOrder)
	write.table(impOder,file=paste(outdir,"/",model,".features.imp.txt",sep=""),quote = F,row.names = F)
}


############################################### training ########################################
set.seed(1024)
trControl <- trainControl(method = 'repeatedcv', number = 7,repeats = 1,classProbs = TRUE, summaryFunction = twoClassSummary,savePredictions = T)

############################ training:RF
myParamGrid <- expand.grid(mtry=c(1:4)) #c(1,2,3,4)
M1 <- train(PE~.,data=gen,method="rf",trControl=trControl,metric = "ROC",tuneGrid = myParamGrid)
M1
M1$pred=M1$pred[which(M1$pred$mtry==M1$bestTune$mtry),]
imp(M1,"RF",outdir)
train2pred("RF",M1,gen,kk,outdir)

############################## training:glmnet
myParamGrid=expand.grid(lambda=c(0,0.001,0.5, 0.01, 0.03, 0.05, 0.1,1),alpha=c(0,0.000001,0.00001,0.0001,0.001,0.01,0.9,1,0.1,0.2,0.3,0.5,0.7,0.8))   #
M1 <- train(PE~.,data=gen,method='glmnet',trControl=trControl,metric = "ROC",tuneGrid=myParamGrid)
imp(M1,"LR",outdir)
lambda_list=M1$bestTune$lambda
alpha_list=M1$bestTune$alpha
M1$pred=M1$pred[which(M1$pred$lambda==lambda_list & M1$pred$alpha== alpha_list),]
train2pred("LR",M1,gen,kk,outdir)


###############################training:SVM
myParamGrid <- expand.grid(C=c(0.1,1,0.3,0.01,0.03,0.004,0.001,0.0001,0.002,0.5,0.8)) #0.0001-10000,greater the value,greater punishment for wrong cases,may result in Model overfitting
M1 <- train(PE~.,data=gen,method="svmLinear",trControl=trControl,tuneGrid=myParamGrid,metric = "ROC")
imp(M1,"SVM",outdir)
M1$pred=M1$pred[which(M1$pred$C==M1$bestTune$C),]
train2pred("SVM",M1,gen,kk,outdir)


############################ training:gbm
myParamGrid <- expand.grid(n.trees=c(100,300,500,800,1000),shrinkage=c(0.01,0.005,0.001,0.0001),interaction.depth=c(1,2,3,4),n.minobsinnode = c(10,20)) M1 <- train(PE~.,data=gen,method="gbm",trControl=trControl,metric = "ROC",tuneGrid = myParamGrid)
M1$pred=M1$pred[which(M1$pred$n.trees==M1$bestTune$n.trees & M1$pred$shrinkage==M1$bestTune$shrinkage & M1$pred$interaction.depth==M1$bestTune$interaction.depth & M1$pred$n.minobsinnode==M1$bestTune$n.minobsinnode),]
imp(M1,"GBM",outdir)
train2pred("GBM",M1,gen,kk,outdir)
