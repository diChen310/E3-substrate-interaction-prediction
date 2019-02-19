############Multiple classifiers based on different negative pairs###################
library(mlbench)
library(caret)
library(DMwR)
library(randomForest)
tumor<-"BRCA"
pairs.r<-read.table(file="E:/project/ubq prediction/R space/R space/pairs.r v-new.txt",sep = '\t')
pairs<-read.table(file="E:/project/ubq prediction/R space/R space/pairs v-new2.txt",sep = '\t')
pairs.in<-read.table(file="E:/project/ubq prediction/R space/R space/pairs.in v-new.txt",sep = '\t')
ppis.other.r<-read.table(file="E:/project/ubq prediction/R space/R space/other ppis v-new.txt",sep = '\t')
pairs<-pairs[pairs$e3.genes != "CDH1",]
all.features<-c("RCT","RCN","RCF","PCT","CNR.CXNR","CNR.CXNP","CCR.PPI","CCP.PPI","CNR.PPI","WCR","WCP","WCRP","WCPP","WCRS","WCRPS","WCPS","WCPPS")
names.set<-c("e3.genes","sub.genes")
#mrna.features<-c("RCT","RCN","RCF","CNR.CXNR","CCR.PPI","CNR.PPI","WCR","WCRP","WCRS","WCRPS")

features<-union(names.set,all.features)
pairs<-pairs[,features]
pairs.r<-pairs.r[,features]
ppis.other.r<-ppis.other.r[,features]
pairs.in<-pairs.in[,features]

pairs.full<-pairs[-manyNAs(pairs,0.001),]
pairs.r.full<-pairs.r[-manyNAs(pairs.r,0.001),]
ppis.full<-ppis.other.r[-manyNAs(ppis.other.r,0.001),]
pairs.in.full<-pairs.in[-manyNAs(pairs.in,0.001),]

features<-all.features

#########5-repeat of crossover validation##########
o<-1:nrow(pairs.full)
train_in_index<-sample(o,500)
m<-length(train_in_index)
rfFuncs$summary = twoClassSummary
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
pairs.tags<-c(rep("T",m),rep("N",m))

for(k in 1:5){
  pairs.r.sub<-pairs.r.full[(1+(k-1)*m):(k*m),]
  all1<-rbind(pairs.full[train_in_index,features],pairs.r.sub[,features])
  XXt<-as.matrix(all1)
  yf<-factor(pairs.tags)
  # # # # run the RFE algorithm
  results <- rfe(XXt, yf, sizes=c(1:ncol(all1)), rfeControl= control ,metric =  "ROC")
  save(results,file=paste0("crossover-randomVSE3U-",k,'.RData'))
  print(results)
 }

for(k in 1:5){
  pairs.r.sub<-ppis.full[(1+(k-1)*m):(k*m),]
  all2<-rbind(pairs.full[train_in_index,features],pairs.r.sub[,features])
  XXt<-as.matrix(all2)
  yf<-factor(pairs.tags)
  results <- rfe(XXt, yf, sizes=c(1:ncol(all2)), rfeControl= control ,metric =  "ROC")
  save(results,file=paste0("crossover-ppisVSE3U-",k,'.RData'))
  print(results)

}
for(k in 1:5){
  pairs.r.sub<-pairs.in.full[(1+(k-1)*m):(k*m),]
  all3<-rbind(pairs.full[train_in_index,features],pairs.r.sub[,features])
  XXt<-as.matrix(all3)
  yf<-factor(pairs.tags)
  # # # # # run the RFE algorithm
  results <- rfe(XXt, yf, sizes=c(1:ncol(all3)), rfeControl= control ,metric =  "ROC")
  print(results)
  save(results,file=paste0("crossover-IndirectVSE3U-",k,'.RData'))
  }



###########################Independent validation################################
k<-6

pairs.r.test<-rbind(pairs.full[o %in% train_in_index == F,],pairs.r.full[(1+(k-1)*m):(520+(k-1)*m),],pairs.in.full[(1+(k-1)*m):(520+(k-1)*m),],ppis.full[(1+(k-1)*m):(520+(k-1)*m),])
pairs.r.test$two<-paste(pairs.r.test$e3.genes,pairs.r.test$sub.genes)
pairs.r.test$from<-c(rep("ESI",nrow(pairs.full[o %in% train_in_index == F,])),rep("r",520),rep("in",520),rep('ppi',520))
pairs.r.test<-pairs.r.test[!duplicated(pairs.r.test$two),]
pairs.r.test.esi<-pairs.r.test[pairs.r.test$from == 'ESI',]
pairs.r.test.r<-pairs.r.test[pairs.r.test$from == 'r',]
pairs.r.test.in<-pairs.r.test[pairs.r.test$from == 'in',]
pairs.r.test.ppi<-pairs.r.test[pairs.r.test$from == 'ppi',]
pairs.r.test<-rbind(pairs.r.test.esi,pairs.r.test.r[sample(1:nrow(pairs.r.test.r),500),],pairs.r.test.in[sample(1:nrow(pairs.r.test.in),500),],pairs.r.test.ppi[sample(1:nrow(pairs.r.test.ppi),500),])


r.test.preds<-data.frame(row.names = pairs.r.test$two)
xx_r_test<-pairs.r.test[,features]
#accs1<-c(0.8847,0.8681,0.8831,0.8727,0.8750)
#accs2<-c(0.8723,0.8607,0.8602,0.8485,0.8490)
#accs3<-c(0.7921,0.7893,0.7999,0.7970,0.7933)

accs1<-c()
accs2<-c()
accs3<-c()


for(k in 1:5){
  load(file=paste0("crossover-randomVSE3U-",k,'.RData'))
  print(results)
  acc1<-max(results$results$ROC)
  accs1<-append(accs1,acc1)
  yf_pred<-predict(results,xx_r_test)
  r.test.preds[,k]<-yf_pred[,'T']
}

for(k in 1:5){
  load(file=paste0("crossover-ppisVSE3U-",k,'.RData'))
  print(results)
  acc2<-max(results$results$ROC)
  accs2<-append(accs2,acc2)

  yf_pred<-predict(results,xx_r_test)
  r.test.preds[,k+5]<-yf_pred[,'T']
}

for(k in 1:5){
  load(file=paste0("crossover-IndirectVSE3U-",k,'.RData'))
  print(results)
  acc3<-max(results$results$ROC)
  accs3<-append(accs3,acc3)

  yf_pred<-predict(results,xx_r_test)
  r.test.preds[,k+10]<-yf_pred[,'T']
}




r.test.preds$mPreds1<-apply(r.test.preds[,1:5],1,mean)
pairs.r.test$mPreds1<-r.test.preds$mPreds1
r.test.preds$mPreds2<-apply(r.test.preds[,6:10],1,mean)
pairs.r.test$mPreds2<-r.test.preds$mPreds2
r.test.preds$mPreds3<-apply(r.test.preds[,11:15],1,mean)
pairs.r.test$mPreds3<-r.test.preds$mPreds3
pairs.r.test$meanwP<-(mean(accs1)*pairs.r.test$mPreds1+mean(accs2)*pairs.r.test$mPreds2+mean(accs3)*pairs.r.test$mPreds3)/3

tags<-numeric(nrow(pairs.r.test))
for(i in 1:nrow(pairs.r.test)){
  tags[i] <- min(pairs.r.test[i,'mPreds1'],pairs.r.test[i,'mPreds2'],pairs.r.test[i,'mPreds3'])
}
pairs.r.test$tags<-tags>0.4

write.table(pairs.r.test,file='Independent test results.txt',sep='\t')
###################################FBXL predictions#############################################
fbxl.pairs<-read.delim(file='fbxl-predictions results-V new.txt',sep = '\t')
o<-1:nrow(pairs.full)
fbxl.genes<-c("SKP2", "FBXL2", "FBXL3", "FBXL4", "FBXL5", "FBXL6", "FBXL7", "FBXL8", "LRRC29", "KDM2B", "KDM2A", "FBXL12", 
              "FBXL13", "FBXL14", "FBXL15", "FBXL16", "FBXL17", "FBXL18", "FBXL19", "FBXL20", "FBXL21", "FBXL22")
pairs.r.full<-pairs.r.full[sample(1:nrow(pairs.r.full),nrow(pairs.r.full)),]
ppis.full<-ppis.full[sample(1:nrow(ppis.full),nrow(ppis.full)),]
pairs.in.full<-pairs.in.full[sample(1:nrow(pairs.in.full),nrow(pairs.in.full)),]
train_in_index<-o[which(as.character(pairs.full$e3.genes) %in% fbxl.genes == F)]####588###
m<-length(train_in_index)
rfFuncs$summary = twoClassSummary
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
pairs.tags<-c(rep("T",m),rep("N",m))


for(k in 1:5){
  pairs.r.sub<-pairs.r.full[(1+(k-1)*m):(k*m),]
  all1<-rbind(pairs.full[train_in_index,features],pairs.r.sub[,features])
  XXt<-as.matrix(all1)
  yf<-factor(pairs.tags)
  # # # # run the RFE algorithm
  results <- rfe(XXt, yf, sizes=c(1:ncol(all1)), rfeControl= control ,metric =  "ROC")
  save(results,file=paste0("fbxl-randomVSE3U-",k,'.RData'))
  print(results)
  xx_fbxl<-fbxl.pairs[,features]
  
  yf_pred<-predict(results,xx_fbxl)
  write.table(yf_pred,file=paste0("fbxl predicts-v1",k,".txt"),sep='\t')
  # 
}

for(k in 1:5){
  pairs.r.sub<-ppis.full[(1+(k-1)*m):(k*m),]
  all2<-rbind(pairs.full[train_in_index,features],pairs.r.sub[,features])
  XXt<-as.matrix(all2)
  yf<-factor(pairs.tags)
  results <- rfe(XXt, yf, sizes=c(1:ncol(all2)), rfeControl= control ,metric =  "ROC")
  save(results,file=paste0("fbxl-ppisVSE3U-",k,'.RData'))
  print(results)
  # summarize the results
  xx_fbxl<-fbxl.pairs[,features]
  #xx_fbxl[is.na(xx_fbxl)]<-0
  yf_pred<-predict(results,xx_fbxl)
  write.table(yf_pred,file=paste0("fbxl predicts-v2",k,".txt"),sep='\t')
  
}
for(k in 1:5){
  pairs.r.sub<-pairs.in.full[(1+(k-1)*m):(k*m),]
  all3<-rbind(pairs.full[train_in_index,features],pairs.r.sub[,features])
  XXt<-as.matrix(all3)
  yf<-factor(pairs.tags)
  # # # # # run the RFE algorithm
  results <- rfe(XXt, yf, sizes=c(1:ncol(all3)), rfeControl= control ,metric =  "ROC")
  print(results)
  save(results,file=paste0("fbxl-IndirectVSE3U-",k,'.RData'))
  xx_fbxl<-fbxl.pairs[,features]
  #xx_fbxl[is.na(xx_fbxl)]<-0
  yf_pred<-predict(results,xx_fbxl)
  write.table(yf_pred,file=paste0("fbxl predicts-v3",k,".txt"),sep='\t')
}
accs1<-c()
accs2<-c()
accs3<-c()

fbxl.pairs$two<-paste(fbxl.pairs$e3.genes,fbxl.pairs$sub.genes)
test.preds<-data.frame(row.names = fbxl.pairs$two)
for( k in 1:5){
  load(file=paste0("fbxl-randomVSE3U-",k,'.RData'))
  acc1<-max(results$results$ROC)
  accs1<-append(accs1,acc1)
  yf_pred<-read.table(file=paste0("fbxl predicts-v1",k,".txt"),sep='\t')
  
  test.preds[,k]<-yf_pred[,'T']
}
for( k in 1:5){
  load(file=paste0("fbxl-ppisVSE3U-",k,'.RData'))
  acc2<-max(results$results$ROC)
  accs2<-append(accs2,acc2)
  
  yf_pred<-read.table(file=paste0("fbxl predicts-v2",k,".txt"),sep='\t')
  test.preds[,k+5]<-yf_pred[,'T']
}
for( k in 1:5){
  load(file=paste0("fbxl-IndirectVSE3U-",k,'.RData'))
  acc3<-max(results$results$ROC)
  accs3<-append(accs3,acc3)
  yf_pred<-read.table(file=paste0("fbxl predicts-v3",k,".txt"),sep='\t')
  test.preds[,k+10]<-yf_pred[,'T']
}







fbxl.pairs$mPreds1<-apply(test.preds[,1:5],1,mean)
fbxl.pairs$mPreds2<-apply(test.preds[,6:10],1,mean)
fbxl.pairs$mPreds3<-apply(test.preds[,11:15],1,mean)
fbxl.pairs$meanwP<-(mean(accs1)*fbxl.pairs$mPreds1+mean(accs2)*fbxl.pairs$mPreds2+mean(accs3)*fbxl.pairs$mPreds3)/3
fbxl.pairs$meanP<-(fbxl.pairs$mPreds1+fbxl.pairs$mPreds2+fbxl.pairs$mPreds3)/3


write.table(fbxl.pairs,file='fbxl-predictions results-V new.txt',sep = '\t')

##############################All ESIs were collected into the training dataset#######################################
candidate.pairs<-read.delim(file=paste0("E:/project/ubq prediction/R space/R space/tumors/",tumor,"/candidate pairs not fbxl-v new.txt"),sep = '\t')
o<-1:nrow(pairs.full)
pairs.r.full<-pairs.r.full[sample(1:nrow(pairs.r.full),nrow(pairs.r.full)),]
ppis.full<-ppis.full[sample(1:nrow(ppis.full),nrow(ppis.full)),]
pairs.in.full<-pairs.in.full[sample(1:nrow(pairs.in.full),nrow(pairs.in.full)),]
train_in_index<-o
m<-length(train_in_index)
rfFuncs$summary = twoClassSummary
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
pairs.tags<-c(rep("T",m),rep("N",m))


for(k in 1:5){
  pairs.r.sub<-pairs.r.full[(1+(k-1)*m):(k*m),]
  all1<-rbind(pairs.full[train_in_index,features],pairs.r.sub[,features])
  XXt<-as.matrix(all1)
  yf<-factor(pairs.tags)
  # # # # run the RFE algorithm
  results <- rfe(XXt, yf, sizes=c(1:ncol(all1)), rfeControl= control ,metric =  "ROC")
  save(results,file=paste0("all-randomVSE3U-",k,'.RData'))
  print(results)
  xx<-candidate.pairs[,features]
  yf_pred<-predict(results,xx)
  write.table(yf_pred,file=paste0("candidate predicts-v1",k,".txt"),sep='\t')
  # 
}

for(k in 1:5){
  pairs.r.sub<-ppis.full[(1+(k-1)*m):(k*m),]
  all2<-rbind(pairs.full[train_in_index,features],pairs.r.sub[,features])
  XXt<-as.matrix(all2)
  yf<-factor(pairs.tags)
  results <- rfe(XXt, yf, sizes=c(1:ncol(all2)), rfeControl= control ,metric =  "ROC")
  save(results,file=paste0("all-ppisVSE3U-",k,'.RData'))
  print(results)
  # summarize the results
  xx<-candidate.pairs[,features]
  #xx_fbxl[is.na(xx_fbxl)]<-0
  yf_pred<-predict(results,xx)
  write.table(yf_pred,file=paste0("candidate predicts-v2",k,".txt"),sep='\t')
  
}
for(k in 1:5){
  pairs.r.sub<-pairs.in.full[(1+(k-1)*m):(k*m),]
  all3<-rbind(pairs.full[train_in_index,features],pairs.r.sub[,features])
  XXt<-as.matrix(all3)
  yf<-factor(pairs.tags)
  # # # # # run the RFE algorithm
  results <- rfe(XXt, yf, sizes=c(1:ncol(all3)), rfeControl= control ,metric =  "ROC")
  print(results)
  save(results,file=paste0("all-IndirectVSE3U-",k,'.RData'))
  xx<-candidate.pairs[,features]
  #xx_fbxl[is.na(xx_fbxl)]<-0
  yf_pred<-predict(results,xx)
  write.table(yf_pred,file=paste0("candidate predicts-v3",k,".txt"),sep='\t')
}
accs1<-c()
accs2<-c()
accs3<-c()

candidate.pairs$two<-paste(candidate.pairs$e3.genes,candidate.pairs$sub.genes)
test.preds<-data.frame(row.names = candidate.pairs$two)
for( k in 1:5){
  load(file=paste0("all-randomVSE3U-",k,'.RData'))
  acc1<-max(results$results$ROC)
  accs1<-append(accs1,acc1)
  yf_pred<-read.table(file=paste0("candidate predicts-v1",k,".txt"),sep='\t')
  
  test.preds[,k]<-yf_pred[,'T']
}
for( k in 1:5){
  load(file=paste0("all-ppisVSE3U-",k,'.RData'))
  acc2<-max(results$results$ROC)
  accs2<-append(accs2,acc2)
  
  yf_pred<-read.table(file=paste0("candidate predicts-v2",k,".txt"),sep='\t')
  test.preds[,k+5]<-yf_pred[,'T']
}
for( k in 1:5){
  load(file=paste0("fbxl-IndirectVSE3U-",k,'.RData'))
  acc3<-max(results$results$ROC)
  accs3<-append(accs3,acc3)
  yf_pred<-read.table(file=paste0("candidate predicts-v3",k,".txt"),sep='\t')
  test.preds[,k+10]<-yf_pred[,'T']
}







candidate.pairs$mPreds1<-apply(test.preds[,1:5],1,mean)
candidate.pairs$mPreds2<-apply(test.preds[,6:10],1,mean)
candidate.pairs$mPreds3<-apply(test.preds[,11:15],1,mean)
candidate.pairs$meanwP<-(mean(accs1)*candidate.pairs$mPreds1+mean(accs2)*candidate.pairs$mPreds2+mean(accs3)*candidate.pairs$mPreds3)/3
candidate.pairs$meanP<-(candidate.pairs$mPreds1+candidate.pairs$mPreds2+candidate.pairs$mPreds3)/3


write.table(candidate.pairs,file='candidate predictions results-V new.txt',sep = '\t')




##########Compared to known ESIs##############
fbxl.pairs<-read.delim(file='E:/project/ubq prediction/R space/R space/classification/final/fbxl-predictions results-V new.txt',sep = '\t')

#fbxl.pairs<-fbxl.pairs[fbxl.pairs$cate == 'Less',]
fbxl.pairs$two<-paste(fbxl.pairs$e3.genes,fbxl.pairs$sub.genes)
pairs.old<-pairs
pairs.old$two<-paste(pairs.old$e3.genes,pairs.old$sub.genes)
fbxl.pairs.confirmed<-read.delim(file="E:/project/ubq prediction/R space/R space/in pairs results/FBXL-UBQPAIRS-TEST.txt",sep='\t',col.names=c("e3.genes",'sub.genes','source','infomation'))
fbxl.pairs.confirmed$two=paste(fbxl.pairs.confirmed$e3.genes,fbxl.pairs.confirmed$sub.genes)
fbxl.pairs.test<-pairs.old[which(pairs.old$e3.genes %in% fbxl.genes),]
fbxl.pairs.test<-fbxl.pairs.test[which(fbxl.pairs.test$sub.genes != "RBX1"),]
fbxl.pairs.test<-fbxl.pairs[which(as.character(fbxl.pairs$two) %in% as.character(union(fbxl.pairs.test$two,fbxl.pairs.confirmed$two))),]
fbxl.pairs.test<-fbxl.pairs.test[as.character(fbxl.pairs.test$e3.genes) != as.character(fbxl.pairs.test$sub.genes),]

tags<-numeric(nrow(fbxl.pairs.test))
for(i in 1:nrow(fbxl.pairs.test)){
  tags[i] <- min(fbxl.pairs.test[i,'mPreds1'],fbxl.pairs.test[i,'mPreds2'],fbxl.pairs.test[i,'mPreds3'])
}
fbxl.pairs.test$tags<-tags>0.4

write.table(fbxl.pairs.test,file="fbxl test pairs-three classes results.txt",sep='\t')
fbxl.pairs.test<-read.delim(file="classification/final/fbxl test pairs-three classes results.txt",sep='\t')
fbxl.pairs.sub<-fbxl.pairs[fbxl.pairs$mPreds1>0.4,]
fbxl.pairs.sub<-fbxl.pairs.sub[fbxl.pairs.sub$mPreds2>0.4,]
fbxl.pairs.sub<-fbxl.pairs.sub[fbxl.pairs.sub$mPreds3>0.4,]

fbxl.pairs.sub<-fbxl.pairs.sub[which(fbxl.pairs.sub$meanwP > 0.65),]

fbxl.pairs.sub<-rbind(fbxl.pairs.sub,fbxl.pairs.test[,colnames(fbxl.pairs.sub)])
fbxl.pairs.sub<-fbxl.pairs.sub[!duplicated(fbxl.pairs.sub$two),]
fbxl.pairs.sub$type<-fbxl.pairs.sub$two %in% fbxl.pairs.test$two
write.table(fbxl.pairs.sub,file = 'fbxl candidate pairs and test pairs-repeated preds-three classes-0.65-all.txt',sep ='\t')

putative.pairs<-read.delim(file = "in pairs results/FBXL-APC pairs.txt",sep = '\t')
putative.pairs$two<-paste(putative.pairs$e3.genes,putative.pairs$sub.genes)
cosmic.genes<-read.delim(file= "E:/data/COSMIC/Census gene list.tsv",sep = '\t')
cosmic.genes.names<-as.character(unique(cosmic.genes[,1]))
fbxl.pairs.sub<-fbxl.pairs.sub[which(fbxl.pairs.sub$sub.genes %in% cosmic.genes.names),]
fbxl.pairs.old<-read.delim(file=paste0("E:/project/ubq prediction/R space/R space/tumors/",tumor,"/fbxl pairs-fill vTP.txt"),sep = '\t')
rownames(fbxl.pairs.old)<-paste(fbxl.pairs.old$e3.genes,fbxl.pairs.old$sub.genes)
output.fbxl.pairs<-rbind(fbxl.pairs.sub[,c('e3.genes','sub.genes','mPres.tm','two')],fbxl.pairs.test[,c('e3.genes','sub.genes','mPres.tm','two')])
output.fbxl.pairs<-output.fbxl.pairs[!duplicated(output.fbxl.pairs$two),]
output.fbxl.pairs$type<-output.fbxl.pairs$two %in% fbxl.pairs.test$two
output.fbxl.pairs$ltm<-output.fbxl.pairs$mPres.tm >=0.5
output.fbxl.pairs$cor.paths<-fbxl.pairs.old[output.fbxl.pairs$two,'WCR.Paths']
output.fbxl.pairs$proP.paths<-fbxl.pairs.old[output.fbxl.pairs$two,'proP.paths']
output.fbxl.pairs$WCRPath.S<-fbxl.pairs.old[output.fbxl.pairs$two,'WCRS.Paths']
output.fbxl.pairs$WCPPath.S<-fbxl.pairs.old[output.fbxl.pairs$two,'WCPPath.S']



write.table(fbxl.pairs.sub,file = "fbxl pairs results/fbxl candidate pairs-repeated preds-three classes-three0.4-cosmic.txt",sep = '\t', quote = F,row.names=F)
fbxl6<-fbxl.pairs[fbxl.pairs$e3.genes == "FBXL6",]
fbxl6<-fbxl6[fbxl6$mPres.tm >0.5,]
fbxl6<-fbxl6[fbxl6$mPreds1 >0.5,]
fbxl6<-fbxl6[fbxl6$mPreds2 >0.5,]
fbxl6<-fbxl6[fbxl6$mPreds3 >0.5,]
fbxl6$pathways<-fbxl.pairs.old[fbxl6$two,'cor.paths']
write.csv(fbxl6,file = "fbxl pairs results/fbxl6-0.5-two classifiers V6.csv", quote = F,row.names=F)

path.count<-count(output.fbxl.pairs,vars = 'cor.paths')
rownames(path.count)<-path.count$cor.paths
o<-order(path.count$freq,decreasing = T)
path.count<-path.count[o,]
path.count.l<-path.count[1:10,]
path.count.l[11,'cor.paths']<-"Others"
path.count.l[11,'freq']<-sum(path.count[11:nrow(path.count),'freq'])
pie(path.count.l$freq,path.count$cor.paths,lwd=2)

###############plot repeats#################

repeats<-read.csv(file= "E:/project/ubq prediction/R space/R space - fbxl prediction/two classes/prediction_five repeats.csv")

matplot(repeats,type='l',pch=1,col=1:5)
legend('topright',legend = colnames(repeats),col=1:5,pch=1)



fbxl.scors<-fbxl.pairs[,c('two','meanwP')]
fbxl.scors$type<-fbxl.pairs$two %in% fbxl.pairs.test$two

ggplot(fbxl.scors, aes(meanwP, fill = type)) + geom_histogram(alpha = 0.2, aes(y = ..density..), position = 'identity',bins = 20)


ggplot(data = fbxl.pairs.sub, aes(x=e3.genes, y=mPres.tm)) + geom_boxplot(col='blue')




