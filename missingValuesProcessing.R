library(DMwR)

features<-c("e3.genes","sub.genes","RCT","RCN","RCF","PCT",
            "CNR.CXNR","CNR.CXNP","CCR.PPI","CCP.PPI","CNR.PPI",
            "WCR","WCP","WCRP","WCPP","WCRS","WCRPS","WCPS","WCPPS")   



for(tumor in tumors){
  pairs.r<-read.delim(file=paste0("E:/project/ubq prediction/R space/R space/tumors/",tumor,"/random pairs v2-TP.txt"),sep = '\t')
  pairs<-read.delim(file=paste0("E:/project/ubq prediction/R space/R space/tumors/",tumor,"/pairs v2-TP.txt"),sep = '\t')
  ppis.other.r<-read.delim(file=paste0("E:/project/ubq prediction/R space/R space/tumors/",tumor,"/random ppis v2-TP.txt"),sep = '\t')
  pairs.in<-read.delim(file="E:/project/ubq prediction/R space/R space/tumors/BRCA/random indirect pairs v2-TP.txt",sep='\t')
  fbxl.pairs<-read.delim(file=paste0("E:/project/ubq prediction/R space/R space/tumors/",tumor,"/fbxl pairs-v2-TPWithIn.txt"),sep = '\t')
  
  pairs<-pairs[,features]
  pairs.r<-pairs.r[,features]
  ppis.other.r<-ppis.other.r[,features]
  pairs.in<-pairs.in[,features]
  fbxl.pairs<-fbxl.pairs[,features]
  
  pairs.full<-pairs[-manyNAs(pairs,0.001),]
  pairs.r.full<-pairs.r[-manyNAs(pairs.r,0.001),]
  ppis.full<-ppis.other.r[-manyNAs(ppis.other.r,0.001),]
  pairs.in.full<-pairs.in[-manyNAs(pairs.in,0.001),]
  
  #fbxl.pairs.na<-fbxl.pairs[manyNAs(fbxl.pairs,0.2),]
  fbxl.pairs<-fbxl.pairs[-manyNAs(fbxl.pairs,0.3),]
  all.pairs.full<-rbind(pairs.full[,features],pairs.r.full[,features],ppis.full[,features])
  save(all.pairs.full,file="D:/all.pairs.full.RData")
  
  
  #fbxl.pairs.na[,features] <- knnImputation(fbxl.pairs.na[,features], k = 10, distData = all.pairs.full[,features])
  fbxl.pairs[,features]<-knnImputation(fbxl.pairs[,features],k = 10,distData = all.pairs.full[,features])
  fbxl.pairs$cate<-"Less"
  fbxl.pairs.na$cate<-"More"
  fbxl.pairs<-rbind(fbxl.pairs,fbxl.pairs.na)
  
  #write.table(fbxl.pairs.na,file=paste0("E:/project/ubq prediction/R space/R space/tumors/",tumor,"/fbxl pairs na-fill.txt"),sep = '\t')
  write.table(fbxl.pairs,file=paste0("E:/project/ubq prediction/R space/R space/tumors/",tumor,"/fbxl pairs-fill v1-TPWithIn.txt"),sep = '\t')
  
}
