##Multidimensional association feature calculation###
library(TCGAbiolinks)
library(igraph)
library(CCA)
library(parallel)
##############Omics-based###########
expcorr<-function(e3.exp,sub.pro,type = "kendall"){
  
  ep.corr <-try(cor.test(e3.exp, sub.pro, method = type),silent = T) 
  if('try-error' %in% class(ep.corr)){
    ep.cor<-NA
    ep.p<-1.1
  }else{
    ep.cor <- ep.corr$estimate
    ep.p <- ep.corr$p.value
  }
  
  # }
  
  if(is.na(ep.cor)){
    ep.cor<-NA
    ep.p<-1.1
  }
  
  return(list(ep.cor,ep.p))
}

##Transcriptomics-based
ubqPairsCorrBasedOnExp <-function(exp.genes, ubq.pairs.in) {
  samples.genes <- colnames(exp.genes)
  samples.genes <- substr(samples.genes, 1, 15)
  colnames(exp.genes) <- samples.genes
  
  samplesTP <-
    TCGAquery_SampleTypes(barcode = colnames(exp.genes),
                          typesample = c("TP"))
  samplesNT <-
    TCGAquery_SampleTypes(barcode = colnames(exp.genes),
                          typesample = c("NT"))
  
  samplesNT.patients<-substr(samplesNT,1,12)
  samplesTP.patients<-substr(samplesTP,1,12)
  barcodesNT<-samplesNT[!duplicated(substr(samplesNT,1,12))]
  barcodesTP<-samplesTP[!duplicated(substr(samplesTP,1,12))]
  barcodesTP<-data.frame(barcodesTP)
  barcodesNT<-data.frame(barcodesNT)
  
  rownames(barcodesNT)<-substr(barcodesNT[,1],1,12)
  rownames(barcodesTP)<-substr(barcodesTP[,1],1,12)
  paired.patients<-intersect(rownames(barcodesNT),rownames(barcodesTP))
  #print(paste("The number of paired samples: ",length(paired.patients)))
  paired.samplesNT<-as.character(barcodesNT[paired.patients,])
  paired.samplesTP<-as.character(barcodesTP[paired.patients,])
 
  e.cors.nt <- c()
  
  e.cors.tp <- c()
  
  e.fold.cors <- c()
 
  pb <- txtProgressBar(min = 0, max = nrow(ubq.pairs.in), style = 3)
  skip = 0
  for (i in 1:nrow(ubq.pairs.in)) {
    e3.gene <- as.character(ubq.pairs.in[i, "e3.genes"])
    sub <- as.character(ubq.pairs.in[i, "sub.genes"])
  
    
    ###############TP and NT#################
    e3.exp.nt <- as.double(exp.genes[e3.gene, samplesNT])
    sub.exp.nt <- as.double(exp.genes[sub, samplesNT])
    e3.exp.nt <- log2(e3.exp.nt + 1e-10)
    sub.exp.nt <- log2(sub.exp.nt + 1e-10)
    e3.exp.tp <- as.double(exp.genes[e3.gene, samplesTP])
    sub.exp.tp <- as.double(exp.genes[sub, samplesTP])
    e3.exp.tp <- log2(e3.exp.tp + 1e-10)
    sub.exp.tp <- log2(sub.exp.tp + 1e-10)
    
    if (length(na.omit(sub.exp.nt)) <= 2 ||
        length(na.omit(e3.exp.nt)) <= 2) {
      e.cor.nt <- 0
      
    } else{
      e.corr <- cor.test(e3.exp.nt, sub.exp.nt)
      e.cor.nt <- e.corr$estimate
      
    }
    
    e.cors.nt <- append(e.cors.nt, e.cor.nt)
    #e.ps.nt <- append(e.ps.nt, e.p.nt)
    
    if (length(na.omit(sub.exp.tp)) <= 2 ||
        length(na.omit(e3.exp.tp)) <= 2) {
      e.cor.tp <- 0
      #e.p.tp <- 1.1
    } else{
      e.corr <- cor.test(e3.exp.tp, sub.exp.tp)
      e.cor.tp <- e.corr$estimate
      #e.p.tp <- e.corr$p.value
    }
    
    e.cors.tp <- append(e.cors.tp, e.cor.tp)
   
    
    #######Fold change corrs#############
    e3.exp.FC<-(exp.genes[e3.gene,paired.samplesTP]+0.0001)/(exp.genes[e3.gene,paired.samplesNT]+0.0001)
    sub.exp.FC<-(exp.genes[sub,paired.samplesTP]+0.0001)/(exp.genes[sub,paired.samplesNT]+0.0001)
    
    e3.exp.FC<-as.double(e3.exp.FC)
    sub.exp.FC<-as.double(sub.exp.FC)
    
    if (length(na.omit(sub.exp.FC)) <= 2 ||
        length(na.omit(e3.exp.FC)) <= 2) {
      e.fold.cor <- 0
      
    } else{
      res <- cor.test(e3.exp.FC, sub.exp.FC,method = 'kendall')
      e.fold.cor <- res$estimate
     
    }
    e.fold.cors <- append(e.fold.cors, e.fold.cor)
  
    setTxtProgressBar(pb, i)
    
  }
  
  ubq.pairs.in$RCN <- e.cors.nt
 
  ubq.pairs.in$RCT <- e.cors.tp
 
  ubq.pairs.in$RCF <- e.fold.cors

  
  return(ubq.pairs.in)
}

##Proteomics-based
ubqPairsCorrBasedOnCPTAC<-function(pairs,exp.pros){
  corrs.gene<-numeric(nrow(pairs))
  ps.gene<-numeric(nrow(pairs))
  for(i in 1:nrow(pairs)){
    
    e3<-as.character(pairs[i,'e3.genes'])
    sub<-as.character(pairs[i,'sub.genes'])
    print(paste(i,e3,sub))
    e3.exp<-as.double(exp.pros[e3,])
    sub.exp<-as.double(exp.pros[sub,])
    res<-expcorr(e3.exp,sub.exp)
    corrs.gene[i]<-res[[1]]
    ps.gene[i]<-res[[2]]
    
  }
  
  pairs$PCT<-corrs.gene

  return(pairs)
}

########Network-based###############

##Network
getPPINetWork<-function(){
  
  ppis<-read.delim(file="E:/data/PPI/BIOGRID.txt",sep='\t')#BioGrid PPIs
  ppis$two1<-paste(ppis$OFFICIAL_SYMBOL_A,ppis$OFFICIAL_SYMBOL_B)
  ppis$two2<-paste(ppis$OFFICIAL_SYMBOL_B,ppis$OFFICIAL_SYMBOL_A)
  
  ppis<-ppis[as.character(ppis$two1) != as.character(ppis$two2),]
  g<-graph.data.frame(ppis,directed = F)
  return(g)
}

PPIBasedSim<-function(g,ubp.pairs.in){
  
  nei.cNs<-c()
  diss<-c()
  for(i in 1:nrow(ubp.pairs.in)){
    e3<-as.character(ubp.pairs.in[i,'e3.genes'])
    sub<-as.character(ubp.pairs.in[i,'sub.genes'])
    if(e3 %in% V(g)$name && sub %in% V(g)$name){
      sp<-shortest_paths(g, from= e3, to = sub)$vpath[[1]]
      dis<-length(sp)-1
      
      e3.neis<-unique(neighbors(g,e3))
      sub.neis<-unique(neighbors(g,sub))
      cc.neis<-intersect(e3.neis,sub.neis)
      r<-length(cc.neis)/sqrt(length(e3.neis)*length(sub.neis))
    }else{
      
      dis<-999
      r<-NA
    }
    
    diss<-append(diss,dis)
    nei.cNs<-append(nei.cNs,r)
  }
  
  pairs$CNR.PPI<-nei.cNs
  pairs$dis<-diss
  return(pairs)
}
removeNAGenes<-function(exp.mat){
  genes.old<-rownames(exp.mat)
  
  genes.new<-c()
  row.has.na <- apply(exp.mat, 1, function(x){any(is.na(x))})
  genes.new<-genes.old[!row.has.na]
  
  return(genes.new)
}
corrsForNetwork<-function(exp.mat,genes){
  samplesTP <-
    TCGAquery_SampleTypes(barcode = colnames(exp.mat),
                          typesample = c("TP"))
  exp.mat<-exp.mat[genes,samplesTP]
  exp.mat<-t(as.matrix(exp.mat))
  corr.mac<-matcor(exp.mat,exp.mat)
  corr.mat<-corr.mac$XYcor
  return(corr.mat)
}


coExpNetwork<-function(corr.mat,genes,th=0.3){
  
  corr.mat[abs(corr.mat)<th]<-0
  corr.mat[abs(corr.mat)>=th]<-1
  network<-graph.adjacency(corr.mat,weighted=NULL)
  return(network)
}

commonNeighRates<-function(g,ubp.pairs.in){
  
  nei.cNs<-c()
  
  for(i in 1:nrow(ubp.pairs.in)){
    e3<-as.character(ubp.pairs.in[i,'e3.genes'])
    sub<-as.character(ubp.pairs.in[i,'sub.genes'])
    if(e3 %in% V(g)$name && sub %in% V(g)$name){
      e3.neis<-unique(neighbors(g,e3))
      sub.neis<-unique(neighbors(g,sub))
      cc.neis<-intersect(e3.neis,sub.neis)
      r<-length(cc.neis)/sqrt(length(e3.neis)*length(sub.neis))
    }else{
      r<-NA
    }
    nei.cNs<-append(nei.cNs,r)
  }
  return(nei.cNs)
}


mergeNeighbors<-function(g,nodes){
  neis<-c()
  for(node in nodes){
    neis.node<-names(unique(neighbors(g,node)))
    neis<-union(neis,neis.node)
  }
  return(neis)
}

commonCoexpressNeighRates<-function(g,corr.mat,ubp.pairs.in,Top=10){
  
  nei.cNs<-c()
  nei.SMs<-c()
  
  for(i in 1:nrow(ubp.pairs.in)){
    e3<-as.character(ubp.pairs.in[i,'e3.genes'])
    sub<-as.character(ubp.pairs.in[i,'sub.genes'])
    if( e3 %in% rownames(corr.mat) && sub %in% rownames(corr.mat)){
      e3.com<-corr.mat[e3,]
      e3.com.genes<-names(e3.com[order(as.double(e3.com),decreasing = T)[1:Top]])
      sub.com<-corr.mat[sub,]
      sub.com.genes<-names(sub.com[order(as.double(sub.com),decreasing = T)[1:Top]])
      
      e3.com.genes.in<-e3.com.genes[e3.com.genes %in% V(g)$name]
      sub.com.genes.in<-sub.com.genes[sub.com.genes %in% V(g)$name]
      e3.com.nei<-mergeNeighbors(g,e3.com.genes.in)
      sub.com.nei<-mergeNeighbors(g,sub.com.genes.in)
      e3.com.nei<-union(e3.com.genes,e3.com.nei)
      sub.com.nei<-union(sub.com.genes,sub.com.nei)
      cc.neis.top<-intersect(e3.com.nei,sub.com.nei)
      r<-length(cc.neis.top)/sqrt(length(e3.com.nei)*length(sub.com.nei))
    }else{
      r<-NA
    }
    
    nei.cNs<-append(nei.cNs,r)
    #print(i)
  }
  return(nei.cNs)
}

##############Pathway -based#############
pathwayExpRelationOnPermutaion<-function(ggc.m,e3.gene,genes.all,genePaths,pathIDToNames, repN=100){
  paths.ids<-unique(as.character(pathIDToNames[,"PathID"]))
  paths.names<-unique(as.character(pathIDToNames[,"PathName"]))
  paths.cors<-numeric(length(paths.ids))
  paths.cor.ps<-numeric(length(paths.ids))
  for(k in 1:length(paths.ids)){
    path.id<-paths.ids[k]
    path.genes<-unique(as.character(genePaths[which(genePaths$paths == path.id),"genes"]))
    path.genes<-intersect(path.genes,genes.all)
    path.genes<-path.genes[path.genes != e3.gene]#############only genes which are not the e3.gene itself#############
    path.genes.corrs<-ggc.m[e3.gene,path.genes]
    path.genes.corrs[is.na(path.genes.corrs)]<-0
    path.cor<-mean(abs(as.double(path.genes.corrs)))
    rep.rs<-numeric(repN)
    for(i in 1:repN){
      r.genes<-sample(genes.all,length(path.genes))
      r.genes.corrs<-ggc.m[e3.gene,r.genes]
      r.cor<-mean(abs(as.double(r.genes.corrs)))
      rep.rs[i]<-r.cor
    }
    prob <- length(rep.rs[rep.rs >= path.cor])/repN
    paths.cors[k]<-path.cor
    paths.cor.ps[k]<-prob
  }
  res<-data.frame(paths.ids,paths.names,paths.cors,paths.cor.ps)
  return(res)
}
e3subPathwayRelations<-function(ubp.pairs.in,e3.path.cors,genePaths,type = 'gene'){
  
  
  paths.cor.max<-numeric(nrow(ubp.pairs.in))
  paths.cor.ps<-numeric(nrow(ubp.pairs.in))
  cor.paths<-character(nrow(ubp.pairs.in))
  for(i in 1:nrow(ubp.pairs.in)){
    
    e3.gene<-as.character(ubp.pairs.in[i,'e3.genes'])
    sub.gene<-as.character(ubp.pairs.in[i,'sub.genes'])
    
    sub.paths.ids<-unique(as.character(genePaths[which(genePaths$genes == sub.gene),"paths"]))
    sub.paths<-as.character(pathIDToNames[which(pathIDToNames$PathID %in% sub.paths.ids),"PathName"])
    
    cor.max<-0
    cor.p<-1
    path.max<-NA
    cor.res<-e3.path.cors[which(e3.path.cors$e3 == e3.gene),]
    for(path in sub.paths){
      
      if( path %in% cor.res$paths ){
        
        cor_<- as.double(cor.res[which(cor.res$paths == path),"corrs"])
        p_<-as.double(cor.res[which(cor.res$paths == path),'ps'])
        if(is.na(cor_)){
          cor_<-0
          p_<-1.1
        }
        if(abs(cor_)>cor.max){
          cor.max<-cor_
          cor.p<-p_
          path.max<-path
        }
      }
      
    }
    if(is.na(path.max)){
      cor.max<-NA
      cor.p<-NA
    }
    paths.cor.max[i]<-cor.max
    paths.cor.ps[i]<-cor.p
    cor.paths[i]<-path.max
    #print(paste(i,e3.gene,sub.gene,cor.max,cor.p))
  }
  if( type == "gene"){
    ubp.pairs.in$WCR<-paths.cor.max
    ubp.pairs.in$WCRP<-1-paths.cor.ps
    ubp.pairs.in$WCR.paths<-cor.paths
  }else if(type == 'pro'){
    ubp.pairs.in$WCP<-paths.cor.max
    ubp.pairs.in$WCPP<-1-paths.cor.ps
    ubp.pairs.in$WCP.paths<-cor.paths
    
  }else if(type == 'gene.gene'){
    ubp.pairs.in$WCRG<-paths.cor.max
    ubp.pairs.in$WCRGP<-1-paths.cor.ps
    ubp.pairs.in$WCRG.paths<-cor.paths
  }else if(type == 'pro.gene'){
    ubp.pairs.in$WCPG<-paths.cor.max
    ubp.pairs.in$WCPGP<-1-paths.cor.ps
    ubp.pairs.in$WCPG.paths<-cor.paths
  }else{
    print("Input wrong type!!!!!!")
  }
  
  return(ubp.pairs.in)
}

newGenesForPathwayRelations<-function(tumor,genes,type = 'gene'){
  if(type == 'gene'){
    tumor.path<-paste0("tumors/",tumor,"/e3 paths TP")
    
  }else if(type == 'pro'){
    tumor.path<-paste0("tumors/",tumor,"/e3 pro paths")
    
  }else{
    print('Wrong type for input.')
  }
  recorded.files<-list.files(path=tumor.path)
  recorded.genes<-c()
  for(f in recorded.files){
    file.content<-strsplit(f," ")[[1]]
    if(file.content[1] == ""){
      gene<-file.content[2]
    }else{
      gene<-file.content[1]
    }
    recorded.genes<-append(recorded.genes,gene)
  }
  recorded.genes<-recorded.genes[!duplicated(recorded.genes)]
  
  genes.new<-genes[genes %in% recorded.genes == F]
  return(genes.new)
}
subE3PathwayRelations<-function(ubp.pairs.in,sub.path.cors,genePaths,type){
  #genePaths<-read.delim(file="E:\\data\\KEGG\\path_genes.txt",sep='\t')
  
  paths.cor.max<-numeric(nrow(ubp.pairs.in))
  paths.cor.ps<-numeric(nrow(ubp.pairs.in))
  cor.paths<-character(nrow(ubp.pairs.in))
  for(i in 1:nrow(ubp.pairs.in)){
    
    e3.gene<-as.character(ubp.pairs.in[i,'e3.genes'])
    sub.gene<-as.character(ubp.pairs.in[i,'sub.genes'])
    
    e3.paths.ids<-unique(as.character(genePaths[which(genePaths$genes == e3.gene),"paths"]))
    e3.paths<-as.character(pathIDToNames[which(pathIDToNames$PathID %in% e3.paths.ids),"PathName"])
    
    cor.max<-0
    cor.p<-1
    path.max<-NA
    cor.res<-sub.path.cors[which(sub.path.cors$e3 == sub.gene),]
    for(path in e3.paths){
      
      if( path %in% cor.res$paths ){
        
        cor_<- as.double(cor.res[which(cor.res$paths == path),"corrs"])
        p_<-as.double(cor.res[which(cor.res$paths == path),'ps'])
        
        if(abs(cor_)>cor.max){
          cor.max<-cor_
          cor.p<-p_
          path.max<-path
        }
      }
      
    }
    
    if(is.na(path.max)){
      cor.max<-NA
      cor.p<-NA
    }
    paths.cor.max[i]<-cor.max
    paths.cor.ps[i]<-cor.p
    cor.paths[i]<-path.max
  }
  
  
  
  if( type == "gene"){
    ubp.pairs.in$WCRS<-paths.cor.max
    ubp.pairs.in$WCRPS<-1-paths.cor.ps
    ubp.pairs.in$WCRS.Paths<-cor.paths
  }else if(type == 'pro'){
    ubp.pairs.in$WCPS<-paths.cor.max
    ubp.pairs.in$WCPPS<-1-paths.cor.ps
    ubp.pairs.in$WCPPS.Paths<-cor.paths
    
  }else if(type == 'gene.gene'){
    ubp.pairs.in$WCRG.S<-paths.cor.max
    ubp.pairs.in$WCRGP.S<-1-paths.cor.ps
    ubp.pairs.in$WCRGpaths.S<-cor.paths
  }else if(type == 'pro.gene'){
    ubp.pairs.in$WCPG.S<-paths.cor.max
    ubp.pairs.in$WCPGP.S<-1-paths.cor.ps
    ubp.pairs.in$WCPGpaths.S<-cor.paths
  }else{
    print("Input wrong type!!!!!!")
  }
  
  return(ubp.pairs.in)
}
###KEGG pathway infomation were downloaded and saved in local folder####
genePaths<-read.delim(file="E:\\data\\KEGG\\path_genes.txt",sep='\t')
genes.path<-unique(as.character(genePaths$genes))
pathIDToNames<-read.delim(file="E:\\data\\KEGG\\Path id to name.txt",sep='\t')

e3.paths<-data.frame(genes = e3.genes,paths = 'path:hsa04120')
genePaths<-rbind(genePaths,e3.paths)

setwd("E:/project/ubq prediction/R space/R space/")
tumor = "BRCA"
g<-getPPINetWork()
#######Begin computation###########
###RNA-seq data were downloaded from TCGA and saved in local folder: "E:\data\TCGA\FRDataFRData
exp.genes.t<-read.delim(file=paste0("E:\\data\\TCGA\\FRData\\",tumor,"ExpressionNormalized.txt"),check.names = F)

data.itraq<-read.delim(file = "E:/data/CPTAC/BRCA/TCGA_Breast_BI_Proteome_CDAP.r2.itraq.tsv" , sep = '\t',check.names = F)
rownames(data.itraq)<-data.itraq$Gene
data.itraq<-data.itraq[,-1]
data.itraq<-data.itraq[c(-1,-2,-3),]###Remove the mean and std values in the first three rows
desc.cols<-c("Description","Organism","Chromosome","Locus",'Descriptio')
data.itraq<-data.itraq[,colnames(data.itraq) %in% desc.cols == F]
data.itraq.lr.cols<-data.itraq[,!grepl("Unshared",colnames(data.itraq))]
colnames(data.itraq.lr.cols)<-substr(colnames(data.itraq.lr.cols),1,10)


pairs <- read.delim(file='pairs.txt')
###Change into random pairs in-direct pairs or FBXL candidate pairs for other types of input pairs
pairs<-ubqPairsCorrBasedOnExp(exp.genes.t,pairs)
pairs<-ubqPairsCorrBasedOnCPTAC(pairs,data.itraq.lr.cols)
#load(file=paste0("tumors/",tumor,"/corr.mat.mRNA.TP.RData"))
corr.mat.mRNA<-corrsForNetwork(exp.genes.t,intersect(genes.path,rownames(exp.genes.t)))
mRNA.g<-coExpNetwork(corr.mat.mRNA,rownames(corr.mat.mRNA),th=0.3)
corr.mat.mRNA<-corrsForNetwork(exp.genes.t,genes.all_)
pairs$CNR.CXNR<-commonNeighRates(mRNA.g,pairs)
pairs$CCR.PPI<-commonCoexpressNeighRates(g,corr.mat.mRNA,pairs)
pairs<-PPIBasedSim(g,pairs)
#load(file=paste0("tumors/",tumor,"/corr.mat.pro.RData"))
genes.g<-rownames(data.itraq.lr.cols)
#pro.g <- coExpNetwork(corr.mat.pro, genes.g, th = 0.3)


pairs$CNR.CXNP <- commonNeighRates(pro.g, pairs)
pairs$CCP.PPI<-commonCoexpressNeighRates(g,corr.mat.pro,pairs)

sub.genes<-unique(as.character(pairs$sub.genes))
e3.genes<-unique(as.character(pairs$e3.genes))

in.genes.all<-union(e3.genes,sub.genes)
in.genes.all<-intersect(in.genes.all,rownames(exp.genes.t))
genes.all_<-intersect(genes.path,rownames(exp.genes.t))
genes.all_<-intersect(genes.all_,rownames(corr.mat.mRNA))
genes.new<-newGenesForPathwayRelations(tumor,in.genes.all)
genes.new<-intersect(genes.new,rownames(corr.mat.mRNA))
library(parallel)
if(length(genes.new)>0){
  cl <- makeCluster(mc <- getOption("cl.cores", 6))
  
  #ggc.m<-cc(t(exp.genes.t[genes.new,]),t(exp.genes.t[genes.all_,]))$XYcor[genes.new,genes.all_]
  ggc.m<-corr.mat.mRNA[genes.new,genes.all_]
  clusterExport(cl=cl, varlist=c("pathwayExpRelationOnPermutaion", "ggc.m", "genes.all_", "genePaths",'pathIDToNames','tumor'))
  
  parLapply(cl,genes.new,function(gene,pathwayExpRelationOnPermutaion,ggc.m,genes.all_,genePaths,pathIDToNames,tumor){
    res<-pathwayExpRelationOnPermutaion(ggc.m,gene,genes.all_,genePaths,pathIDToNames)
    write.table(res,file=paste0("E:/project/ubq prediction/R space/R space/tumors/",tumor,"/e3 paths TP/",gene," related pathways.txt"),sep='\t')
  },pathwayExpRelationOnPermutaion,ggc.m,genes.all_,genePaths,pathIDToNames,tumor)
  
  stopCluster(cl)
  
  rm(ggc.m)
  gc()
}

cl <- makeCluster(mc <- getOption("cl.cores", 6))
clusterExport(cl=cl, varlist=c("tumor"))
results<-parLapply(cl,in.genes.all,function(gene,tumor){
  e3.path.corrs<-read.delim(file=paste0("E:/project/ubq prediction/R space/R space/tumors/",tumor,"/e3 paths TP/",gene," related pathways.txt"),sep = '\t')
  e3.path.corrs<-e3.path.corrs[,c(2,3,4)]
  e3.path.corrs$e3<-gene
  return(e3.path.corrs)
},tumor)

e3.path.cors <- do.call('rbind',results)
colnames(e3.path.cors)<-c('paths','corrs','ps','e3')
e3.path.cors<-e3.path.cors[,c('e3','paths','corrs','ps')]
stopCluster(cl)
pairs<-e3subPathwayRelations(pairs,e3.path.cors,genePaths,type = 'gene')
pairs<-subE3PathwayRelations(pairs,e3.path.cors,genePaths,type = 'gene')



genes.g<-rownames(data.itraq.lr.cols)

in.genes.all<-union(e3.genes,sub.genes)
in.genes.all<-intersect(in.genes.all,genes.g)
genes.new<-newGenesForPathwayRelations(tumor,in.genes.all,type='pro')
genes.all_<-intersect(genes.path,genes.g)
genes.all_<-intersect(genes.all_,rownames(corr.mat.pro))
if(length(genes.new)>0){
  genes.all_<-intersect(genes.path,rownames(data.itraq.lr.cols))
  
  ggc.m<-corr.mat.pro[genes.new,genes.all_]
  
  cl <- makeCluster(mc <- getOption("cl.cores", 6))
  clusterExport(cl=cl, varlist=c("pathwayExpRelationOnPermutaion", "ggc.m", "genes.all_", "genePaths",'pathIDToNames','tumor'))
  parLapply(cl,genes.new,function(gene,pathwayExpRelationOnPermutaion,ggc.m,genes.all_,genePaths,pathIDToNames,tumor){
    res<-pathwayExpRelationOnPermutaion(ggc.m,gene,genes.all_,genePaths,pathIDToNames)
    write.table(res,file=paste0("E:/project/ubq prediction/R space/R space/tumors/",tumor,"/e3 pro paths/",gene," related pathways.txt"),sep='\t')
  },pathwayExpRelationOnPermutaion,ggc.m,genes.all_,genePaths,pathIDToNames,tumor)
  stopCluster(cl)
  rm(ggc.m)
  gc()
}

cl <- makeCluster(mc <- getOption("cl.cores", 6))
clusterExport(cl=cl, varlist=c("tumor"))
in.genes.all<-intersect(rownames(data.itraq.lr.cols),in.genes.all)
results<-parLapply(cl,in.genes.all,function(gene,tumor){
  e3.path.corrs<-read.delim(file=paste0("tumors/",tumor,"/e3 pro paths/",gene," related pathways.txt"),sep = '\t')
  e3.path.corrs<-e3.path.corrs[,c(2,3,4)]
  e3.path.corrs$e3<-gene
  return(e3.path.corrs)
},tumor)

e3.path.cors <- do.call('rbind',results)
colnames(e3.path.cors)<-c('paths','corrs','ps','e3')
e3.path.cors<-e3.path.cors[,c('e3','paths','corrs','ps')]
stopCluster(cl)

pairs<-e3subPathwayRelations(pairs,e3.path.cors,genePaths,type = 'pro')
pairs<-subE3PathwayRelations(pairs,e3.path.cors,genePaths,type= 'pro')
write.table(pairs,file = "predictions-final.txt",sep = '\t')

