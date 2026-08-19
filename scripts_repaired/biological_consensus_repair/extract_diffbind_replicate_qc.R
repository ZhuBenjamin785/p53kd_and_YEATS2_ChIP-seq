#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DiffBind))
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2) stop("Usage: Rscript extract_diffbind_replicate_qc.R DBA.rds OUT.csv")
db<-readRDS(args[[1]]); gr<-dba.peakset(db,bRetrieve=TRUE)
x<-as.data.frame(mcols(gr)); x[]<-lapply(x,as.numeric)
logx<-log2(as.matrix(x)+1)
p<-cor(logx,use="pairwise.complete.obs",method="pearson")
s<-cor(logx,use="pairwise.complete.obs",method="spearman")
pairs<-t(combn(colnames(logx),2))
out<-data.frame(sample_1=pairs[,1],sample_2=pairs[,2],
 pearson_log2_score=apply(pairs,1,function(z)p[z[1],z[2]]),
 spearman_score=apply(pairs,1,function(z)s[z[1],z[2]]),tested_regions=nrow(logx))
write.csv(out,args[[2]],row.names=FALSE)
print(out)
