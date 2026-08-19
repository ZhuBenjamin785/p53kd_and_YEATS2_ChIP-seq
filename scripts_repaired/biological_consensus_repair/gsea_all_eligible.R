#!/usr/bin/env Rscript

# Preranked GSEA using every RNA gene with a finite Wald statistic.
suppressPackageStartupMessages({
  library(AnnotationDbi); library(clusterProfiler); library(org.Hs.eg.db)
  library(ReactomePA); library(dplyr)
})

argv <- commandArgs(trailingOnly=TRUE)
getarg <- function(flag, default=NULL) { i<-match(flag,argv); if(is.na(i)) return(default); argv[[i+1L]] }
input <- getarg("--results"); outdir <- getarg("--outdir")
if (is.null(input)||is.null(outdir)) stop("Required: --results FILE --outdir DIR")
min_size <- as.integer(getarg("--min-size","15")); max_size <- as.integer(getarg("--max-size","500"))
fdr <- as.numeric(getarg("--fdr","0.05")); seed <- as.integer(getarg("--seed","20260818"))
set.seed(seed); dir.create(file.path(outdir,"tables"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(outdir,"mapping"),recursive=TRUE,showWarnings=FALSE)

x <- read.csv(input,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(all(c("gene_id","stat","baseMean") %in% names(x)))
x$ENSEMBL <- sub("\\..*$","",x$gene_id); x$stat <- suppressWarnings(as.numeric(x$stat))
x$baseMean <- suppressWarnings(as.numeric(x$baseMean)); x <- x[is.finite(x$stat)&nzchar(x$ENSEMBL),]
mapping <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,keys=unique(x$ENSEMBL),keytype="ENSEMBL",
                                                  columns=c("ENSEMBL","ENTREZID","SYMBOL"))) |>
  filter(!is.na(ENTREZID)) |> distinct(ENSEMBL,ENTREZID,.keep_all=TRUE)
ranked <- inner_join(x,mapping,by="ENSEMBL") |>
  group_by(ENTREZID) |> summarise(stat=median(stat),SYMBOL=paste(sort(unique(na.omit(SYMBOL))),collapse=";"),
                                  Ensembl_rows=n(),.groups="drop") |>
  arrange(desc(stat),ENTREZID)
geneList <- ranked$stat; names(geneList) <- ranked$ENTREZID; geneList <- sort(geneList,decreasing=TRUE)
write.csv(ranked,file.path(outdir,"mapping","ranked_entrez_genes.csv"),row.names=FALSE)

empty <- function() data.frame(Pathway_ID=character(),Pathway=character(),Set_size=integer(),
  Enrichment_score=numeric(),NES=numeric(),P_value=numeric(),Adjusted_P_value_FDR=numeric(),
  Rank_at_max=integer(),Leading_edge=character(),Core_genes=character())
format_result <- function(obj,db) {
  if(is.null(obj)||!nrow(as.data.frame(obj))) return(empty())
  if(db!="GO_BP") obj <- setReadable(obj,OrgDb=org.Hs.eg.db,keyType="ENTREZID")
  as.data.frame(obj) |> transmute(Pathway_ID=ID,Pathway=Description,Set_size=setSize,
    Enrichment_score=enrichmentScore,NES=NES,P_value=pvalue,Adjusted_P_value_FDR=p.adjust,
    Rank_at_max=rank,Leading_edge=leading_edge,Core_genes=gsub("/",";",core_enrichment,fixed=TRUE)) |>
    arrange(Adjusted_P_value_FDR,desc(abs(NES)))
}
run_db <- function(db) tryCatch({
  obj <- switch(db,
    GO_BP=gseGO(geneList=geneList,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont="BP",exponent=1,
      minGSSize=min_size,maxGSSize=max_size,eps=0,pvalueCutoff=1,pAdjustMethod="BH",verbose=FALSE,seed=TRUE,by="fgsea"),
    KEGG=gseKEGG(geneList=geneList,organism="hsa",keyType="ncbi-geneid",exponent=1,
      minGSSize=min_size,maxGSSize=max_size,eps=0,pvalueCutoff=1,pAdjustMethod="BH",verbose=FALSE,seed=TRUE,by="fgsea"),
    Reactome=gsePathway(geneList=geneList,organism="human",exponent=1,minGSSize=min_size,
      maxGSSize=max_size,eps=0,pvalueCutoff=1,pAdjustMethod="BH",verbose=FALSE,seed=TRUE))
  list(data=format_result(obj,db),status="completed")
},error=function(e) list(data=empty(),status=paste("ERROR",conditionMessage(e))))

summary <- list()
for(db in c("GO_BP","KEGG","Reactome")) {
  ans<-run_db(db); out<-ans$data |> mutate(Database=db,.before=1)
  write.csv(out,file.path(outdir,"tables",paste0(db,"__all_pathways.csv")),row.names=FALSE)
  summary[[db]]<-data.frame(Database=db,Input_rows=nrow(x),Ranked_Entrez_genes=length(geneList),
    Tested_pathways=nrow(out),Significant_pathways=sum(out$Adjusted_P_value_FDR<fdr,na.rm=TRUE),
    Positive_NES_significant=sum(out$Adjusted_P_value_FDR<fdr & out$NES>0,na.rm=TRUE),
    Negative_NES_significant=sum(out$Adjusted_P_value_FDR<fdr & out$NES<0,na.rm=TRUE),
    FDR_threshold=fdr,Min_gene_set_size=min_size,Max_gene_set_size=max_size,Seed=seed,Status=ans$status)
}
write.csv(bind_rows(summary),file.path(outdir,"tables","gsea_run_summary.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
cat("All-gene GSEA complete\n"); print(bind_rows(summary))

