#!/usr/bin/env Rscript

# Direction-specific RNA ORA using the actual RNA testing universe.
suppressPackageStartupMessages({
  library(AnnotationDbi); library(clusterProfiler); library(org.Hs.eg.db)
  library(ReactomePA); library(dplyr)
})

argv <- commandArgs(trailingOnly=TRUE)
getarg <- function(flag, default=NULL) {
  i <- match(flag, argv); if (is.na(i)) return(default)
  if (i == length(argv)) stop("Missing value for ", flag); argv[[i+1L]]
}
background_file <- getarg("--background")
up_file <- getarg("--up")
down_file <- getarg("--down")
outdir <- getarg("--outdir")
if (any(vapply(list(background_file,up_file,down_file,outdir), is.null, logical(1))))
  stop("Required: --background FILE --up FILE --down FILE --outdir DIR")
fdr <- as.numeric(getarg("--fdr", "0.05")); min_size <- as.integer(getarg("--min-size", "10"))
max_size <- as.integer(getarg("--max-size", "500"))
dir.create(file.path(outdir,"tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir,"mapping"), recursive=TRUE, showWarnings=FALSE)

read_genes <- function(path) {
  x <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  hit <- names(x)[tolower(names(x)) %in% c("gene","symbol","gene_name","gene_symbol")]
  if (!length(hit)) stop("No gene-symbol column in ", path)
  unique(toupper(trimws(x[[hit[[1L]]]])))[nzchar(unique(toupper(trimws(x[[hit[[1L]]]]))))]
}
map_symbols <- function(x, label) {
  if (!length(x)) return(data.frame(Set=character(),SYMBOL=character(),ENTREZID=character()))
  suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=x, keytype="SYMBOL",
                                         columns=c("SYMBOL","ENTREZID"))) |>
    filter(!is.na(ENTREZID)) |> distinct(SYMBOL,ENTREZID) |> mutate(Set=label,.before=1)
}
symbols <- list(background=read_genes(background_file), RNA_up=read_genes(up_file), RNA_down=read_genes(down_file))
maps <- Map(map_symbols, symbols, names(symbols))
universe <- unique(maps$background$ENTREZID)
genes <- lapply(maps[c("RNA_up","RNA_down")], function(x) intersect(unique(x$ENTREZID),universe))
write.csv(bind_rows(maps), file.path(outdir,"mapping","symbol_to_entrez.csv"), row.names=FALSE)

empty <- function() data.frame(Pathway_ID=character(),Pathway=character(),GeneRatio=character(),
  BgRatio=character(),Count=integer(),P_value=numeric(),Adjusted_P_value_FDR=numeric(),Genes=character())
format_result <- function(x, db) {
  if (is.null(x) || !nrow(as.data.frame(x))) return(empty())
  if (db != "GO_BP") x <- setReadable(x, OrgDb=org.Hs.eg.db, keyType="ENTREZID")
  as.data.frame(x) |> transmute(Pathway_ID=ID,Pathway=Description,GeneRatio=GeneRatio,
    BgRatio=BgRatio,Count=Count,P_value=pvalue,Adjusted_P_value_FDR=p.adjust,
    Genes=gsub("/",";",geneID,fixed=TRUE)) |> arrange(Adjusted_P_value_FDR,P_value)
}
run_db <- function(g, db) {
  if (!length(g)) return(list(data=empty(),status="empty foreground"))
  tryCatch({
    x <- switch(db,
      GO_BP=enrichGO(gene=g,universe=universe,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont="BP",
        pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1,minGSSize=min_size,maxGSSize=max_size,readable=TRUE),
      KEGG=enrichKEGG(gene=g,universe=universe,organism="hsa",keyType="ncbi-geneid",
        pAdjustMethod="BH",pvalueCutoff=1,qvalueCutoff=1,minGSSize=min_size,maxGSSize=max_size,use_internal_data=FALSE),
      Reactome=enrichPathway(gene=g,universe=universe,organism="human",pAdjustMethod="BH",
        pvalueCutoff=1,qvalueCutoff=1,minGSSize=min_size,maxGSSize=max_size,readable=FALSE))
    list(data=format_result(x,db),status="completed")
  }, error=function(e) list(data=empty(),status=paste("ERROR",conditionMessage(e))))
}

summary <- list()
for (direction in names(genes)) for (db in c("GO_BP","KEGG","Reactome")) {
  ans <- run_db(genes[[direction]],db)
  out <- ans$data |> mutate(Direction=direction,Database=db,.before=1)
  write.csv(out,file.path(outdir,"tables",paste0(direction,"__",db,"__all_pathways.csv")),row.names=FALSE)
  summary[[paste(direction,db)]] <- data.frame(Direction=direction,Database=db,
    Foreground_symbols=length(symbols[[direction]]),Foreground_Entrez_in_universe=length(genes[[direction]]),
    Background_symbols=length(symbols$background),Background_Entrez=length(universe),
    Tested_pathways=nrow(out),Significant_pathways=sum(out$Adjusted_P_value_FDR < fdr,na.rm=TRUE),
    FDR_threshold=fdr,Status=ans$status)
}
write.csv(bind_rows(summary),file.path(outdir,"tables","fea_run_summary.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
cat("Directional FEA complete\n"); print(bind_rows(summary))

