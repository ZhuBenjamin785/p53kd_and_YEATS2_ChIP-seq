#!/usr/bin/env Rscript

# Provisional sensitivity analysis on existing BAMs.  It fixes broad-peak
# handling and duplicate policy but cannot fix the upstream separate-species
# alignment; outputs must remain labelled provisional.
suppressPackageStartupMessages({
  library(DiffBind); library(ChIPseeker); library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db); library(GenomicRanges); library(BiocParallel)
})
args <- commandArgs(trailingOnly=TRUE)
outdir <- if(length(args)) args[[1]] else stop("Output directory required")
remove_duplicates <- if(length(args)>=2) as.logical(args[[2]]) else TRUE
root <- "/gpfs/projects/b1042/LauberthLab/BenFolder"
human_bam_dir <- if(length(args)>=3) args[[3]] else file.path(root,"p53kdH4K16ac","p53kdbamfiles","human")
dm6_bam_dir <- if(length(args)>=4) args[[4]] else file.path(root,"p53kdH4K16ac","p53kdbamfiles","dm6")
peak_root <- if(length(args)>=5) args[[5]] else file.path(root,"p53kdH4K16ac","macs3_results_p53kd","peaks")
result_status <- Sys.getenv("CHIP_RESULT_STATUS",
  unset="PROVISIONAL — upstream hg38/dm6 species disambiguation is unresolved")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
threads <- max(1L,as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","1")))
register(MulticoreParam(workers=threads,progressbar=FALSE),default=TRUE)
p <- function(...) file.path(root,...)
samples <- data.frame(
 SampleID=c("Scramble_rep1","Scramble_rep2","p53KD_rep1","p53KD_rep2"),
 Condition=c("Scramble","Scramble","p53KD","p53KD"),Replicate=c(1,2,1,2),
 bamReads=file.path(human_bam_dir,c("Scr_H4K16ac_1_S0_L001.sorted.bam","Scr_H4K16ac_2_S0_L001.sorted.bam","P53_H4K16ac_1_S0_L001.sorted.bam","P53_H4K16ac_2_S0_L001.sorted.bam")),
 bamControl=file.path(human_bam_dir,c("Scr_Input1_S0_L001.sorted.bam","Scr_Input2_S0_L001.sorted.bam","MutP53_Input1_S0_L001.sorted.bam","MutP53_Input2_S0_L001.sorted.bam")),
 Spikein=file.path(dm6_bam_dir,c("Scr_H4K16ac_1_S0_L001.sorted.bam","Scr_H4K16ac_2_S0_L001.sorted.bam","P53_H4K16ac_1_S0_L001.sorted.bam","P53_H4K16ac_2_S0_L001.sorted.bam")),
 Peaks=file.path(peak_root,c("Scramble_H4K16ac_rep1/Scramble_H4K16ac_rep1_peaks.broadPeak","Scramble_H4K16ac_rep2/Scramble_H4K16ac_rep2_peaks.broadPeak","p53KD_H4K16ac_rep1/p53KD_H4K16ac_rep1_peaks.broadPeak","p53KD_H4K16ac_rep2/p53KD_H4K16ac_rep2_peaks.broadPeak")),
 PeakCaller="macs",stringsAsFactors=FALSE)
stopifnot(all(file.exists(unlist(samples[c("bamReads","bamControl","Spikein","Peaks")]))))
write.csv(samples,file.path(outdir,"sample_sheet.csv"),row.names=FALSE)
db <- dba(sampleSheet=samples); db$config$cores<-threads; db$config$RunParallel<-threads>1
db <- dba.blacklist(db,blacklist=DBA_BLACKLIST_HG38,greylist=FALSE,cores=threads)
db$config$doGreylist<-FALSE; db$config$doBlacklist<-FALSE
db <- dba.count(db,summits=FALSE,bRemoveDuplicates=remove_duplicates,bParallel=threads>1,bSubControl=FALSE)
pdf(file.path(outdir,"PCA.pdf")); dba.plotPCA(db,attributes=DBA_CONDITION); dev.off()
pdf(file.path(outdir,"correlation_heatmap.pdf")); plot(db); dev.off()
db <- dba.normalize(db,normalize=DBA_NORM_LIB,spikein=TRUE)
db <- dba.contrast(db,categories=DBA_CONDITION,minMembers=2)
db <- dba.analyze(db,method=DBA_DESEQ2,bParallel=threads>1,bBlacklist=FALSE,bGreylist=FALSE)
all <- as.data.frame(dba.report(db,method=DBA_DESEQ2,th=1))
write.csv(all,file.path(outdir,"DiffBind_all_peaks.csv"),row.names=FALSE)
write.csv(all[!is.na(all$FDR)&all$FDR<0.05,],file.path(outdir,"DiffBind_significant_peaks.csv"),row.names=FALSE)
gr <- makeGRangesFromDataFrame(all,keep.extra.columns=TRUE)
anno <- as.data.frame(annotatePeak(gr,TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                                   tssRegion=c(-3000,3000),annoDb="org.Hs.eg.db"))
write.csv(anno,file.path(outdir,"diffbind_peak_annotations.csv"),row.names=FALSE)
saveRDS(db,file.path(outdir,"DiffBind_object.rds"))
capture.output(sessionInfo(),file=file.path(outdir,"sessionInfo.txt"))
writeLines(c(paste("STATUS:",result_status),
 paste("summits=FALSE; bRemoveDuplicates=",remove_duplicates),"blacklist=hg38; greylist=FALSE",
 paste("human_bam_dir=",human_bam_dir),paste("dm6_bam_dir=",dm6_bam_dir),paste("peak_root=",peak_root),
 "normalization=DBA_NORM_LIB with spikein=TRUE","contrast=p53KD - Scramble; DESeq2; FDR<0.05"),
 file.path(outdir,"PARAMETERS.txt"))
