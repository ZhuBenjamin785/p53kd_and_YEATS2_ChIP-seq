library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(BiocParallel)

results_dir <- "mof_macs3_results/diffbind_results"
dir.create(results_dir, showWarnings=FALSE, recursive=TRUE)


n_threads <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
register(MulticoreParam(workers=n_threads, progressbar=FALSE), default=TRUE)



samples <- data.frame(
  SampleID = c("Scr_MOF_rep1", "Scr_MOF_rep2", "p53sh_MOF_rep1", "p53sh_MOF_rep2"),
  Condition = c("Scr_MOF", "Scr_MOF", "p53sh_MOF", "p53sh_MOF"),
  Replicate = c(1, 2, 1, 2),
  bamReads = c(
    "MOF_bamfiles/human/Scr_MOF_1.sorted.bam",
    "MOF_bamfiles/human/Scr_MOF_2.sorted.bam",
    "MOF_bamfiles/human/p53sh_MOF_1.sorted.bam",
    "MOF_bamfiles/human/p53sh_MOF_2.sorted.bam"
  ),
  bamControl = c(
    "MOF_bamfiles/human/Scr_Input1_S0_L001.sorted.bam",
    "MOF_bamfiles/human/Scr_Input2_S0_L001.sorted.bam",
    "MOF_bamfiles/human/MutP53_Input1_S0_L001.sorted.bam",
    "MOF_bamfiles/human/MutP53_Input2_S0_L001.sorted.bam"
  ),
  Spikein = c(
    "MOF_bamfiles/dm6/Scr_MOF_1.sorted.bam",
    "MOF_bamfiles/dm6/Scr_MOF_2.sorted.bam",
    "MOF_bamfiles/dm6/p53sh_MOF_1.sorted.bam",
    "MOF_bamfiles/dm6/p53sh_MOF_2.sorted.bam"
  ),
  Peaks = c(
    "mof_macs3_results/peaks/Scr_MOF_rep1/Scr_MOF_rep1_peaks.broadPeak",
    "mof_macs3_results/peaks/Scr_MOF_rep2/Scr_MOF_rep2_peaks.broadPeak",
    "mof_macs3_results/peaks/p53sh_MOF_rep1/p53sh_MOF_rep1_peaks.broadPeak",
    "mof_macs3_results/peaks/p53sh_MOF_rep2/p53sh_MOF_rep2_peaks.broadPeak"
  ),
  PeakCaller = "macs",
  stringsAsFactors = FALSE
)

stopifnot(all(file.exists(samples$bamReads)))
stopifnot(all(file.exists(samples$bamControl)))
stopifnot(all(file.exists(samples$Spikein)))
stopifnot(all(file.exists(samples$Peaks)))


write.csv(samples, file.path(results_dir, "diffbind_mof_samples.csv"), row.names = FALSE)
db <- dba(sampleSheet = samples)
db$config$cores <- n_threads
db$config$RunParallel <- n_threads > 1L

db <- dba.blacklist(
    db,
    blacklist=DBA_BLACKLIST_HG38,
    greylist=FALSE,
    cores=n_threads
)

db$config$doGreylist <- FALSE
db$config$doBlacklist <- FALSE

db <- dba.count(
    db,
    summits=250,
    bParallel=n_threads > 1L,
    bSubControl=FALSE
)

pdf(file.path(results_dir, "DiffBind_MOF_PCA.pdf"))
dba.plotPCA(db, attributes=DBA_CONDITION)
dev.off()

pdf(file.path(results_dir, "DiffBind_MOF_correlation_heatmap.pdf"))
plot(db)
dev.off()


db <- dba.normalize(
    db,
    normalize=DBA_NORM_LIB,
    spikein=TRUE
)

db <- dba.contrast(
    db,
    categories=DBA_CONDITION,
    minMembers=2
)

db <- dba.analyze(
    db,
    method=DBA_DESEQ2,
    bParallel=n_threads > 1L,
    bBlacklist=FALSE,
    bGreylist=FALSE
)

all_db <- dba.report(
    db,
    method=DBA_DESEQ2,
    th=1
)

all_df <- as.data.frame(all_db)

write.csv(
    all_df,
    file.path(results_dir, "DiffBind_all_peaks.csv"),
    row.names=FALSE
)



sig <- all_df[
    !is.na(all_df$FDR) &
    all_df$FDR < 0.05,
]

write.csv(
    sig,
    file.path(results_dir, "DiffBind_significant_peaks.csv"),
    row.names=FALSE
)

gained <- sig[sig$Fold > 0, ]
lost   <- sig[sig$Fold < 0, ]

write.csv(gained, file.path(results_dir, "MOF_gained_peaks.csv"), row.names=FALSE)
write.csv(lost,   file.path(results_dir, "MOF_lost_peaks.csv"), row.names=FALSE)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene



gain_gr <- makeGRangesFromDataFrame(
    gained,
    seqnames.field="seqnames",
    start.field="start",
    end.field="end",
    keep.extra.columns=TRUE
)

loss_gr <- makeGRangesFromDataFrame(
    lost,
    seqnames.field="seqnames",
    start.field="start",
    end.field="end",
    keep.extra.columns=TRUE
)

gain_anno <- annotatePeak(
    gain_gr,
    TxDb=txdb,
    tssRegion=c(-3000, 3000),
    annoDb="org.Hs.eg.db"
)

loss_anno <- annotatePeak(
    loss_gr,
    TxDb=txdb,
    tssRegion=c(-3000, 3000),
    annoDb="org.Hs.eg.db"
)

write.csv(
    as.data.frame(gain_anno),
    file.path(results_dir, "MOF_gained_annotated.csv"),
    row.names=FALSE
)

write.csv(
    as.data.frame(loss_anno),
    file.path(results_dir, "MOF_lost_annotated.csv"),
    row.names=FALSE
)

pdf(file.path(results_dir, "MOF_gain_peak_annotation.pdf"))
plotAnnoPie(gain_anno)
dev.off()

pdf(file.path(results_dir, "MOF_loss_peak_annotation.pdf"))
plotAnnoPie(loss_anno)
dev.off()

gain_df <- as.data.frame(gain_anno)
loss_df <- as.data.frame(loss_anno)

gain_genes <- unique(na.omit(gain_df$geneId))
loss_genes <- unique(na.omit(loss_df$geneId))

if (requireNamespace("clusterProfiler", quietly=TRUE) &&
    requireNamespace("enrichplot", quietly=TRUE)) {
    library(clusterProfiler)
    library(enrichplot)

    gain_go <- enrichGO(
        gene          = gain_genes,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05
    )

    loss_go <- enrichGO(
        gene          = loss_genes,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05
    )

    write.csv(as.data.frame(gain_go), file.path(results_dir, "MOF_gain_GO_BP.csv"), row.names=FALSE)
    write.csv(as.data.frame(loss_go), file.path(results_dir, "MOF_loss_GO_BP.csv"), row.names=FALSE)

    pdf(file.path(results_dir, "MOF_gain_GO.pdf"))
    if (nrow(as.data.frame(gain_go)) > 0) {
        print(dotplot(gain_go, showCategory=15) +
              ggtitle("Processes associated with gained MOF"))
    } else {
        plot.new()
        text(0.5, 0.5, "No significant GO terms for gained MOF")
    }
    dev.off()

    pdf(file.path(results_dir, "MOF_loss_GO.pdf"))
    if (nrow(as.data.frame(loss_go)) > 0) {
        print(dotplot(loss_go, showCategory=15) +
              ggtitle("Processes associated with lost MOF"))
    } else {
        plot.new()
        text(0.5, 0.5, "No significant GO terms for lost MOF")
    }
    dev.off()
} else {
    message("clusterProfiler/enrichplot not installed; skipping GO enrichment and plots.")
}
