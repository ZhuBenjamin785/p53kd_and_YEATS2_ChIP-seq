
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
})

out_root <- file.path(getwd(), "chipseq_summary_plots")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
args <- commandArgs(trailingOnly = TRUE)
requested_dataset <- if (length(args)) args[[1]] else "all"

datasets <- list(
  MOF = list(
    peaks = c(
      `Scramble MOF 1` = "mof_macs3_results/peaks/Scr_MOF_rep1/Scr_MOF_rep1_peaks.broadPeak",
      `Scramble MOF 2` = "mof_macs3_results/peaks/Scr_MOF_rep2/Scr_MOF_rep2_peaks.broadPeak",
      `p53sh MOF 1` = "mof_macs3_results/peaks/p53sh_MOF_rep1/p53sh_MOF_rep1_peaks.broadPeak",
      `p53sh MOF 2` = "mof_macs3_results/peaks/p53sh_MOF_rep2/p53sh_MOF_rep2_peaks.broadPeak"
    ),
    diff = "mof_macs3_results/diffbind_results/DiffBind_all_peaks.csv"
  ),
  p53KD = list(
    peaks = c(
      `Scramble 1` = "macs3_results_p53kd/peaks/Scramble_H4K16ac_rep1/Scramble_H4K16ac_rep1_peaks.broadPeak",
      `Scramble 2` = "macs3_results_p53kd/peaks/Scramble_H4K16ac_rep2/Scramble_H4K16ac_rep2_peaks.broadPeak",
      `p53 KD 1` = "macs3_results_p53kd/peaks/p53KD_H4K16ac_rep1/p53KD_H4K16ac_rep1_peaks.broadPeak",
      `p53 KD 2` = "macs3_results_p53kd/peaks/p53KD_H4K16ac_rep2/p53KD_H4K16ac_rep2_peaks.broadPeak"
    ),
    diff = "diffbind_results/DiffBind_all_peaks.csv"
  ),
  YEATS2KD = list(
    peaks = c(
      `Scramble 1` = "macs3_results_yeats2/peaks/Scramble_H4K16ac_rep1/Scramble_H4K16ac_rep1_peaks.broadPeak",
      `Scramble 2` = "macs3_results_yeats2/peaks/Scramble_H4K16ac_rep2/Scramble_H4K16ac_rep2_peaks.broadPeak",
      `YEATS2 KD 1` = "macs3_results_yeats2/peaks/YEATS2KD_H4K16ac_rep1/YEATS2KD_H4K16ac_rep1_peaks.broadPeak",
      `YEATS2 KD 2` = "macs3_results_yeats2/peaks/YEATS2KD_H4K16ac_rep2/YEATS2KD_H4K16ac_rep2_peaks.broadPeak"
    ),
    diff = "diffbind_results_yeats2/DiffBind_YEATS2_all_peaks.csv"
  ),
  p53_0hr = list(
    peaks = c(`p53, 0 hr` = "fastqchip_macs3_results/peaks/p53_0hr/p53_0hr_peaks.narrowPeak"),
    diff = NA_character_
  )
)

detail_category <- function(annotation) {
  result <- as.character(annotation)
  result[grepl("^Promoter \\(<=1kb\\)", result)] <- "Promoter (<=1 kb)"
  result[grepl("^Promoter \\(1-2kb\\)", result)] <- "Promoter (1-2 kb)"
  result[grepl("^Promoter \\(2-3kb\\)", result)] <- "Promoter (2-3 kb)"
  result[grepl("^5' UTR", result)] <- "5' UTR"
  result[grepl("^3' UTR", result)] <- "3' UTR"
  result[grepl("^1st Exon", result)] <- "First exon"
  result[grepl("^Exon", result)] <- "Other exon"
  result[grepl("^1st Intron", result)] <- "First intron"
  result[grepl("^Intron", result)] <- "Other intron"
  result[grepl("^Downstream", result)] <- "Downstream"
  result[grepl("^Distal Intergenic", result)] <- "Distal intergenic"
  result
}

if (requested_dataset != "all" && !requested_dataset %in% names(datasets)) {
  stop("Unknown dataset: ", requested_dataset, ". Choose: ", paste(names(datasets), collapse = ", "))
}
datasets <- if (requested_dataset == "all") datasets else datasets[requested_dataset]

for (dataset in names(datasets)) {
  cfg <- datasets[[dataset]]
  peak_files <- cfg$peaks
  bad <- peak_files[!file.exists(peak_files) | file.info(peak_files)$size == 0]
  if (length(bad)) stop(dataset, " missing peaks:\n", paste(bad, collapse = "\n"))
  outdir <- file.path(out_root, dataset)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  annotations <- lapply(peak_files, function(path) suppressWarnings(annotatePeak(
    path, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE
  )))
  rows <- do.call(rbind, Map(function(x, sample) {
    data.frame(dataset = dataset, sample = sample, as.data.frame(x), stringsAsFactors = FALSE)
  }, annotations, names(annotations)))
  rows$annotation_category <- detail_category(rows$annotation)
  write.csv(rows, file.path(outdir, "chipseeker_peak_annotations.csv"), row.names = FALSE)

  summary <- as.data.frame(table(rows$sample, rows$annotation_category), stringsAsFactors = FALSE)
  names(summary) <- c("Sample", "Category", "Count")
  summary <- summary[summary$Count > 0, ]
  totals <- aggregate(Count ~ Sample, summary, sum)
  names(totals)[2] <- "Total"
  summary <- merge(summary, totals, by = "Sample")
  summary$Percentage <- 100 * summary$Count / summary$Total
  write.csv(summary, file.path(outdir, "peak_annotation_summary.csv"), row.names = FALSE)

  if (!is.na(cfg$diff) && file.exists(cfg$diff)) {
    diff <- read.csv(cfg$diff, check.names = FALSE)
    required <- c("seqnames", "start", "end")
    if (all(required %in% names(diff)) && nrow(diff)) {
      gr <- makeGRangesFromDataFrame(diff, keep.extra.columns = TRUE,
        seqnames.field = "seqnames", start.field = "start", end.field = "end")
      diff_anno <- suppressWarnings(annotatePeak(gr, TxDb = txdb,
        annoDb = "org.Hs.eg.db", verbose = FALSE))
      diff_df <- as.data.frame(diff_anno)
      diff_df$annotation_category <- detail_category(diff_df$annotation)
      write.csv(diff_df, file.path(outdir, "diffbind_peak_annotations.csv"), row.names = FALSE)
    }
  }
}
message("ChIPseeker annotation tables written to: ", out_root)
