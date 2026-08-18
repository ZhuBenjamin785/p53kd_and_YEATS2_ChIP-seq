#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: split_diffbind_genes.R DIFFBIND_CSV OUTPUT_DIR FDR_CUTOFF LOGFC_CUTOFF")
}

input_csv <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
fdr_cutoff <- as.numeric(args[[3]])
logfc_cutoff <- as.numeric(args[[4]])
if (!is.finite(fdr_cutoff) || fdr_cutoff <= 0 || fdr_cutoff >= 1) {
  stop("FDR_CUTOFF must be > 0 and < 1")
}
if (!is.finite(logfc_cutoff) || logfc_cutoff < 0) {
  stop("LOGFC_CUTOFF must be >= 0")
}

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
diffbind <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("seqnames", "start", "end", "Fold", "FDR")
missing <- setdiff(required, names(diffbind))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

diffbind$Fold <- as.numeric(diffbind$Fold)
diffbind$FDR <- as.numeric(diffbind$FDR)
significant <- diffbind[is.finite(diffbind$Fold) & is.finite(diffbind$FDR) &
                          diffbind$FDR <= fdr_cutoff &
                          abs(diffbind$Fold) >= logfc_cutoff, , drop = FALSE]

if (!nrow(significant)) {
  stop("No peaks passed FDR <= ", fdr_cutoff, " and abs(Fold) >= ", logfc_cutoff)
}
write.csv(significant, file.path(output_dir, "filtered_diffbind_peaks.csv"), row.names = FALSE)

annotate_diffbind_peaks <- function(peaks) {
  gr <- makeGRangesFromDataFrame(
    peaks,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )
  annotated <- as.data.frame(suppressWarnings(annotatePeak(
    gr,
    tssRegion = c(-3000, 3000),
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb = "org.Hs.eg.db",
    verbose = FALSE
  )))
  annotated$direction <- ifelse(annotated$Fold > 0, "up", "down")
  annotated
}

annotated <- annotate_diffbind_peaks(significant)
up <- annotated[annotated$direction == "up", , drop = FALSE]
down <- annotated[annotated$direction == "down", , drop = FALSE]

write.csv(annotated, file.path(output_dir, "significant_annotated.csv"), row.names = FALSE)
write.csv(up, file.path(output_dir, "significant_up_peaks.csv"), row.names = FALSE)
write.csv(down, file.path(output_dir, "significant_down_peaks.csv"), row.names = FALSE)

write_gene_table <- function(annotation, filename) {
  columns <- intersect(c("geneId", "ENSEMBL", "SYMBOL", "GENENAME"), names(annotation))
  if (!length(columns)) {
    write.table(data.frame(), filename, sep = "\t", quote = FALSE, row.names = FALSE)
    return(invisible(NULL))
  }
  genes <- annotation[, columns, drop = FALSE]
  for (column in columns) {
    genes[[column]] <- as.character(genes[[column]])
    genes[[column]][is.na(genes[[column]]) | genes[[column]] %in% c("", "NA", ".")] <- NA
  }
  genes <- genes[!ifelse("geneId" %in% columns, is.na(genes$geneId),
                         apply(is.na(genes), 1, all)), , drop = FALSE]
  genes <- unique(genes)
  write.table(genes, filename, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_gene_table(up, file.path(output_dir, "up_genes.tsv"))
write_gene_table(down, file.path(output_dir, "down_genes.tsv"))

gene_value <- function(values) {
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values) & !values %in% c("NA", ".")]
  if (length(values)) values[[1]] else NA_character_
}

make_gene_list <- function(annotation) {
  if (!"geneId" %in% names(annotation)) return(data.frame())
  gene_ids <- as.character(annotation$geneId)
  valid <- !is.na(gene_ids) & nzchar(gene_ids) & !gene_ids %in% c("NA", ".")
  annotation <- annotation[valid, , drop = FALSE]
  if (!nrow(annotation)) return(data.frame())

  rows <- lapply(split(annotation, as.character(annotation$geneId)), function(peaks) {
    data.frame(
      geneId = gene_value(peaks$geneId),
      ENSEMBL = if ("ENSEMBL" %in% names(peaks)) gene_value(peaks$ENSEMBL) else NA_character_,
      SYMBOL = if ("SYMBOL" %in% names(peaks)) gene_value(peaks$SYMBOL) else NA_character_,
      GENENAME = if ("GENENAME" %in% names(peaks)) gene_value(peaks$GENENAME) else NA_character_,
      direction = if (length(unique(peaks$direction)) == 1) unique(peaks$direction) else "mixed",
      peak_count = nrow(peaks),
      strongest_abs_fold = max(abs(peaks$Fold), na.rm = TRUE),
      minimum_FDR = min(peaks$FDR, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

gene_list <- make_gene_list(annotated)
up_gene_list <- gene_list[gene_list$direction == "up", , drop = FALSE]
down_gene_list <- gene_list[gene_list$direction == "down", , drop = FALSE]
write.csv(gene_list, file.path(output_dir, "significant_gene_list.csv"), row.names = FALSE)
write.csv(up_gene_list, file.path(output_dir, "up_gene_list.csv"), row.names = FALSE)
write.csv(down_gene_list, file.path(output_dir, "down_gene_list.csv"), row.names = FALSE)

make_gene_peak_ranking <- function(annotation) {
  if (!"geneId" %in% names(annotation)) return(data.frame())
  gene_ids <- as.character(annotation$geneId)
  valid <- !is.na(gene_ids) & nzchar(gene_ids) & !gene_ids %in% c("NA", ".")
  annotation <- annotation[valid, , drop = FALSE]
  if (!nrow(annotation)) return(data.frame())

  rows <- lapply(split(annotation, as.character(annotation$geneId)), function(peaks) {
    highest <- peaks[which.max(peaks$Fold), , drop = FALSE][1, ]
    lowest <- peaks[which.min(peaks$Fold), , drop = FALSE][1, ]
    data.frame(
      geneId = gene_value(peaks$geneId),
      ENSEMBL = if ("ENSEMBL" %in% names(peaks)) gene_value(peaks$ENSEMBL) else NA_character_,
      SYMBOL = if ("SYMBOL" %in% names(peaks)) gene_value(peaks$SYMBOL) else NA_character_,
      GENENAME = if ("GENENAME" %in% names(peaks)) gene_value(peaks$GENENAME) else NA_character_,
      peak_count = nrow(peaks),
      highest_logFC = highest$Fold,
      highest_logFC_chrom = highest$seqnames,
      highest_logFC_start = highest$start,
      highest_logFC_end = highest$end,
      lowest_logFC = lowest$Fold,
      lowest_logFC_chrom = lowest$seqnames,
      lowest_logFC_start = lowest$start,
      lowest_logFC_end = lowest$end,
      stringsAsFactors = FALSE
    )
  })
  ranking <- do.call(rbind, rows)
  ranking$rank_highest_logFC <- rank(-ranking$highest_logFC, ties.method = "min")
  ranking$rank_lowest_logFC <- rank(ranking$lowest_logFC, ties.method = "min")
  ranking$rank_highest_abs_logFC <- rank(-pmax(abs(ranking$highest_logFC),
                                               abs(ranking$lowest_logFC)),
                                         ties.method = "min")
  ranking[order(ranking$rank_highest_abs_logFC, ranking$rank_highest_logFC), ]
}

gene_peak_ranking <- make_gene_peak_ranking(annotated)
write.csv(gene_peak_ranking,
          file.path(output_dir, "gene_peak_logFC_ranking.csv"),
          row.names = FALSE)

summary <- data.frame(
  input_table = input_csv,
  fdr_cutoff = fdr_cutoff,
  logfc_cutoff = logfc_cutoff,
  input_peaks = nrow(diffbind),
  significant_peaks = nrow(annotated),
  up_peaks = nrow(up),
  down_peaks = nrow(down),
  up_gene_ids = nrow(up_gene_list),
  down_gene_ids = nrow(down_gene_list)
)
write.csv(summary, file.path(output_dir, "split_summary.csv"), row.names = FALSE)
message("Wrote filtered, annotated, up/down, and ranked gene tables to: ", output_dir)
