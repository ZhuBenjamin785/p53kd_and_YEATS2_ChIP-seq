
args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: Rscript chipseeker_fasterqchip.r [--output-dir DIR] <peak-file> [<peak-file> ...]",
  "Annotates hg38 MACS3 narrowPeak/broadPeak files and writes tables and summary plots.",
  sep = "\n"
)

output_dir <- NULL
peak_files <- character()
index <- 1L
while (index <= length(args)) {
  if (args[[index]] == "--output-dir") {
    if (index == length(args)) stop("--output-dir requires a directory.\n", usage, call. = FALSE)
    output_dir <- args[[index + 1L]]
    index <- index + 2L
  } else if (args[[index]] %in% c("-h", "--help")) {
    cat(usage, "\n")
    quit(status = 0L)
  } else {
    peak_files <- c(peak_files, args[[index]])
    index <- index + 1L
  }
}

if (!length(peak_files)) stop(usage, call. = FALSE)
missing_files <- peak_files[!file.exists(peak_files)]
if (length(missing_files)) {
  stop("Peak file(s) not found:\n", paste(missing_files, collapse = "\n"), call. = FALSE)
}
empty_files <- peak_files[file.info(peak_files)$size == 0]
if (length(empty_files)) {
  stop("Peak file(s) are empty:\n", paste(empty_files, collapse = "\n"), call. = FALSE)
}

required_packages <- c(
  "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "ggplot2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Required R package(s) are missing: ", paste(missing_packages, collapse = ", "),
    "\nActivate the chipseeker Conda environment before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
})

if (is.null(output_dir)) output_dir <- file.path(dirname(peak_files[[1L]]), "chipseeker")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_name <- function(peak_file) {
  name <- basename(peak_file)
  name <- sub("_peaks\\.(narrowPeak|broadPeak)$", "", name)
  sub("\\.(narrowPeak|broadPeak)$", "", name)
}

annotations <- list()
summary_rows <- list()
for (peak_file in peak_files) {
  sample <- sample_name(peak_file)
  annotation <- annotatePeak(
    peak = peak_file,
    tssRegion = c(-3000, 3000),
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb = "org.Hs.eg.db",
    verbose = FALSE
  )
  annotations[[sample]] <- annotation

  annotated_table <- as.data.frame(annotation)
  table_file <- file.path(output_dir, paste0(sample, "_peaks_annotated.tsv"))
  write.table(annotated_table, table_file, sep = "\t", quote = FALSE, row.names = FALSE)

  annotation_counts <- as.data.frame(table(annotated_table$annotation), stringsAsFactors = FALSE)
  names(annotation_counts) <- c("annotation", "peak_count")
  annotation_counts$sample <- sample
  summary_rows[[sample]] <- annotation_counts[, c("sample", "annotation", "peak_count")]
  message("Wrote: ", table_file)
}

summary_table <- do.call(rbind, summary_rows)
write.table(
  summary_table,
  file.path(output_dir, "annotation_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

pdf(file.path(output_dir, "chipseeker_annotation_plots.pdf"), width = 11, height = 8.5)
print(plotAnnoBar(annotations))
for (sample in names(annotations)) {
  print(plotAnnoPie(annotations[[sample]], main = sample))
  print(plotDistToTSS(annotations[[sample]]) + ggtitle(sample))
}
dev.off()
message("Wrote ChIPseeker results to: ", normalizePath(output_dir))
