
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    "Usage: Rscript chipseekerannotation.r <peak-file> [<peak-file> ...]",
    call. = FALSE
  )
}

missing_files <- args[!file.exists(args)]
if (length(missing_files) > 0) {
  stop(
    "Peak file(s) not found:\n", paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

empty_files <- args[file.info(args)$size == 0]
if (length(empty_files) > 0) {
  stop(
    "Peak file(s) are empty:\n", paste(empty_files, collapse = "\n"),
    call. = FALSE
  )
}

required_packages <- c(
  "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Required R package(s) are missing: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them in the chipseeker Conda environment before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})

for (peak_file in args) {
  annotation <- annotatePeak(
    peak = peak_file,
    tssRegion = c(-3000, 3000),
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb = "org.Hs.eg.db"
  )

  output_file <- file.path(
    dirname(peak_file),
    paste0(tools::file_path_sans_ext(basename(peak_file)), "_annotated.tsv")
  )

  write.table(
    as.data.frame(annotation),
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  message("Wrote: ", output_file)
}
