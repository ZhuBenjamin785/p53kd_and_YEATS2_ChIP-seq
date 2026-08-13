
packages <- c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "ggplot2")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(ggplot2)
})

peak_files <- c(
  Scramble1 = "macs3_results/peaks/Scramble_H4K16ac_rep1/Scramble_H4K16ac_rep1_peaks.broadPeak",
  Scramble2 = "macs3_results/peaks/Scramble_H4K16ac_rep2/Scramble_H4K16ac_rep2_peaks.broadPeak",
  KD1 = "macs3_results/peaks/YEATS2KD_H4K16ac_rep1/YEATS2KD_H4K16ac_rep1_peaks.broadPeak",
  KD2 = "macs3_results/peaks/YEATS2KD_H4K16ac_rep2/YEATS2KD_H4K16ac_rep2_peaks.broadPeak"
)
bad_files <- peak_files[!file.exists(peak_files) | file.info(peak_files)$size == 0]
if (length(bad_files)) stop("Missing or empty peak files:\n", paste(bad_files, collapse = "\n"))

output_dir <- "macs3_results/chipseeker_visualizations"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
input_counts <- vapply(peak_files, function(x) length(readLines(x, warn = FALSE)), integer(1))

annotations <- lapply(peak_files, function(peak_file) {
  message("Annotating: ", peak_file)
  suppressWarnings(annotatePeak(
    peak_file,
    tssRegion = c(-3000, 3000),
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb = "org.Hs.eg.db",
    verbose = FALSE
  ))
})

pdf_file <- file.path(output_dir, "chipseeker_annotation_plots.pdf")
pdf(pdf_file, width = 11, height = 8.5, onefile = TRUE)
print(plotAnnoBar(annotations))
for (sample_name in names(annotations)) {
  plotAnnoPie(annotations[[sample_name]], main = sample_name)
  print(plotDistToTSS(annotations[[sample_name]]) + ggtitle(sample_name))
}
dev.off()

gene_ids <- lapply(annotations, function(x) {
  ids <- as.character(as.data.frame(x)$geneId)
  sort(unique(ids[!is.na(ids) & nzchar(ids)]))
})
write_genes <- function(ids, filename) write.table(
  data.frame(geneId = ids), filename, sep = "\t", quote = FALSE, row.names = FALSE
)
write_genes(unique(c(gene_ids$Scramble1, gene_ids$Scramble2)),
            file.path(output_dir, "scramble_gene_ids.txt"))
write_genes(unique(c(gene_ids$KD1, gene_ids$KD2)),
            file.path(output_dir, "yeats2kd_gene_ids.txt"))

summary_table <- data.frame(
  sample = names(annotations),
  input_peak_count = unname(input_counts),
  annotated_peak_count = vapply(annotations, function(x) nrow(as.data.frame(x)), integer(1)),
  unique_gene_count = vapply(gene_ids, length, integer(1))
)
write.table(summary_table, file.path(output_dir, "annotation_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote plots and tables to: ", output_dir)


combined_plot <- (
  p1 | p2 | p3 | p4
) / (
  t1 | t2 | t3 | t4
)

ggsave(
  "chipseeker_annotation_combined.pdf",
  plot = combined_plot,
  width = 24,
  height = 12,
  units = "in"
)