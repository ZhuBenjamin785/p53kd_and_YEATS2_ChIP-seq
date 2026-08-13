
packages <- c(
  "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db",
  "ggplot2", "patchwork"
)
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(ggplot2)
  library(patchwork)
})

peak_files <- c(
  p53kd2 = "p53kd_consensus.bed",
  scramble = "macs3_results_p53kd/peaks/Scramble_H4K16ac_rep1/Scramble_H4K16ac_rep1_peaks.broadPeak"
)
bad_files <- peak_files[!file.exists(peak_files) | file.info(peak_files)$size == 0]
if (length(bad_files)) stop("Missing or empty peak files:\n", paste(bad_files, collapse = "\n"))

output_dir <- "macs3_results_p53kd/chipseeker_visualizations"
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

gene_ids <- lapply(annotations, function(x) {
  ids <- as.character(as.data.frame(x)$geneId)
  sort(unique(ids[!is.na(ids) & nzchar(ids)]))
})
write_genes <- function(ids, filename) write.table(
  data.frame(geneId = ids), filename, sep = "\t", quote = FALSE, row.names = FALSE
)
for (sample_name in names(gene_ids)) {
  safe_name <- gsub("[^[:alnum:]_-]+", "_", tolower(sample_name))
  write_genes(
    gene_ids[[sample_name]],
    file.path(output_dir, paste0(safe_name, "_gene_ids.txt"))
  )
}

summary_table <- data.frame(
  sample = names(annotations),
  input_peak_count = unname(input_counts),
  annotated_peak_count = vapply(annotations, function(x) nrow(as.data.frame(x)), integer(1)),
  unique_gene_count = vapply(gene_ids, length, integer(1))
)
write.table(summary_table, file.path(output_dir, "annotation_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
pie_plot <- function(annotation, sample_name) {
  wrap_elements(full = ~plotAnnoPie(annotation, main = sample_name))
}

pie_plots <- Map(pie_plot, annotations, names(annotations))
tss_plots <- Map(
  function(annotation, sample_name) {
    plotDistToTSS(annotation) + ggtitle(sample_name)
  },
  annotations,
  names(annotations)
)

combined_plot <- wrap_plots(pie_plots, nrow = 1) /
  wrap_plots(tss_plots, nrow = 1)

combined_pdf <- file.path(output_dir, "chipseeker_annotation_combined.pdf")
ggsave(
  filename = combined_pdf,
  plot = combined_plot,
  width = max(8, 6 * length(annotations)),
  height = 10,
  units = "in"
)

pdf_file <- file.path(output_dir, "chipseeker_annotation_plots.pdf")
pdf(pdf_file, width = 11, height = 8.5, onefile = TRUE)
print(combined_plot)
print(plotAnnoBar(annotations))
for (sample_name in names(annotations)) {
  print(pie_plots[[sample_name]])
  print(tss_plots[[sample_name]])
}
dev.off()

output_files <- c(pdf_file, combined_pdf, file.path(output_dir, "annotation_summary.tsv"))
bad_outputs <- output_files[!file.exists(output_files) | file.info(output_files)$size == 0]
if (length(bad_outputs)) {
  stop("Failed to create output files:\n", paste(bad_outputs, collapse = "\n"))
}
message("Wrote combined plot to: ", combined_pdf)
message("Wrote combined and individual plots to: ", pdf_file)
