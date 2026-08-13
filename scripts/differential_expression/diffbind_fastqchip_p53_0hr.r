
suppressPackageStartupMessages({
  library(DiffBind)
  library(BiocParallel)
})

results_dir <- "diffbind_results_fastqchip_p53_0hr"
plots_dir <- file.path(results_dir, "plots")
tables_dir <- file.path(results_dir, "tables")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

n_threads <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
register(MulticoreParam(workers = n_threads, progressbar = FALSE), default = TRUE)

peak_file <- "fastqchip_macs3_results/peaks/p53_0hr/p53_0hr_peaks.narrowPeak"

samples <- data.frame(
  SampleID = c("p53_0hr_rep1", "p53_0hr_rep2"),
  Condition = c("p53_0hr", "p53_0hr"),
  Replicate = c(1, 2),
  bamReads = c(
    "fastqchip_bamfiles/human/SRR5944063.sorted.bam",
    "fastqchip_bamfiles/human/SRR5944064.sorted.bam"
  ),
  bamControl = c(
    "fastqchip_bamfiles/human/SRR5944081.sorted.bam",
    "fastqchip_bamfiles/human/SRR5944082.sorted.bam"
  ),
  Spikein = c(
    "fastqchip_bamfiles/dm6/SRR5944063.sorted.bam",
    "fastqchip_bamfiles/dm6/SRR5944064.sorted.bam"
  ),
  Peaks = c(peak_file, peak_file),
  PeakCaller = c("macs", "macs"),
  stringsAsFactors = FALSE
)

write.csv(samples, file.path(results_dir, "metadata.csv"), row.names = FALSE)
write.table(
  samples, file.path(results_dir, "metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

required_files <- unique(c(
  samples$bamReads, samples$bamControl, samples$Spikein, samples$Peaks
))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

for (bam in c(samples$bamReads, samples$bamControl, samples$Spikein)) {
  if (!file.exists(paste0(bam, ".bai"))) {
    stop("Missing BAM index: ", paste0(bam, ".bai"))
  }
}

db <- dba(sampleSheet = samples)
db$config$cores <- n_threads
db$config$RunParallel <- n_threads > 1L

db <- dba.blacklist(
  db,
  blacklist = DBA_BLACKLIST_HG38,
  greylist = FALSE,
  cores = n_threads
)
db$config$doGreylist <- FALSE
db$config$doBlacklist <- FALSE

db <- dba.count(
  db,
  summits = 250,
  bParallel = n_threads > 1L,
  bSubControl = FALSE
)

db <- dba.normalize(db, normalize = DBA_NORM_LIB, spikein = TRUE)

write.csv(
  as.data.frame(dba.show(db)),
  file.path(tables_dir, "sample_summary.csv"),
  row.names = FALSE
)

pdf(file.path(plots_dir, "DiffBind_PCA.pdf"), width = 8, height = 6)
print(dba.plotPCA(db, attributes = DBA_ID))
dev.off()

pdf(file.path(plots_dir, "DiffBind_correlation_heatmap.pdf"), width = 8, height = 8)
print(plot(db))
dev.off()

heatmap_file <- file.path(plots_dir, "DiffBind_binding_heatmap.pdf")
pdf(heatmap_file, width = 9, height = 8)
heatmap_ok <- tryCatch({
  print(dba.plotHeatmap(db, correlations = FALSE))
  TRUE
}, error = function(e) {
  plot.new()
  text(0.5, 0.5, paste("Binding heatmap unavailable:", conditionMessage(e)))
  FALSE
})
dev.off()
writeLines(
  if (heatmap_ok) "dba.plotHeatmap completed" else "dba.plotHeatmap was unavailable; see PDF message",
  file.path(tables_dir, "heatmap_status.txt")
)

saveRDS(db, file.path(results_dir, "diffbind_p53_0hr_normalized.rds"))

message("DiffBind results written to: ", results_dir)
message("No differential contrast was run: this dataset contains one biological condition.")
