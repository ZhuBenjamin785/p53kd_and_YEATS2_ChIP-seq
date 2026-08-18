
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: run_go_enrichment.R PROJECT_DIR ANALYSIS_DIR")
}
project_dir <- normalizePath(args[[1]])
analysis_dir <- normalizePath(args[[2]])
tables_dir <- file.path(analysis_dir, "tables")

membership <- read.csv(file.path(tables_dir, "annotated_gene_overlap_membership.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
membership$annotated_gene_id <- sub("\\.[0-9]+$", "", membership$annotated_gene_id)
for (column in c("in_YEATS2KD", "in_p53KD", "shared")) {
  membership[[column]] <- tolower(as.character(membership[[column]])) == "true"
}

annotation_paths <- c(
  file.path(project_dir, "shared/chipseq_summary_plots/YEATS2KD/diffbind_peak_annotations.csv"),
  file.path(project_dir, "shared/chipseq_summary_plots/p53KD/diffbind_peak_annotations.csv")
)
universe <- unique(unlist(lapply(annotation_paths, function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  sub("\\.[0-9]+$", "", na.omit(x$ENSEMBL))
})))
universe <- universe[nzchar(universe) & universe != "NA"]

gene_sets <- list(
  YEATS2KD_only = membership$annotated_gene_id[membership$in_YEATS2KD & !membership$in_p53KD],
  p53KD_only = membership$annotated_gene_id[!membership$in_YEATS2KD & membership$in_p53KD],
  shared = membership$annotated_gene_id[membership$shared]
)

empty_result <- data.frame(
  ID = character(), Description = character(), GeneRatio = character(),
  BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
  qvalue = numeric(), geneID = character(), Count = integer(),
  check.names = FALSE
)
status_rows <- list()

for (set_name in names(gene_sets)) {
  genes <- unique(gene_sets[[set_name]])
  genes <- genes[nzchar(genes) & genes != "NA"]
  output_path <- file.path(tables_dir, paste0("GO_BP_", set_name, ".csv"))
  if (length(genes) < 5) {
    write.csv(empty_result, output_path, row.names = FALSE)
    status_rows[[set_name]] <- data.frame(
      gene_set = set_name, input_genes = length(genes), mapped_genes = length(genes),
      significant_terms = 0,
      status = "Too few genes for a reliable over-representation test (minimum 5)"
    )
    next
  }
  result <- enrichGO(
    gene = genes, universe = universe, OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5,
    maxGSSize = 500, readable = TRUE
  )
  result_df <- as.data.frame(result)
  write.csv(result_df, output_path, row.names = FALSE)
  mapped <- unique(unlist(strsplit(result_df$geneID, "/", fixed = TRUE)))
  mapped <- mapped[nzchar(mapped)]
  status_rows[[set_name]] <- data.frame(
    gene_set = set_name, input_genes = length(genes), mapped_genes = length(mapped),
    significant_terms = sum(result_df$p.adjust < 0.05, na.rm = TRUE), status = "Test completed"
  )
}

write.csv(do.call(rbind, status_rows),
          file.path(tables_dir, "GO_enrichment_status.csv"), row.names = FALSE)
