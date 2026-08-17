#!/usr/bin/env Rscript

# Integrate p53KD RNA-seq differential expression with H4K16ac DiffBind peaks.
# Contrast convention for both datasets: p53KD/shp53 relative to Scramble/shLacZ.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- "/gpfs/projects/b1042/LauberthLab/BenFolder"
rna_file <- if (length(args) >= 1) args[[1]] else file.path(
  project_dir, "shared/rna_seq_dea/shp53_vs_shLacZ_0hr/significant_results.csv"
)
chip_file <- if (length(args) >= 2) args[[2]] else file.path(
  project_dir, "p53kdH4K16ac/diffbind_results/split_genes/significant_annotated.csv"
)
output_dir <- if (length(args) >= 3) args[[3]] else file.path(
  project_dir, "shared/rna_chip_integration/p53KD_H4K16ac_vs_RNAseq"
)

rna_padj_cutoff <- 0.05
chip_fdr_cutoff <- 0.05
category_levels <- c("H4K16ac loss + RNA down", "H4K16ac gain + RNA up",
                     "H4K16ac loss + RNA up", "H4K16ac gain + RNA down")
category_colors <- c(
  "H4K16ac loss + RNA down" = "#2166AC",
  "H4K16ac gain + RNA up" = "#B2182B",
  "H4K16ac loss + RNA up" = "#67A9CF",
  "H4K16ac gain + RNA down" = "#EF8A62"
)

find_column <- function(data, candidates, label, required = TRUE) {
  hit <- names(data)[tolower(names(data)) %in% tolower(candidates)]
  if (length(hit) == 0) {
    if (required) stop("Missing ", label, " column. Expected one of: ", paste(candidates, collapse = ", "))
    return(NA_character_)
  }
  hit[[1]]
}

write_csv <- function(data, path) {
  write.csv(data, path, row.names = FALSE, na = "")
}

theme_publication <- function() {
  theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey30"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
}

if (!file.exists(rna_file)) stop("RNA-seq file not found: ", rna_file)
if (!file.exists(chip_file)) stop("H4K16ac file not found: ", chip_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rna_raw <- read.csv(rna_file, check.names = FALSE, stringsAsFactors = FALSE)
chip_raw <- read.csv(chip_file, check.names = FALSE, stringsAsFactors = FALSE)

rna_gene_col <- find_column(rna_raw, c("gene_name", "Gene", "Gene symbol", "SYMBOL"), "RNA gene symbol")
rna_fc_col <- find_column(rna_raw, c("log2FoldChange", "log2FC", "Fold"), "RNA log2 fold-change")
rna_padj_col <- find_column(rna_raw, c("padj", "FDR", "adjusted_pvalue"), "RNA adjusted P-value")

chip_gene_col <- find_column(chip_raw, c("SYMBOL", "Gene", "Gene symbol", "gene_name"), "ChIP gene symbol")
chip_fc_col <- find_column(chip_raw, c("Fold", "log2FoldChange", "log2FC"), "H4K16ac Fold")
chip_fdr_col <- find_column(chip_raw, c("FDR", "padj", "adjusted_pvalue"), "H4K16ac FDR")
chip_annotation_col <- find_column(chip_raw, c("annotation", "Peak annotation", "annotation_category"), "peak annotation")
chip_distance_col <- find_column(chip_raw, c("distanceToTSS", "Distance to TSS", "distance_to_tss"), "distance to TSS")

optional_chip_col <- function(candidates) find_column(chip_raw, candidates, paste(candidates, collapse = "/"), FALSE)
seq_col <- optional_chip_col(c("seqnames", "Chr", "chromosome"))
start_col <- optional_chip_col(c("start", "Start"))
end_col <- optional_chip_col(c("end", "End"))

rna_sig <- rna_raw %>%
  transmute(
    Gene = toupper(trimws(as.character(.data[[rna_gene_col]]))),
    RNA_log2FC = suppressWarnings(as.numeric(.data[[rna_fc_col]])),
    RNA_padj = suppressWarnings(as.numeric(.data[[rna_padj_col]]))
  ) %>%
  filter(!is.na(Gene), Gene != "", !is.na(RNA_log2FC), !is.na(RNA_padj), RNA_padj < rna_padj_cutoff) %>%
  arrange(RNA_padj, desc(abs(RNA_log2FC))) %>%
  distinct(Gene, .keep_all = TRUE)

chip_sig <- chip_raw %>%
  transmute(
    Gene = toupper(trimws(as.character(.data[[chip_gene_col]]))),
    H4K16ac_Fold = suppressWarnings(as.numeric(.data[[chip_fc_col]])),
    H4K16ac_FDR = suppressWarnings(as.numeric(.data[[chip_fdr_col]])),
    Peak_annotation = as.character(.data[[chip_annotation_col]]),
    Distance_to_TSS = suppressWarnings(as.numeric(.data[[chip_distance_col]])),
    Peak_chr = if (!is.na(seq_col)) as.character(.data[[seq_col]]) else NA_character_,
    Peak_start = if (!is.na(start_col)) suppressWarnings(as.numeric(.data[[start_col]])) else NA_real_,
    Peak_end = if (!is.na(end_col)) suppressWarnings(as.numeric(.data[[end_col]])) else NA_real_
  ) %>%
  filter(!is.na(Gene), Gene != "", !is.na(H4K16ac_Fold), !is.na(H4K16ac_FDR), H4K16ac_FDR < chip_fdr_cutoff)

collapse_strongest_peak <- function(data) {
  data %>%
    arrange(Gene, H4K16ac_FDR, desc(abs(H4K16ac_Fold))) %>%
    group_by(Gene) %>%
    slice(1) %>%
    ungroup()
}

safe_correlations <- function(data) {
  if (nrow(data) < 3 || sd(data$H4K16ac_Fold) == 0 || sd(data$RNA_log2FC) == 0) {
    return(tibble(Method = c("Pearson", "Spearman"), N = nrow(data),
                  Correlation = NA_real_, P_value = NA_real_))
  }
  pearson <- cor.test(data$H4K16ac_Fold, data$RNA_log2FC, method = "pearson")
  spearman <- suppressWarnings(cor.test(data$H4K16ac_Fold, data$RNA_log2FC,
                                        method = "spearman", exact = FALSE))
  tibble(
    Method = c("Pearson", "Spearman"),
    N = nrow(data),
    Correlation = c(unname(pearson$estimate), unname(spearman$estimate)),
    P_value = c(pearson$p.value, spearman$p.value)
  )
}

make_venn_plot <- function(rna_genes, chip_genes, path_prefix, title) {
  overlap <- length(intersect(rna_genes, chip_genes))
  rna_only <- length(setdiff(rna_genes, chip_genes))
  chip_only <- length(setdiff(chip_genes, rna_genes))
  theta <- seq(0, 2 * pi, length.out = 300)
  circles <- bind_rows(
    tibble(x = -0.65 + cos(theta), y = sin(theta), set = "RNA-seq"),
    tibble(x = 0.65 + cos(theta), y = sin(theta), set = "H4K16ac")
  )
  plot <- ggplot(circles, aes(x, y, group = set, fill = set)) +
    geom_polygon(alpha = 0.32, color = "grey25", linewidth = 0.6) +
    annotate("text", x = -0.85, y = 0, label = rna_only, size = 6, fontface = "bold") +
    annotate("text", x = 0, y = 0, label = overlap, size = 6, fontface = "bold") +
    annotate("text", x = 0.85, y = 0, label = chip_only, size = 6, fontface = "bold") +
    annotate("text", x = -0.9, y = 1.12, label = "RNA-seq", fontface = "bold") +
    annotate("text", x = 0.9, y = 1.12, label = "H4K16ac", fontface = "bold") +
    scale_fill_manual(values = c("RNA-seq" = "#4C78A8", "H4K16ac" = "#E45756")) +
    coord_equal(xlim = c(-1.8, 1.8), ylim = c(-1.35, 1.4), clip = "off") +
    labs(title = title, subtitle = "Significant gene-set overlap") +
    theme_void(base_size = 12) + theme(legend.position = "none", plot.title = element_text(face = "bold"))
  ggsave(paste0(path_prefix, ".pdf"), plot, width = 6, height = 5)
  ggsave(paste0(path_prefix, ".png"), plot, width = 6, height = 5, dpi = 300)
}

run_integration <- function(chip_peaks, analysis_name, promoter_only = FALSE) {
  analysis_dir <- file.path(output_dir, analysis_name)
  table_dir <- file.path(analysis_dir, "tables")
  figure_dir <- file.path(analysis_dir, "figures")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  strongest <- collapse_strongest_peak(chip_peaks)
  integrated <- inner_join(rna_sig, strongest, by = "Gene") %>%
    mutate(
      Biological_category = case_when(
        H4K16ac_Fold < 0 & RNA_log2FC < 0 ~ "H4K16ac loss + RNA down",
        H4K16ac_Fold > 0 & RNA_log2FC > 0 ~ "H4K16ac gain + RNA up",
        H4K16ac_Fold < 0 & RNA_log2FC > 0 ~ "H4K16ac loss + RNA up",
        H4K16ac_Fold > 0 & RNA_log2FC < 0 ~ "H4K16ac gain + RNA down",
        TRUE ~ "No directional change"
      ),
      Biological_category = factor(Biological_category, levels = category_levels)
    ) %>%
    filter(!is.na(Biological_category)) %>%
    arrange(Biological_category, H4K16ac_FDR, RNA_padj)

  final_table <- integrated %>%
    select(Gene, RNA_log2FC, RNA_padj, H4K16ac_Fold, H4K16ac_FDR,
           Peak_annotation, Distance_to_TSS, Peak_chr, Peak_start, Peak_end,
           Biological_category)

  counts <- final_table %>%
    count(Biological_category, .drop = FALSE, name = "Gene_count") %>%
    mutate(Concordance = ifelse(Biological_category %in%
      c("H4K16ac loss + RNA down", "H4K16ac gain + RNA up"), "Concordant", "Discordant"))
  correlations <- safe_correlations(final_table)

  write_csv(rna_sig, file.path(table_dir, "significant_rna_genes.csv"))
  write_csv(strongest, file.path(table_dir, "significant_h4k16ac_strongest_peak_per_gene.csv"))
  write_csv(final_table, file.path(table_dir, "integrated_gene_summary.csv"))
  write_csv(counts, file.path(table_dir, "category_counts.csv"))
  write_csv(correlations, file.path(table_dir, "correlation_statistics.csv"))
  for (category in category_levels) {
    filename <- c(
      "H4K16ac loss + RNA down" = "loss_down.csv",
      "H4K16ac gain + RNA up" = "gain_up.csv",
      "H4K16ac loss + RNA up" = "loss_up.csv",
      "H4K16ac gain + RNA down" = "gain_down.csv"
    )[[category]]
    write_csv(filter(final_table, Biological_category == category), file.path(table_dir, filename))
  }

  pearson_r <- correlations$Correlation[correlations$Method == "Pearson"]
  spearman_rho <- correlations$Correlation[correlations$Method == "Spearman"]
  subtitle <- sprintf("n = %d; Pearson r = %.3f; Spearman rho = %.3f",
                      nrow(final_table), pearson_r, spearman_rho)
  label_data <- final_table %>%
    mutate(label_score = abs(H4K16ac_Fold) + abs(RNA_log2FC)) %>%
    arrange(desc(label_score)) %>%
    slice_head(n = 15)

  scatter <- ggplot(final_table, aes(H4K16ac_Fold, RNA_log2FC, color = Biological_category)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(alpha = 0.72, size = 2) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
    geom_text_repel(data = label_data, aes(label = Gene), size = 3, max.overlaps = Inf,
                    box.padding = 0.35, min.segment.length = 0) +
    scale_color_manual(values = category_colors, drop = FALSE) +
    labs(
      title = if (promoter_only) "Promoter H4K16ac change versus RNA expression" else "H4K16ac change versus RNA expression",
      subtitle = subtitle,
      x = "H4K16ac DiffBind Fold (p53KD - Scramble)",
      y = "RNA-seq log2 fold-change (shp53 - shLacZ)",
      color = "Biological category"
    ) + theme_publication()
  ggsave(file.path(figure_dir, "h4k16ac_vs_rna_scatter.pdf"), scatter, width = 9, height = 7)
  ggsave(file.path(figure_dir, "h4k16ac_vs_rna_scatter.png"), scatter, width = 9, height = 7, dpi = 300)

  heatmap_data <- final_table %>%
    mutate(score = abs(H4K16ac_Fold) + abs(RNA_log2FC)) %>%
    arrange(desc(score)) %>%
    slice_head(n = 100) %>%
    select(Gene, RNA_log2FC, H4K16ac_Fold) %>%
    pivot_longer(c(RNA_log2FC, H4K16ac_Fold), names_to = "Dataset", values_to = "log2FC") %>%
    mutate(Dataset = recode(Dataset, RNA_log2FC = "RNA-seq", H4K16ac_Fold = "H4K16ac"))
  gene_order <- unique(heatmap_data$Gene)
  heatmap_data$Gene <- factor(heatmap_data$Gene, levels = rev(gene_order))
  heatmap_plot <- ggplot(heatmap_data, aes(Dataset, Gene, fill = log2FC)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    labs(title = "Top overlapping genes", subtitle = "Ranked by combined absolute change",
         x = NULL, y = NULL, fill = "log2FC") +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 6))
  heatmap_height <- max(7, min(24, length(gene_order) * 0.18))
  ggsave(file.path(figure_dir, "overlap_heatmap_top100.pdf"), heatmap_plot,
         width = 5.5, height = heatmap_height, limitsize = FALSE)
  ggsave(file.path(figure_dir, "overlap_heatmap_top100.png"), heatmap_plot,
         width = 5.5, height = heatmap_height, dpi = 300, limitsize = FALSE)

  make_venn_plot(rna_sig$Gene, strongest$Gene,
                 file.path(figure_dir, "significant_gene_overlap_venn"),
                 if (promoter_only) "RNA-seq and promoter H4K16ac overlap" else "RNA-seq and H4K16ac overlap")

  list(counts = mutate(counts, Analysis = analysis_name),
       correlations = mutate(correlations, Analysis = analysis_name),
       n_chip_genes = nrow(strongest), n_overlap = nrow(final_table))
}

all_result <- run_integration(chip_sig, "all_peaks", promoter_only = FALSE)
promoter_peaks <- chip_sig %>% filter(grepl("promoter", Peak_annotation, ignore.case = TRUE))
promoter_result <- run_integration(promoter_peaks, "promoter_peaks", promoter_only = TRUE)

summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(bind_rows(all_result$counts, promoter_result$counts) %>%
            select(Analysis, everything()), file.path(summary_dir, "all_category_counts.csv"))
write_csv(bind_rows(all_result$correlations, promoter_result$correlations) %>%
            select(Analysis, everything()), file.path(summary_dir, "all_correlation_statistics.csv"))
write_csv(tibble(
  Analysis = c("all_peaks", "promoter_peaks"),
  Significant_RNA_genes = nrow(rna_sig),
  Significant_H4K16ac_genes = c(all_result$n_chip_genes, promoter_result$n_chip_genes),
  Overlapping_genes = c(all_result$n_overlap, promoter_result$n_overlap)
), file.path(summary_dir, "integration_overview.csv"))

writeLines(c(
  "RNA-seq and H4K16ac integration",
  paste("RNA input:", normalizePath(rna_file)),
  paste("H4K16ac input:", normalizePath(chip_file)),
  paste("RNA threshold: padj <", rna_padj_cutoff),
  paste("H4K16ac threshold: FDR <", chip_fdr_cutoff),
  "Duplicate peak rule: lowest FDR per gene; largest absolute Fold breaks ties.",
  "Positive Fold/log2FC = gain/up after p53 knockdown; negative = loss/down.",
  "Promoter analysis includes annotations containing the word 'Promoter'."
), file.path(summary_dir, "analysis_notes.txt"))

message("Integration complete: ", output_dir)
