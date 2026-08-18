#!/usr/bin/env Rscript

# Hypergeometric over-representation analysis (ORA) for integrated H4K16ac/RNA
# gene sets. All foregrounds and the user-defined universe are mapped from human
# gene symbols to Entrez IDs before testing. GO BP, KEGG, and Reactome are tested
# independently with Benjamini-Hochberg correction.

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(org.Hs.eg.db)
  library(ReactomePA)
})

options(stringsAsFactors = FALSE)

usage <- paste(
  "Usage:",
  "  Rscript h4k16ac_rna_ora.R --background BACKGROUND.csv \\",
  "    --loss-down LOSS_DOWN.csv --gain-up GAIN_UP.csv \\",
  "    --loss-up LOSS_UP.csv --gain-down GAIN_DOWN.csv \\",
  "    --outdir OUTPUT_DIR [--fdr 0.05] [--show 15] \\",
  "    [--min-gs-size 10] [--max-gs-size 500]",
  "",
  "Input tables must contain a gene-symbol column named one of:",
  "Gene, SYMBOL, gene_name, gene_symbol, or Symbol.",
  sep = "\n"
)

parse_args <- function(x) {
  if (length(x) == 0L || any(x %in% c("-h", "--help"))) {
    cat(usage, "\n")
    quit(save = "no", status = if (length(x) == 0L) 1L else 0L)
  }
  if (length(x) %% 2L != 0L || any(!grepl("^--", x[seq(1L, length(x), 2L)]))) {
    stop("Arguments must be supplied as --name value pairs.\n\n", usage)
  }
  keys <- sub("^--", "", x[seq(1L, length(x), 2L)])
  values <- x[seq(2L, length(x), 2L)]
  stats::setNames(as.list(values), keys)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("background", "loss-down", "gain-up", "loss-up", "gain-down", "outdir")
missing_args <- setdiff(required, names(args))
if (length(missing_args) > 0L) {
  stop("Missing required argument(s): ", paste(paste0("--", missing_args), collapse = ", "),
       "\n\n", usage)
}

outdir <- normalizePath(args$outdir, mustWork = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

# Move the operator definition ahead of its first runtime use after parsing.
fdr_cutoff <- as.numeric(args$fdr %||% "0.05")
n_show <- as.integer(args$show %||% "15")
min_gs_size <- as.integer(args$`min-gs-size` %||% "10")
max_gs_size <- as.integer(args$`max-gs-size` %||% "500")

if (!is.finite(fdr_cutoff) || fdr_cutoff <= 0 || fdr_cutoff > 1) stop("--fdr must be in (0, 1].")
if (is.na(n_show) || n_show < 1L) stop("--show must be a positive integer.")
if (is.na(min_gs_size) || is.na(max_gs_size) || min_gs_size < 1L || max_gs_size < min_gs_size) {
  stop("Gene-set size limits are invalid.")
}

input_files <- c(
  loss_down = args$`loss-down`,
  gain_up = args$`gain-up`,
  loss_up = args$`loss-up`,
  gain_down = args$`gain-down`
)
display_names <- c(
  loss_down = "H4K16ac loss + RNA down",
  gain_up = "H4K16ac gain + RNA up",
  loss_up = "H4K16ac loss + RNA up",
  gain_down = "H4K16ac gain + RNA down"
)
all_input_files <- c(background = args$background, input_files)
not_found <- names(all_input_files)[!file.exists(all_input_files)]
if (length(not_found) > 0L) {
  stop("Input file(s) not found: ",
       paste(paste0(not_found, "=", all_input_files[not_found]), collapse = ", "))
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(outdir, "tables")
figure_dir <- file.path(outdir, "figures")
mapping_dir <- file.path(outdir, "mapping")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mapping_dir, recursive = TRUE, showWarnings = FALSE)
# Headless compute nodes often lack a writable fontconfig cache.
font_cache <- file.path(tempdir(), "fontconfig-cache")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache)

message("H4K16ac/RNA pathway ORA")
message("Output: ", outdir)
message("FDR threshold: ", fdr_cutoff)

read_gene_symbols <- function(path) {
  first_line <- readLines(path, n = 1L, warn = FALSE)
  sep <- if (grepl("\\t", first_line)) "\t" else ","
  x <- read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "",
                  check.names = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  candidates <- c("Gene", "SYMBOL", "gene_name", "gene_symbol", "Symbol")
  hit <- names(x)[tolower(names(x)) %in% tolower(candidates)]
  if (length(hit) == 0L) {
    stop("No recognized gene-symbol column in ", path, ". Found: ", paste(names(x), collapse = ", "))
  }
  symbols <- trimws(as.character(x[[hit[[1L]]]]))
  unique(symbols[!is.na(symbols) & nzchar(symbols)])
}

map_symbols <- function(symbols, set_name) {
  symbols <- unique(trimws(symbols))
  symbols <- symbols[nzchar(symbols)]
  if (length(symbols) == 0L) {
    return(data.frame(Set = character(), SYMBOL = character(), ENTREZID = character()))
  }
  mapped <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db, keys = symbols, keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID")
  )) %>%
    filter(!is.na(ENTREZID), nzchar(ENTREZID)) %>%
    distinct(SYMBOL, ENTREZID) %>%
    mutate(Set = set_name, .before = 1L)
  mapped
}

background_symbols <- read_gene_symbols(args$background)
foreground_symbols <- lapply(input_files, read_gene_symbols)
background_map <- map_symbols(background_symbols, "background")
foreground_maps <- Map(map_symbols, foreground_symbols, names(foreground_symbols))

background_entrez <- unique(background_map$ENTREZID)
if (length(background_entrez) == 0L) stop("No background symbols mapped to Entrez IDs; ORA cannot run.")

# Only genes in the declared universe may enter a hypergeometric foreground.
foreground_entrez <- lapply(foreground_maps, function(x) intersect(unique(x$ENTREZID), background_entrez))

mapping_detail <- bind_rows(background_map, bind_rows(foreground_maps))
write.csv(mapping_detail, file.path(mapping_dir, "symbol_to_entrez_mapping.csv"), row.names = FALSE, na = "")

unmapped_detail <- bind_rows(
  data.frame(
    Set = rep("background", length(setdiff(background_symbols, background_map$SYMBOL))),
    SYMBOL = setdiff(background_symbols, background_map$SYMBOL)
  ),
  bind_rows(lapply(names(foreground_symbols), function(nm) {
    missing_symbols <- setdiff(foreground_symbols[[nm]], foreground_maps[[nm]]$SYMBOL)
    data.frame(Set = rep(nm, length(missing_symbols)), SYMBOL = missing_symbols)
  }))
)
write.csv(unmapped_detail, file.path(mapping_dir, "unmapped_symbols.csv"), row.names = FALSE, na = "")

outside_universe <- bind_rows(lapply(names(foreground_maps), function(nm) {
  x <- foreground_maps[[nm]] %>% filter(!ENTREZID %in% background_entrez)
  if (nrow(x) == 0L) return(x)
  x$Set <- nm
  x
}))
write.csv(outside_universe,
          file.path(mapping_dir, "mapped_foreground_genes_outside_universe.csv"),
          row.names = FALSE, na = "")

mapping_summary <- bind_rows(
  data.frame(
    Gene_set = "background", Display_name = "Background universe",
    Input_symbols = length(background_symbols), Mapped_symbols = length(unique(background_map$SYMBOL)),
    Unique_Entrez_IDs = length(background_entrez), Entrez_IDs_in_universe = length(background_entrez)
  ),
  bind_rows(lapply(names(foreground_symbols), function(nm) {
    data.frame(
      Gene_set = nm, Display_name = unname(display_names[[nm]]),
      Input_symbols = length(foreground_symbols[[nm]]),
      Mapped_symbols = length(unique(foreground_maps[[nm]]$SYMBOL)),
      Unique_Entrez_IDs = length(unique(foreground_maps[[nm]]$ENTREZID)),
      Entrez_IDs_in_universe = length(foreground_entrez[[nm]])
    )
  }))
)
write.csv(mapping_summary, file.path(table_dir, "gene_mapping_summary.csv"), row.names = FALSE)

empty_result <- function() {
  data.frame(
    Pathway_ID = character(), Pathway = character(), GeneRatio = character(), BgRatio = character(),
    Count = integer(), P_value = numeric(), Adjusted_P_value_FDR = numeric(),
    Genes = character(), stringsAsFactors = FALSE
  )
}

format_result <- function(enrichment, database) {
  if (is.null(enrichment)) return(empty_result())
  x <- as.data.frame(enrichment)
  if (nrow(x) == 0L) return(empty_result())
  if (database != "GO_BP") {
    enrichment <- setReadable(enrichment, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    x <- as.data.frame(enrichment)
  }
  x %>%
    transmute(
      Pathway_ID = ID,
      Pathway = Description,
      GeneRatio = GeneRatio,
      BgRatio = BgRatio,
      Count = Count,
      P_value = pvalue,
      Adjusted_P_value_FDR = p.adjust,
      Genes = gsub("/", ";", geneID, fixed = TRUE)
    ) %>%
    arrange(Adjusted_P_value_FDR, P_value)
}

run_one_database <- function(genes, database) {
  if (length(genes) == 0L) return(list(result = empty_result(), status = "empty foreground"))
  tryCatch({
    enrichment <- switch(
      database,
      GO_BP = enrichGO(
        gene = genes, universe = background_entrez, OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = min_gs_size,
        maxGSSize = max_gs_size, readable = TRUE
      ),
      KEGG = enrichKEGG(
        gene = genes, universe = background_entrez, organism = "hsa",
        keyType = "ncbi-geneid", pAdjustMethod = "BH", pvalueCutoff = 1,
        qvalueCutoff = 1, minGSSize = min_gs_size, maxGSSize = max_gs_size,
        use_internal_data = FALSE
      ),
      Reactome = enrichPathway(
        gene = genes, universe = background_entrez, organism = "human",
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = min_gs_size, maxGSSize = max_gs_size, readable = FALSE
      )
    )
    list(result = format_result(enrichment, database), status = "completed")
  }, error = function(e) {
    warning(database, " failed: ", conditionMessage(e), call. = FALSE)
    list(result = empty_result(), status = paste0("ERROR: ", conditionMessage(e)))
  })
}

ratio_to_numeric <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2L) return(NA_real_)
    as.numeric(z[[1L]]) / as.numeric(z[[2L]])
  }, numeric(1L))
}

plot_placeholder <- function(title, subtitle) {
  ggplot() +
    annotate("text", x = 0, y = 0.08, label = "No significant pathways", size = 6, fontface = "bold") +
    annotate("text", x = 0, y = -0.08, label = subtitle, size = 3.6, color = "grey35") +
    xlim(-1, 1) + ylim(-0.4, 0.4) + labs(title = title) + theme_void(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

make_plots <- function(result, set_id, database, status) {
  sig <- result %>% filter(!is.na(Adjusted_P_value_FDR), Adjusted_P_value_FDR < fdr_cutoff)
  title_suffix <- paste0(display_names[[set_id]], " — ", database)
  if (nrow(sig) == 0L) {
    subtitle <- if (status == "completed") paste0("BH FDR ≥ ", fdr_cutoff) else status
    dot <- plot_placeholder(paste("ORA dot plot:", title_suffix), subtitle)
    bar <- plot_placeholder(paste("ORA bar plot:", title_suffix), subtitle)
  } else {
    top <- sig %>% slice_head(n = n_show) %>%
      mutate(GeneRatio_numeric = ratio_to_numeric(GeneRatio),
             minus_log10_FDR = -log10(pmax(Adjusted_P_value_FDR, .Machine$double.xmin)),
             Pathway_plot = factor(Pathway, levels = rev(Pathway)))
    dot <- ggplot(top, aes(GeneRatio_numeric, Pathway_plot, size = Count,
                           color = Adjusted_P_value_FDR)) +
      geom_point(alpha = 0.9) +
      scale_color_viridis_c(option = "magma", direction = -1, trans = "log10",
                            name = "BH FDR") +
      scale_size_continuous(range = c(3, 9)) +
      labs(title = "Pathway over-representation", subtitle = title_suffix,
           x = "Gene ratio", y = NULL) + theme_publication()
    bar <- ggplot(top, aes(minus_log10_FDR, Pathway_plot, fill = Count)) +
      geom_col(width = 0.72) +
      scale_fill_viridis_c(option = "cividis", name = "Gene count") +
      labs(title = "Pathway over-representation", subtitle = title_suffix,
           x = expression(-log[10]("BH FDR")), y = NULL) + theme_publication()
  }
  prefix <- file.path(figure_dir, paste(set_id, database, sep = "__"))
  ggsave(paste0(prefix, "__dotplot.pdf"), dot, width = 9.5, height = 6.5)
  ggsave(paste0(prefix, "__dotplot.png"), dot, width = 9.5, height = 6.5, dpi = 300)
  ggsave(paste0(prefix, "__barplot.pdf"), bar, width = 9.5, height = 6.5)
  ggsave(paste0(prefix, "__barplot.png"), bar, width = 9.5, height = 6.5, dpi = 300)
}

theme_publication <- function() {
  theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(color = "grey25", margin = margin(b = 8)),
      axis.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 18, 10, 10)
    )
}

databases <- c("GO_BP", "KEGG", "Reactome")
all_results <- list()
run_summary <- list()

for (set_id in names(foreground_entrez)) {
  message("\n", display_names[[set_id]])
  message("  Foreground symbols: ", length(foreground_symbols[[set_id]]))
  message("  Successfully mapped symbols: ", length(unique(foreground_maps[[set_id]]$SYMBOL)))
  message("  Unique Entrez IDs in background: ", length(foreground_entrez[[set_id]]))

  for (database in databases) {
    ans <- run_one_database(foreground_entrez[[set_id]], database)
    result <- ans$result %>%
      mutate(Gene_set = rep(set_id, n()), Database = rep(database, n())) %>%
      select(Gene_set, Database, everything())
    output_csv <- file.path(table_dir, paste0(set_id, "__", database, "__all_pathways.csv"))
    write.csv(result, output_csv, row.names = FALSE, na = "")

    n_sig <- sum(!is.na(result$Adjusted_P_value_FDR) & result$Adjusted_P_value_FDR < fdr_cutoff)
    message("  ", database, ": ", n_sig, " significant pathways (", ans$status, ")")
    run_summary[[paste(set_id, database, sep = "__")]] <- data.frame(
      Gene_set = set_id, Display_name = display_names[[set_id]], Database = database,
      Foreground_symbols = length(foreground_symbols[[set_id]]),
      Mapped_symbols = length(unique(foreground_maps[[set_id]]$SYMBOL)),
      Foreground_Entrez_in_universe = length(foreground_entrez[[set_id]]),
      Background_symbols = length(background_symbols),
      Background_mapped_symbols = length(unique(background_map$SYMBOL)),
      Background_Entrez_IDs = length(background_entrez),
      Tested_pathways = nrow(result), Significant_pathways = n_sig,
      FDR_threshold = fdr_cutoff, Status = ans$status
    )
    all_results[[paste(set_id, database, sep = "__")]] <- result
    make_plots(result, set_id, database, ans$status)
  }
}

run_summary <- bind_rows(run_summary)
write.csv(run_summary, file.path(table_dir, "ora_run_summary.csv"), row.names = FALSE, na = "")

all_results_df <- bind_rows(all_results)
significant_summary <- all_results_df %>%
  filter(!is.na(Adjusted_P_value_FDR), Adjusted_P_value_FDR < fdr_cutoff) %>%
  group_by(Gene_set, Database) %>%
  slice_head(n = n_show) %>%
  ungroup()
write.csv(significant_summary, file.path(table_dir, "significant_pathway_summary.csv"),
          row.names = FALSE, na = "")

# A compact, publication-ready table containing the leading significant terms.
if (nrow(significant_summary) > 0L) {
  table_data <- significant_summary %>%
    mutate(
      `Gene set` = display_names[Gene_set],
      FDR = formatC(Adjusted_P_value_FDR, format = "e", digits = 2),
      `P value` = formatC(P_value, format = "e", digits = 2)
    ) %>%
    select(`Gene set`, Database, Pathway, GeneRatio, BgRatio, Count, `P value`, FDR) %>%
    slice_head(n = 40)
  table_plot <- gridExtra::tableGrob(table_data, rows = NULL,
                                     theme = gridExtra::ttheme_minimal(base_size = 8))
  table_height <- max(4, 0.27 * nrow(table_data) + 1.2)
  ggsave(file.path(figure_dir, "significant_pathway_summary_table.pdf"), table_plot,
         width = 14, height = table_height, limitsize = FALSE)
  ggsave(file.path(figure_dir, "significant_pathway_summary_table.png"), table_plot,
         width = 14, height = table_height, dpi = 300, limitsize = FALSE)
} else {
  table_plot <- plot_placeholder("Significant pathway summary", paste0("No terms at BH FDR < ", fdr_cutoff))
  ggsave(file.path(figure_dir, "significant_pathway_summary_table.pdf"), table_plot, width = 9, height = 4)
  ggsave(file.path(figure_dir, "significant_pathway_summary_table.png"), table_plot,
         width = 9, height = 4, dpi = 300)
}

capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
writeLines(c(
  "Hypergeometric ORA of integrated H4K16ac/RNA gene sets",
  paste0("Completed: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("Background input symbols: ", length(background_symbols)),
  paste0("Background mapped symbols: ", length(unique(background_map$SYMBOL))),
  paste0("Background unique Entrez IDs: ", length(background_entrez)),
  paste0("BH FDR threshold: ", fdr_cutoff),
  paste0("Gene-set size range: ", min_gs_size, "-", max_gs_size),
  "See tables/ora_run_summary.csv for per-list/per-database counts and status."
), file.path(outdir, "README_results.txt"))

message("\nBackground")
message("  Input symbols: ", length(background_symbols))
message("  Successfully mapped symbols: ", length(unique(background_map$SYMBOL)))
message("  Unique Entrez IDs: ", length(background_entrez))
message("\nORA complete: ", outdir)
