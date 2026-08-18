# H4K16ac/RNA hypergeometric ORA

This workflow tests each integrated gene list separately for over-representation
of GO Biological Process, KEGG, and Reactome terms with `clusterProfiler` and
`ReactomePA`. It maps human symbols to Entrez IDs, restricts every foreground to
the supplied background universe, and applies Benjamini-Hochberg correction.

The background file must be the complete set of genes that could have entered
the integrated analysis—not the union of significant genes. Each CSV/TSV input
needs a symbol column named `Gene`, `SYMBOL`, `gene_name`, `gene_symbol`, or
`Symbol`.

The preferred entry point is Python. It validates all inputs and then invokes R
for the clusterProfiler-specific statistical work:

```bash
python \
  /gpfs/projects/b1042/LauberthLab/BenFolder/shared/scripts/over_representation_analysis/run_h4k16ac_rna_ora.py \
  --background /path/to/background_gene_universe.csv \
  --loss-down /path/to/loss_down.csv \
  --gain-up /path/to/gain_up.csv \
  --loss-up /path/to/loss_up.csv \
  --gain-down /path/to/gain_down.csv \
  --outdir /path/to/ora_results
```

The underlying R workflow can also be called directly when desired.

## Run both scopes with matched universes

For this dataset, use the two-scope runner. It constructs each universe from
complete, nonsignificant-inclusive source results and runs each foreground only
against its matching universe:

```bash
python \
  /gpfs/projects/b1042/LauberthLab/BenFolder/shared/scripts/over_representation_analysis/run_both_matched_ora.py
```

The all-peak universe is the intersection of RNA-eligible genes and genes linked
to any tested DiffBind peak. The promoter universe is independently constructed
as the intersection of RNA-eligible genes and genes linked to a tested peak whose
annotation contains `promoter`, exactly matching the upstream integration rule.
Genes with missing RNA `log2FoldChange`/`padj` or ChIP `Fold`/`FDR` are not
eligible. The exact universes and an audit table are saved under `universes/`.

Optional controls are `--fdr 0.05`, `--show 15`, `--min-gs-size 10`, and
`--max-gs-size 500`. KEGG uses the current online KEGG annotation; GO and
Reactome use the versions recorded in `sessionInfo.txt`.

For every gene list/database pair, the script writes an `all_pathways.csv` even
when the foreground is empty, no pathways pass FDR, or a database call fails.
The `Status` column in `tables/ora_run_summary.csv` distinguishes these cases.
All tested pathways are retained because enrichment calls use p- and q-value
cutoffs of 1; significance filtering is only applied to plots and summary files.

Outputs include:

- 12 complete pathway CSV files with ratios, counts, raw P values, BH FDR, and
  contributing genes
- mapping details, unmapped symbols, and foreground/background mapping counts
- per-analysis status and significant-pathway counts
- combined significant-pathway summary CSV and formatted PDF/PNG table
- dot and bar plots in PDF and 300-dpi PNG, including explicit placeholders for
  empty or nonsignificant analyses
- `sessionInfo.txt` for package-version provenance
