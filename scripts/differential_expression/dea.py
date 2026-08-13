#!/usr/bin/env python3
                      
"""Differential expression analysis for the RNA-seq featureCounts matrix."""

from pathlib import Path
import argparse
import pickle
import sys

import pandas as pd


project_dir = Path(__file__).resolve().parent
default_counts = project_dir / "rna_seq_featurecounts" / "rna_seq_featureCounts_cleaned.txt"
default_metadata = project_dir / "rna_seq_metadata.tsv"
default_output = project_dir / "rna_seq_dea"


def write_metadata_template(counts_path, metadata_path):
    counts = pd.read_csv(counts_path, sep="\t", nrows=0)
    samples = [column for column in counts.columns if column not in {"gene_name", "Geneid", "Length"}]
    metadata = pd.DataFrame({"sample": samples, "condition": ["TODO"] * len(samples)})
    metadata.to_csv(metadata_path, sep="\t", index=False)


def load_counts_and_metadata(counts_path, metadata_path, min_cpm=0.3, min_samples=2):
    if not counts_path.exists():
        raise SystemExit("Count matrix not found: {}".format(counts_path))
    if not metadata_path.exists():
        write_metadata_template(counts_path, metadata_path)
        raise SystemExit(
            "Metadata template created at {}. Fill in the condition column and rerun.".format(metadata_path)
        )

    table = pd.read_csv(counts_path, sep="\t")
    gene_column = "gene_name" if "gene_name" in table.columns else "Geneid"
    if gene_column not in table.columns:
        raise SystemExit("Count matrix must contain a gene_name or Geneid column.")
    if "Length" in table.columns:
        table = table.drop(columns="Length")

    sample_columns = [column for column in table.columns if column != gene_column]
    if not sample_columns:
        raise SystemExit("No sample columns found in {}.".format(counts_path))
    table[sample_columns] = table[sample_columns].apply(pd.to_numeric, errors="raise")
    if (table[sample_columns] < 0).any().any():
        raise SystemExit("Counts must be nonnegative integers.")

                                                                              
                                                   
    counts = table.groupby(gene_column, sort=False)[sample_columns].sum()
    counts.index = counts.index.astype(str)
    counts = counts.loc[~counts.index.isin({"", "NA", "nan", "None"})]

    library_sizes = counts.sum(axis=0)
    cpm = counts.div(library_sizes.replace(0, 1), axis=1) * 1_000_000
    keep = (cpm > min_cpm).sum(axis=1) >= min_samples
    counts = counts.loc[keep]
    if counts.empty:
        raise SystemExit("No genes passed the CPM filter.")

    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str).fillna("")
    required = {"sample", "condition"}
    missing = required - set(metadata.columns)
    if missing:
        raise SystemExit("Metadata is missing required column(s): {}.".format(", ".join(sorted(missing))))
    metadata = metadata[["sample", "condition"]]
    if metadata["sample"].duplicated().any():
        raise SystemExit("Metadata contains duplicate sample names.")
    if metadata["condition"].str.strip().isin({"", "TODO", "NA"}).any():
        raise SystemExit("Every sample must have a real condition in {}.".format(metadata_path))

    metadata = metadata.set_index("sample")
    missing_samples = sorted(set(sample_columns) - set(metadata.index))
    extra_samples = sorted(set(metadata.index) - set(sample_columns))
    if missing_samples or extra_samples:
        message = []
        if missing_samples:
            message.append("missing metadata for: {}".format(", ".join(missing_samples)))
        if extra_samples:
            message.append("metadata has samples absent from counts: {}".format(", ".join(extra_samples)))
        raise SystemExit("Sample mismatch: " + "; ".join(message))

    metadata = metadata.loc[sample_columns]
    metadata.index.name = "sample"
    return counts.T.astype(int), metadata


def run_dea(counts_path, metadata_path, output_dir, case=None, control=None, min_cpm=0.3, min_samples=2):
    counts, metadata = load_counts_and_metadata(counts_path, metadata_path, min_cpm, min_samples)
    conditions = list(metadata["condition"].unique())
    if case is None or control is None:
        if len(conditions) != 2:
            raise SystemExit(
                "Specify --case and --control; metadata contains conditions: {}".format(", ".join(conditions))
            )
        control, case = conditions[0], conditions[1]
    if case not in conditions or control not in conditions:
        raise SystemExit("Both --case and --control must occur in the metadata condition column.")
    if case == control:
        raise SystemExit("--case and --control must be different conditions.")

    return run_deseq2(counts, metadata, output_dir, case, control)


def run_exploratory(counts_path, metadata_path, output_dir, case, control, min_cpm=0.3, min_samples=2, min_log2fc=1.0):
    counts, metadata = load_counts_and_metadata(counts_path, metadata_path, min_cpm, min_samples)
    conditions = set(metadata["condition"])
    if case not in conditions or control not in conditions:
        raise SystemExit("Both --case and --control must occur in the metadata condition column.")

    case_samples = metadata.index[metadata["condition"] == case]
    control_samples = metadata.index[metadata["condition"] == control]
    if len(case_samples) == 0 or len(control_samples) == 0:
        raise SystemExit("Both exploratory comparison groups must contain samples.")

                                                                              
                                                                         
                                                   
    library_sizes = counts.sum(axis=1)
    cpm = counts.div(library_sizes.replace(0, 1), axis=0) * 1_000_000
    case_mean = cpm.loc[case_samples].mean(axis=0)
    control_mean = cpm.loc[control_samples].mean(axis=0)
    results = pd.DataFrame(
        {
            "baseMean_CPM": (case_mean + control_mean) / 2,
            "case_mean_CPM": case_mean,
            "control_mean_CPM": control_mean,
            "log2FoldChange": (case_mean + 1).div(control_mean + 1).map(lambda value: __import__("math").log2(value)),
            "pvalue": float("nan"),
            "padj": float("nan"),
        }
    )
    results.index.name = "gene_name"
    results = results.sort_values("log2FoldChange", ascending=False)
    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "exploratory_results.csv")
    selected = results[results["log2FoldChange"].abs() >= min_log2fc]
    selected.to_csv(output_dir / "exploratory_fold_change_genes.csv")
    metadata.to_csv(output_dir / "metadata.csv")
    print("Exploratory fold-change results written to {}".format(output_dir))
    print("No p-values or adjusted p-values were calculated.")
    return results, output_dir


def run_deseq2(counts, metadata, output_dir, case, control):

    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.default_inference import DefaultInference
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:
        raise SystemExit(
            "PyDESeq2 is unavailable. Run this script in the pydeseq2 conda environment."
        ) from exc

    print("Running DESeq2 on {} samples and {} genes.".format(counts.shape[0], counts.shape[1]))
    inference = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design="~condition",
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", case, control], inference=inference)
    stats.summary()
    results = stats.results_df.sort_values("padj", na_position="last")

    output_dir.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(output_dir / "metadata.csv")
    results.to_csv(output_dir / "results.csv")
    with open(output_dir / "result_adata.pkl", "wb") as handle:
        pickle.dump(dds.to_picklable_anndata(), handle)
    print("Results written to {}".format(output_dir / "results.csv"))
    return results, output_dir


def write_significant(results, output_dir, padj_threshold=0.05, min_log2fc=1.0):
    significant = results.loc[
        (results["padj"] <= padj_threshold)
        & (results["log2FoldChange"].abs() >= min_log2fc)
    ].sort_values("padj")
    significant.to_csv(output_dir / "significant_results.csv")
    print("Significant genes written to {} ({})".format(output_dir / "significant_results.csv", len(significant)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, default=default_counts)
    parser.add_argument("--metadata", type=Path, default=default_metadata)
    parser.add_argument("--output-dir", type=Path, default=default_output)
    parser.add_argument("--case")
    parser.add_argument("--control")
    parser.add_argument("--min-cpm", type=float, default=0.3)
    parser.add_argument("--min-samples", type=int, default=2)
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--min-log2fc", type=float, default=1.0)
    parser.add_argument("--exploratory", action="store_true", help="Compatibility flag; exploratory CPM/log2 fold-change analysis is the default.")
    parser.add_argument("--deseq2", action="store_true", help="Run statistical DESeq2 instead; this requires biological replicates.")
    args = parser.parse_args()
    if not args.case or not args.control:
        raise SystemExit("--case and --control are required.")
    if args.deseq2:
        results, output_dir = run_dea(
            args.counts, args.metadata, args.output_dir, args.case, args.control, args.min_cpm, args.min_samples
        )
        write_significant(results, output_dir, args.padj, args.min_log2fc)
    else:
        run_exploratory(
            args.counts, args.metadata, args.output_dir, args.case, args.control,
            args.min_cpm, args.min_samples, args.min_log2fc,
        )


if __name__ == "__main__":
    main()


"""to run


python3 dea.py --exploratory \
  --case shp53_0hr \
  --control shLacZ_0hr \
  --output-dir rna_seq_dea/my_comparison"""
