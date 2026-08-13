#!/usr/bin/env python3
                      
"""Compare YEATS2 and p53 H4K16ac differential peaks and nearby genes."""

import argparse
import csv
import os
import random
from collections import defaultdict


root = os.path.dirname(os.path.abspath(__file__))


def read_peaks(path, direction):
    peaks = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            peaks.append(
                {
                    "chrom": row["seqnames"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "direction": direction,
                    "fold": row.get("Fold", ""),
                    "fdr": row.get("FDR", ""),
                }
            )
    return peaks


def read_annotations(path):
    annotations = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["seqnames"], int(row["start"]), int(row["end"]))
            annotations[key] = {
                "ensembl": row.get("ENSEMBL", ""),
                "symbol": row.get("SYMBOL", ""),
                "gene_name": row.get("GENENAME", ""),
                "annotation": row.get("annotation", ""),
                "distance_to_tss": row.get("distanceToTSS", ""),
            }
    return annotations


def peak_key(peak):
    return peak["chrom"], peak["start"], peak["end"]


def index_by_chrom(peaks):
    result = defaultdict(list)
    for peak in peaks:
        result[peak["chrom"]].append(peak)
    for chrom in result:
        result[chrom].sort(key=lambda p: (p["start"], p["end"]))
    return result


def find_overlaps(query, target_index):
    matches = []
    for target in target_index.get(query["chrom"], []):
        if target["start"] > query["end"]:
            break
        if target["end"] >= query["start"]:
            matches.append(target)
    return matches


def strip_version(gene_id):
    if not gene_id or gene_id.strip() in {"NA", "NaN", "nan", ".", "None"}:
        return ""
    return gene_id.split(".", 1)[0]


def annotations_match(left, right):
    left_id = strip_version(left.get("ensembl", ""))
    right_id = strip_version(right.get("ensembl", ""))
    if left_id and right_id:
        return left_id == right_id
    left_symbol = left.get("symbol", "")
    right_symbol = right.get("symbol", "")
    return bool(left_symbol and right_symbol and left_symbol == right_symbol)


def annotation_for(peak, annotation_map):
    return annotation_map.get(peak_key(peak), {})


def gene_set(peaks, annotation_map):
    genes = set()
    for peak in peaks:
        annotation = annotation_for(peak, annotation_map)
        gene = strip_version(annotation.get("ensembl", ""))
        if gene:
            genes.add(gene)
    return genes


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def metric_flags(peak, p53_index):
    directions = {match["direction"] for match in find_overlaps(peak, p53_index)}
    own = peak["direction"]
    return {
        "any_overlap": bool(directions),
        "gain_gain": own == "gain" and "gain" in directions,
        "gain_loss": own == "gain" and "loss" in directions,
        "loss_gain": own == "loss" and "gain" in directions,
        "loss_loss": own == "loss" and "loss" in directions,
        "concordant": (own == "gain" and "gain" in directions)
        or (own == "loss" and "loss" in directions),
        "discordant": (own == "gain" and "loss" in directions)
        or (own == "loss" and "gain" in directions),
    }


def permutation_test(background, n_gain, n_loss, p53_index, observed, iterations, seed):
    background_flags = []
    for peak in background:
        directions = {match["direction"] for match in find_overlaps(peak, p53_index)}
        background_flags.append(("gain" in directions, "loss" in directions))

    metrics = list(observed)
    totals = {metric: 0.0 for metric in metrics}
    exceed = {metric: 0 for metric in metrics}
    rng = random.Random(seed)
    sample_size = n_gain + n_loss
    if sample_size > len(background):
        raise SystemExit("Error: more differential YEATS2 peaks than tested background peaks.")

    for _ in range(iterations):
        selected = rng.sample(range(len(background)), sample_size)
        counts = {metric: 0 for metric in metrics}
        for position, index in enumerate(selected):
            has_gain, has_loss = background_flags[index]
            direction = "gain" if position < n_gain else "loss"
            any_overlap = has_gain or has_loss
            concordant = (direction == "gain" and has_gain) or (direction == "loss" and has_loss)
            discordant = (direction == "gain" and has_loss) or (direction == "loss" and has_gain)
            values = {
                "any_overlap": any_overlap,
                "gain_gain": direction == "gain" and has_gain,
                "gain_loss": direction == "gain" and has_loss,
                "loss_gain": direction == "loss" and has_gain,
                "loss_loss": direction == "loss" and has_loss,
                "concordant": concordant,
                "discordant": discordant,
            }
            for metric in metrics:
                counts[metric] += int(values[metric])
        for metric in metrics:
            totals[metric] += counts[metric]
            exceed[metric] += int(counts[metric] >= observed[metric])

    rows = []
    for metric in metrics:
        expected = totals[metric] / iterations
        enrichment = observed[metric] / expected if expected else "Inf" if observed[metric] else "NA"
        rows.append(
            {
                "metric": metric,
                "observed": observed[metric],
                "random_mean": "{:.4f}".format(expected),
                "fold_enrichment": "{:.4f}".format(enrichment) if isinstance(enrichment, float) else enrichment,
                "empirical_p_value": "{:.6g}".format((exceed[metric] + 1) / (iterations + 1)),
                "permutations": iterations,
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260811)
    parser.add_argument("--output-dir", default=os.path.join(root, "yeats2_p53_overlap"))
    args = parser.parse_args()
    if args.permutations < 1:
        raise SystemExit("Error: --permutations must be at least 1.")

    p53_dir = os.path.join(root, "diffbind_results")
    yeats_dir = os.path.join(root, "diffbind_results_yeats2")
    p53_gain = read_peaks(os.path.join(p53_dir, "H4K16ac_gained_peaks.csv"), "gain")
    p53_loss = read_peaks(os.path.join(p53_dir, "H4K16ac_lost_peaks.csv"), "loss")
    yeats_gain = read_peaks(os.path.join(yeats_dir, "YEATS2_H4K16ac_gained_peaks.csv"), "gain")
    yeats_loss = read_peaks(os.path.join(yeats_dir, "YEATS2_H4K16ac_lost_peaks.csv"), "loss")
    p53_peaks = p53_gain + p53_loss
    yeats_peaks = yeats_gain + yeats_loss
    p53_index = index_by_chrom(p53_peaks)

    p53_annotations = {}
    p53_annotations.update(read_annotations(os.path.join(p53_dir, "H4K16ac_gained_annotated.csv")))
    p53_annotations.update(read_annotations(os.path.join(p53_dir, "H4K16ac_lost_annotated.csv")))
    yeats_annotations = {}
    yeats_annotations.update(read_annotations(os.path.join(yeats_dir, "YEATS2_H4K16ac_gained_annotated.csv")))
    yeats_annotations.update(read_annotations(os.path.join(yeats_dir, "YEATS2_H4K16ac_lost_annotated.csv")))

    overlap_rows = []
    unique_overlap_flags = []
    same_gene_yeats = set()
    for yeats_peak in yeats_peaks:
        matches = find_overlaps(yeats_peak, p53_index)
        unique_overlap_flags.append(metric_flags(yeats_peak, p53_index))
        for p53_peak in matches:
            y_annotation = annotation_for(yeats_peak, yeats_annotations)
            p_annotation = annotation_for(p53_peak, p53_annotations)
            same_gene = annotations_match(y_annotation, p_annotation)
            if same_gene:
                same_gene_yeats.add(peak_key(yeats_peak))
            overlap_bp = min(yeats_peak["end"], p53_peak["end"]) - max(yeats_peak["start"], p53_peak["start"]) + 1
            overlap_rows.append(
                {
                    "chrom": yeats_peak["chrom"],
                    "yeats2_start": yeats_peak["start"],
                    "yeats2_end": yeats_peak["end"],
                    "yeats2_direction": yeats_peak["direction"],
                    "yeats2_fold": yeats_peak["fold"],
                    "yeats2_fdr": yeats_peak["fdr"],
                    "p53_start": p53_peak["start"],
                    "p53_end": p53_peak["end"],
                    "p53_direction": p53_peak["direction"],
                    "p53_fold": p53_peak["fold"],
                    "p53_fdr": p53_peak["fdr"],
                    "overlap_bp": overlap_bp,
                    "direction_match": yeats_peak["direction"] == p53_peak["direction"],
                    "yeats2_ensembl": y_annotation.get("ensembl", ""),
                    "yeats2_symbol": y_annotation.get("symbol", ""),
                    "p53_ensembl": p_annotation.get("ensembl", ""),
                    "p53_symbol": p_annotation.get("symbol", ""),
                    "same_nearby_gene": same_gene,
                }
            )

    metric_names = ["any_overlap", "gain_gain", "gain_loss", "loss_gain", "loss_loss", "concordant", "discordant"]
    observed = {metric: sum(int(flags[metric]) for flags in unique_overlap_flags) for metric in metric_names}

    background = read_peaks(os.path.join(yeats_dir, "DiffBind_YEATS2_all_peaks.csv"), "background")
    deduplicated = {}
    for peak in background:
        deduplicated[peak_key(peak)] = peak
    background = list(deduplicated.values())
    permutation_rows = permutation_test(
        background, len(yeats_gain), len(yeats_loss), p53_index, observed, args.permutations, args.seed
    )

    p53_gene_sets = {"gain": gene_set(p53_gain, p53_annotations), "loss": gene_set(p53_loss, p53_annotations)}
    yeats_gene_sets = {"gain": gene_set(yeats_gain, yeats_annotations), "loss": gene_set(yeats_loss, yeats_annotations)}
    gene_symbols = {}
    for annotation in list(p53_annotations.values()) + list(yeats_annotations.values()):
        gene = strip_version(annotation.get("ensembl", ""))
        symbol = annotation.get("symbol", "")
        if gene and symbol not in {"", "NA", "."}:
            gene_symbols[gene] = symbol
    shared_gene_rows = []
    shared_gene_counts = {}
    for yeats_direction in ("gain", "loss"):
        for p53_direction in ("gain", "loss"):
            shared = sorted(yeats_gene_sets[yeats_direction] & p53_gene_sets[p53_direction])
            shared_gene_counts[(yeats_direction, p53_direction)] = len(shared)
            for gene in shared:
                shared_gene_rows.append(
                    {
                        "ensembl": gene,
                        "symbol": gene_symbols.get(gene, ""),
                        "yeats2_direction": yeats_direction,
                        "p53_direction": p53_direction,
                        "direction_match": yeats_direction == p53_direction,
                    }
                )

    os.makedirs(args.output_dir, exist_ok=True)
    write_tsv(
        os.path.join(args.output_dir, "overlapping_peak_pairs.tsv"),
        [
            "chrom", "yeats2_start", "yeats2_end", "yeats2_direction", "yeats2_fold", "yeats2_fdr",
            "p53_start", "p53_end", "p53_direction", "p53_fold", "p53_fdr", "overlap_bp",
            "direction_match", "yeats2_ensembl", "yeats2_symbol", "p53_ensembl", "p53_symbol",
            "same_nearby_gene",
        ],
        overlap_rows,
    )
    write_tsv(
        os.path.join(args.output_dir, "permutation_enrichment.tsv"),
        ["metric", "observed", "random_mean", "fold_enrichment", "empirical_p_value", "permutations"],
        permutation_rows,
    )
    write_tsv(
        os.path.join(args.output_dir, "shared_nearby_genes.tsv"),
        ["ensembl", "symbol", "yeats2_direction", "p53_direction", "direction_match"],
        shared_gene_rows,
    )

    annotated_pairs = sum(
        bool(row["yeats2_ensembl"] or row["yeats2_symbol"])
        and bool(row["p53_ensembl"] or row["p53_symbol"])
        for row in overlap_rows
    )
    same_gene_pairs = sum(row["same_nearby_gene"] for row in overlap_rows)
    summary_rows = [
        ("p53_gain_peaks", len(p53_gain)),
        ("p53_loss_peaks", len(p53_loss)),
        ("yeats2_gain_peaks", len(yeats_gain)),
        ("yeats2_loss_peaks", len(yeats_loss)),
        ("yeats2_peaks_overlapping_any_p53_diff_peak", observed["any_overlap"]),
        ("yeats2_overlap_percent", "{:.2f}".format(100 * observed["any_overlap"] / len(yeats_peaks))),
        ("gain_gain_unique_yeats2_peaks", observed["gain_gain"]),
        ("gain_loss_unique_yeats2_peaks", observed["gain_loss"]),
        ("loss_gain_unique_yeats2_peaks", observed["loss_gain"]),
        ("loss_loss_unique_yeats2_peaks", observed["loss_loss"]),
        ("concordant_unique_yeats2_peaks", observed["concordant"]),
        ("discordant_unique_yeats2_peaks", observed["discordant"]),
        ("overlapping_peak_pairs", len(overlap_rows)),
        ("overlap_pairs_with_both_gene_annotations", annotated_pairs),
        ("overlap_pairs_with_same_nearby_gene", same_gene_pairs),
        ("unique_overlapping_yeats2_peaks_with_same_nearby_gene", len(same_gene_yeats)),
        ("shared_nearby_genes_gain_gain", shared_gene_counts[("gain", "gain")]),
        ("shared_nearby_genes_gain_loss", shared_gene_counts[("gain", "loss")]),
        ("shared_nearby_genes_loss_gain", shared_gene_counts[("loss", "gain")]),
        ("shared_nearby_genes_loss_loss", shared_gene_counts[("loss", "loss")]),
        ("shared_nearby_gene_direction_combinations", len(shared_gene_rows)),
    ]
    write_tsv(
        os.path.join(args.output_dir, "summary.tsv"),
        ["metric", "value"],
        ({"metric": metric, "value": value} for metric, value in summary_rows),
    )
    print("Wrote YEATS2/p53 overlap results to {}".format(args.output_dir))


if __name__ == "__main__":
    main()
