#!/usr/bin/env python3
"""Align the shP53 single-end RNA-seq FASTQs with HISAT2."""

import argparse
import os
import re
import shutil
import subprocess


PROJECT = "/gpfs/projects/b1042/LauberthLab/BenFolder"
FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def validate_hisat2_index(prefix):
    suffixes = (".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8")
    for extension in (".ht2", ".ht2l"):
        if all(os.path.isfile(prefix + suffix + extension) for suffix in suffixes):
            return
    raise SystemExit(f"ERROR: incomplete HISAT2 index: {prefix}")


def sample_name(path):
    """Use the SRR accession, even when the FASTQ has a descriptive filename."""
    match = re.search(r"(SRR\d+)", os.path.basename(path))
    if match:
        return match.group(1)
    name = os.path.basename(path)
    for suffix in FASTQ_EXTENSIONS:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def find_fastqs(inputdir):
    if not os.path.isdir(inputdir):
        raise SystemExit(f"ERROR: FASTQ directory not found: {inputdir}")
    fastqs = sorted(
        os.path.join(inputdir, name)
        for name in os.listdir(inputdir)
        if name.endswith(FASTQ_EXTENSIONS)
    )
    if not fastqs:
        raise SystemExit(f"ERROR: no single-end FASTQs found in {inputdir}")
    return fastqs


def align_sort_index(fastq, bamdir, index, threads, sort_threads):
    sample = sample_name(fastq)
    os.makedirs(bamdir, exist_ok=True)
    output = os.path.join(bamdir, f"{sample}.sorted.bam")
    hisat2_cmd = [
        "hisat2", "-p", str(threads), "--very-sensitive",
        "-x", index, "-U", fastq,
    ]
    sort_cmd = ["samtools", "sort", "-@", str(sort_threads), "-o", output, "-"]
    print(f"Aligning {sample} (single-end)", flush=True)
    hisat2_proc = subprocess.Popen(hisat2_cmd, stdout=subprocess.PIPE)
    sort_proc = subprocess.Popen(sort_cmd, stdin=hisat2_proc.stdout)
    hisat2_proc.stdout.close()
    sort_rc = sort_proc.wait()
    hisat2_rc = hisat2_proc.wait()
    if hisat2_rc or sort_rc:
        raise SystemExit(f"ERROR: alignment failed for {fastq}")
    subprocess.check_call(["samtools", "index", "-@", str(sort_threads), output])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()

    for command in ("hisat2", "samtools"):
        if not shutil.which(command):
            raise SystemExit(f"ERROR: {command} is not available in PATH")

    inputdir = os.environ.get(
        "SHP53_FASTQ_DIR", f"{PROJECT}/shared/rnaseqtrimmeddir/shP53_SRA_single_end"
    )
    bamdir = os.environ.get(
        "SHP53_BAM_DIR", f"{PROJECT}/shared/BAMfiles/rnaseq_bamfiles/human"
    )
    index = os.environ.get(
        "SHP53_HISAT2_INDEX", f"{PROJECT}/../Genome/hg38_with_rDNA_all"
    )
    validate_hisat2_index(index)
    fastqs = find_fastqs(inputdir)
    samples = [sample_name(path) for path in fastqs]
    if len(samples) != len(set(samples)):
        raise SystemExit("ERROR: duplicate SRR accessions would overwrite BAM files")

    total_threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "16"))
    sort_threads = min(4, max(1, (total_threads - 2) // 4))
    align_threads = max(1, total_threads - sort_threads - 1)
    if args.check_only:
        print(f"Preflight OK: {len(fastqs)} single-end FASTQs")
        return

    for fastq in fastqs:
        align_sort_index(fastq, bamdir, index, align_threads, sort_threads)


if __name__ == "__main__":
    main()
