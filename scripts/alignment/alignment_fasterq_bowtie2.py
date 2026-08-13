#!/usr/bin/env python3
                      
"""Align single-end FASTQs produced by fasterq.sh with Bowtie2."""

import argparse
import os
import shutil
import subprocess
from pathlib import Path


project_dir = Path(__file__).resolve().parent
fastq_dir = project_dir / "fastqchip"
bam_dir = project_dir / "fastqchip_bamfiles" / "human"
default_index = project_dir / "index" / "human" / "hgenome_index"
fastq_suffixes = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
target_samples = {"SRR5944063", "SRR5944064", "SRR5944081", "SRR5944082"}


def validate_index(prefix):
    prefix = str(prefix)
    suffixes = (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")
    for extension in (".bt2", ".bt2l"):
        if all(Path(prefix + suffix + extension).is_file() for suffix in suffixes):
            return
    raise SystemExit("Incomplete Bowtie2 index: {}".format(prefix))


def sample_name(path):
    name = path.name
    for suffix in fastq_suffixes:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.endswith("_1"):
        name = name[:-2]
    return name


def find_fastqs(input_dir):
    fastqs = sorted(
        path for path in input_dir.iterdir()
        if path.is_file()
        and path.name.endswith(fastq_suffixes)
        and sample_name(path) in target_samples
    )
    if not fastqs:
        raise SystemExit("No FASTQ files found in {}".format(input_dir))
    samples = [sample_name(path) for path in fastqs]
    if len(samples) != len(set(samples)):
        raise SystemExit("Duplicate sample names would overwrite BAM outputs.")
    return fastqs


def align(fastq, index, threads, sort_threads):
    bam_dir.mkdir(parents=True, exist_ok=True)
    sample = sample_name(fastq)
    bam = bam_dir / "{}.sorted.bam".format(sample)
    bowtie2_cmd = [
        "bowtie2", "--very-sensitive", "-p", str(threads), "-x", str(index), "-U", str(fastq),
    ]
    sort_cmd = ["samtools", "sort", "-@", str(sort_threads), "-o", str(bam), "-"]
    print("Aligning {}".format(sample), flush=True)
    bowtie2 = subprocess.Popen(bowtie2_cmd, stdout=subprocess.PIPE)
    sorter = subprocess.Popen(sort_cmd, stdin=bowtie2.stdout)
    bowtie2.stdout.close()
    sort_code = sorter.wait()
    bowtie2_code = bowtie2.wait()
    if bowtie2_code or sort_code:
        raise SystemExit("Alignment failed for {}".format(fastq))
    subprocess.check_call(["samtools", "index", "-@", str(sort_threads), str(bam)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=Path, default=Path(os.environ.get("HUMAN_BOWTIE2_INDEX", default_index)))
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if not fastq_dir.is_dir():
        raise SystemExit("FASTQ directory not found: {}".format(fastq_dir))
    if not shutil.which("bowtie2") or not shutil.which("samtools"):
        raise SystemExit("bowtie2 and samtools must be available in PATH.")
    validate_index(args.index)
    fastqs = find_fastqs(fastq_dir)
    total_threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "18"))
    sort_threads = min(4, max(1, (total_threads - 2) // 4))
    bowtie2_threads = max(1, total_threads - sort_threads - 1)
    if args.check_only:
        print("Preflight OK: {} single-end FASTQs; index: {}".format(len(fastqs), args.index))
        return
    for fastq in fastqs:
        align(fastq, args.index, bowtie2_threads, sort_threads)


if __name__ == "__main__":
    main()
