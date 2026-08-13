#!/usr/bin/env python3
                      
"""Align fastqchip single-end FASTQs to the dm6 Bowtie2 index."""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path


project_dir = Path(__file__).resolve().parent
fastq_dir = project_dir / "fastqchip"
bam_dir = project_dir / "fastqchip_bamfiles" / "dm6"
default_index = project_dir / "index" / "dm6" / "dm6genome_index"
fastq_suffixes = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
target_samples = {"SRR5944063", "SRR5944064", "SRR5944081", "SRR5944082"}


def validate_index(prefix):
    prefix = str(prefix)
    suffixes = (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")
    for extension in (".bt2", ".bt2l"):
        if all(Path(prefix + suffix + extension).is_file() for suffix in suffixes):
            return
    raise SystemExit(f"Incomplete Bowtie2 index: {prefix}")


def sample_name(path):
    name = path.name
    for suffix in fastq_suffixes:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.endswith("_1"):
        name = name[:-2]
    return name


def main():
    index = Path(os.environ.get("DM6_BOWTIE2_INDEX", default_index))
    if not fastq_dir.is_dir():
        raise SystemExit(f"FASTQ directory not found: {FASTQ_DIR}")
    if not shutil.which("bowtie2") or not shutil.which("samtools"):
        raise SystemExit("bowtie2 and samtools must be available in PATH.")
    validate_index(index)

    fastqs = sorted(
        p for p in fastq_dir.iterdir()
        if p.is_file()
        and p.name.endswith(fastq_suffixes)
        and sample_name(p) in target_samples
    )
    if not fastqs:
        raise SystemExit(f"No FASTQs found in {FASTQ_DIR}")
    samples = [sample_name(p) for p in fastqs]
    if len(samples) != len(set(samples)):
        raise SystemExit("Duplicate sample names would overwrite BAM outputs.")

    bam_dir.mkdir(parents=True, exist_ok=True)
    total_threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "16"))
    sort_threads = min(4, max(1, (total_threads - 2) // 4))
    bowtie2_threads = max(1, total_threads - sort_threads - 1)

    for fastq, sample in zip(fastqs, samples):
        bam = bam_dir / f"{sample}.sorted.bam"
        bai = Path(str(bam) + ".bai")
                                                                         
                                                                     
        if bam.is_file() and bai.is_file():
            check = subprocess.run(
                ["samtools", "quickcheck", "-v", str(bam)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if check.returncode == 0:
                print(f"Reusing complete dm6 BAM for {sample}", flush=True)
                continue
            bam.unlink(missing_ok=True)
            bai.unlink(missing_ok=True)

                                                                             
                                                                     
        token = f"{os.getpid()}-{next(tempfile._get_candidate_names())}"
        partial_bam = bam_dir / f".{sample}.{token}.partial.bam"
        temp_prefix = bam_dir / f".{sample}.{token}.sorttmp"
        print(f"Aligning {sample} to dm6", flush=True)
        bowtie2_cmd = [
            "bowtie2", "--very-sensitive", "-p", str(bowtie2_threads),
            "-x", str(index), "-U", str(fastq),
        ]
        sort_cmd = [
            "samtools", "sort", "-@", str(sort_threads),
            "-T", str(temp_prefix), "-o", str(partial_bam), "-",
        ]
        bowtie2 = subprocess.Popen(bowtie2_cmd, stdout=subprocess.PIPE)
        sorter = subprocess.Popen(sort_cmd, stdin=bowtie2.stdout)
        bowtie2.stdout.close()
        sort_code = sorter.wait()
        bowtie2_code = bowtie2.wait()
        if bowtie2_code or sort_code:
            partial_bam.unlink(missing_ok=True)
            raise SystemExit(f"Alignment failed for {fastq}")
        subprocess.check_call(["samtools", "quickcheck", "-v", str(partial_bam)])
        os.replace(partial_bam, bam)
        subprocess.check_call(["samtools", "index", "-@", str(sort_threads), str(bam)])


if __name__ == "__main__":
    main()
