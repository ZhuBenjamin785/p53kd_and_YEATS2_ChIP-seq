import argparse
import os
import shutil
import subprocess

fastq_extensions = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
script_dir = os.path.dirname(os.path.abspath(__file__))


def validate_hisat2_index(prefix):
    suffixes = (".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8")
    for extension in (".ht2", ".ht2l"):
        if all(os.path.isfile(prefix + suffix + extension) for suffix in suffixes):
            return
    raise SystemExit("Error: incomplete HISAT2 index: {}.".format(prefix))


def sample_name(path):
    name = os.path.basename(path)
    for suffix in fastq_extensions:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    for marker in ("_val_1", "_R1", "_1", "_trimmed"):
        if marker in name:
            name = name.split(marker, 1)[0]
            break
    return name


def find_fastqs(inputdir):
    if not os.path.isdir(inputdir):
        raise SystemExit("Error: trimmed FASTQ directory not found: {}.".format(inputdir))
    fastqs = sorted(
        os.path.join(inputdir, name)
        for name in os.listdir(inputdir)
        if name.endswith(fastq_extensions) and "_val_" not in name
    )
    if not fastqs:
        raise SystemExit("Error: no single-end trimmed FASTQ files found in {}.".format(inputdir))
    return fastqs


def align_sort_index(fastq, bamdir, index, genome, alignment_threads, sort_threads):
    outdir = os.path.join(bamdir, genome)
    os.makedirs(outdir, exist_ok=True)
    output = os.path.join(outdir, sample_name(fastq) + ".sorted.bam")
    hisat2 = ["hisat2", "-p", str(alignment_threads), "--very-sensitive", "-x", index, "-U", fastq]
    sort = ["samtools", "sort", "-@", str(sort_threads), "-o", output, "-"]
    print("Aligning {} against {}".format(sample_name(fastq), genome), flush=True)
    hisat2_proc = subprocess.Popen(hisat2, stdout=subprocess.PIPE)
    sort_proc = subprocess.Popen(sort, stdin=hisat2_proc.stdout)
    hisat2_proc.stdout.close()
    sort_rc = sort_proc.wait()
    hisat2_rc = hisat2_proc.wait()
    if hisat2_rc or sort_rc:
        raise SystemExit("Error: alignment failed for {}.".format(fastq))
    subprocess.check_call(["samtools", "index", "-@", str(sort_threads), output])


def main():
    parser = argparse.ArgumentParser(description="Align single-end RNA-seq reads to the human genome.")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if not shutil.which("hisat2") or not shutil.which("samtools"):
        raise SystemExit("Error: hisat2 and samtools must be available in PATH.")
    genomes = [
        (
            os.environ.get(
                "HUMAN_HISAT2_INDEX",
                "/gpfs/projects/b1042/LauberthLab/Genome/hg38/genome",
            ),
            "human",
        ),
    ]
    for index, _ in genomes:
        validate_hisat2_index(index)
    fastqs = find_fastqs("rnaseqtrimmeddir")
    samples = [sample_name(path) for path in fastqs]
    if len(samples) != len(set(samples)):
        raise SystemExit("Error: duplicate RNA-seq sample names would overwrite BAM files.")
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "16"))
    sort_threads = min(4, max(1, (threads - 2) // 4))
    alignment_threads = max(1, threads - sort_threads - 1)
    if args.check_only:
        print("Preflight OK: {} single-end samples, {} genomes.".format(len(fastqs), len(genomes)))
        return
    for fastq in fastqs:
        for index, genome in genomes:
            align_sort_index(fastq, "rnaseq_bamfiles", index, genome, alignment_threads, sort_threads)


if __name__ == "__main__":
    main()
