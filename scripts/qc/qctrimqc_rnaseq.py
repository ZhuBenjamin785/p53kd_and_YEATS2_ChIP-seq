import os
import subprocess


fastq_extensions = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def find_fastqs(inputdir):
    if not os.path.isdir(inputdir):
        raise SystemExit("Error: FASTQ input directory not found: {}.".format(inputdir))
    fastqs = sorted(
        os.path.join(inputdir, name)
        for name in os.listdir(inputdir)
        if name.endswith(fastq_extensions)
    )
    if not fastqs:
        raise SystemExit("Error: no FASTQ files found in {}.".format(inputdir))
    return fastqs


def runtrimgalore(fastq, trimmeddir, threads=4):
    os.makedirs(trimmeddir, exist_ok=True)
    cmd = [
        "trim_galore", "--illumina", "--cores", str(threads),
        "--output_dir", trimmeddir, fastq,
    ]
    print("Trimming {}".format(os.path.basename(fastq)), flush=True)
    subprocess.check_call(cmd)


def runfastqc(fastqdir, fastqcdir, threads=4):
    os.makedirs(fastqcdir, exist_ok=True)
    fastqs = find_fastqs(fastqdir)
    print("Running FastQC on {} FASTQ files".format(len(fastqs)), flush=True)
    subprocess.check_call(["fastqc", "--threads", str(threads), "--outdir", fastqcdir] + fastqs)


def runmultiqc(fastqcdir, multiqcdir):
    os.makedirs(multiqcdir, exist_ok=True)
    print("Running MultiQC", flush=True)
    subprocess.check_call(["multiqc", fastqcdir, "--outdir", multiqcdir, "--force"])


def main():
    inputdir = "fastq_rnaseq"
    trimmeddir = "rnaseqtrimmeddir"
    fastqcdir = "rnaseqfastqcdir"
    multiqcdir = "rnaseqresults"
    fastqs = find_fastqs(inputdir)
    for fastq in fastqs:
        runtrimgalore(fastq, trimmeddir, threads=8)
    runfastqc(trimmeddir, fastqcdir, threads=10)
    runmultiqc(fastqcdir, multiqcdir)
    print("RNA-seq QC reports are in {}/ and {}/.".format(fastqcdir, multiqcdir))


if __name__ == "__main__":
    main()
