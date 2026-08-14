import os
import subprocess


fastqvariations = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def findpairs(inputdir):
    filenames = set(
        name for name in os.listdir(inputdir) if name.endswith(fastqvariations)
    )
    pairs = []
    for r1name in sorted(name for name in filenames if "_R1_" in name):
        r2name = r1name.replace("_R1_", "_R2_", 1)
        if r2name not in filenames:
            raise SystemExit("error: could not find reverse for {}".format(r1name))
        pairs.append((os.path.join(inputdir, r1name), os.path.join(inputdir, r2name)))
    if not pairs:
        raise SystemExit(
            "Error: no paired files named with _R1_/_R2_ found in {}.".format(inputdir)
        )
    return pairs


def runtrimgalore(r1path, r2path, trimmeddir, threads=4):
    os.makedirs(trimmeddir, exist_ok=True)

    cmd = ["trim_galore","--paired","--illumina", "--cores", str(threads), "--output_dir", trimmeddir, r1path, r2path]
   
    print(
        "trimming {} and {}".format(os.path.basename(r1path), os.path.basename(r2path))
    )
    subprocess.check_call(cmd)


def runfastqc(ogdir, fastqcdir, threads=4):
    fastqs = sorted(
        os.path.join(ogdir, name)
        for name in os.listdir(ogdir)
        if name.endswith(fastqvariations)
    )
    if not fastqs:
        raise SystemExit("Error: no FASTQ files found in {}.".format(ogdir))

    os.makedirs(fastqcdir, exist_ok=True)

    print("Running FastQC on FASTQ files")
    cmd = ["fastqc", "--threads", str(threads), "--outdir", fastqcdir] + fastqs
    subprocess.check_call(cmd)


def runmultiqc(fastqcdir, multiqcdir):
    os.makedirs(multiqcdir, exist_ok=True)

    print("Running MultiQC")
    subprocess.check_call(["multiqc", fastqcdir, "--outdir", multiqcdir, "--force"])


def main():
    inputdir = os.environ.get("MOF_FASTQ_DIR", "MOF")
    trimmeddir = os.environ.get("MOF_TRIMMED_DIR", "MOF_trimmed")
    fastqcdir = os.environ.get("MOF_FASTQC_DIR", "MOF_fastqcdir")
    multiqcdir = os.environ.get("MOF_MULTIQC_DIR", "MOF_results")

    os.makedirs(trimmeddir, exist_ok=True)
    os.makedirs(fastqcdir, exist_ok=True)
    os.makedirs(multiqcdir, exist_ok=True)

    for r1path, r2path in findpairs(inputdir):
        runtrimgalore(r1path, r2path, trimmeddir, threads=16)
    runfastqc(trimmeddir, fastqcdir, threads=10)
    runmultiqc(fastqcdir, multiqcdir)
    print("MOF QC reports are in {} and {}.".format(fastqcdir, multiqcdir))


if __name__ == "__main__":
    main()
