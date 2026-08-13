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
    if not os.path.isdir(trimmeddir):
        os.makedirs(trimmeddir)

    cmd = ["trim_galore","--paired","--illumina", "--cores", str(threads), "--output_dir", trimmeddir, r1path, r2path]
   
    print(
        "trimming {} and {}".format(os.path.basename(r1path), os.path.basename(r2path))
    )
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
        print("Trim Galore failed for {}, skipping.".format(os.path.basename(r1path)))


def runfastqc(OGdir, fastqcdir, threads=4):
    fastqs = sorted(
        os.path.join(OGdir, name)
        for name in os.listdir(OGdir)
        if name.endswith(fastqvariations)
    )
    if not fastqs:
        raise SystemExit("Error: no FASTQ files found in {}.".format(OGdir))

    if not os.path.isdir(fastqcdir):
        os.makedirs(fastqcdir)

    print("Running FastQC on FASTQ files")
    cmd = ["fastqc", "--threads", str(threads), "--outdir", fastqcdir] + fastqs
    subprocess.check_call(cmd)


def runmultiqc(fastqcdir, multiqcdir):
    if not os.path.isdir(multiqcdir):
        os.makedirs(multiqcdir)

    print("Running MultiQC")
    subprocess.check_call(["multiqc", fastqcdir, "--outdir", multiqcdir, "--force"])


def main():
    inputdir = "p53kd"
    trimmeddir = "p53kdtrimmeddir"
    fastqcdir = "p53kdfastqcdir"
    multiqcdir = "p53kdresults"



    if not os.path.isdir(trimmeddir):
        os.makedirs(trimmeddir)
    if not os.path.isdir(fastqcdir):
        os.makedirs(fastqcdir)
    if not os.path.isdir(multiqcdir):
        os.makedirs(multiqcdir)

    for r1path, r2path in findpairs(inputdir):
        runtrimgalore(r1path, r2path, trimmeddir, threads=16)
    runfastqc(trimmeddir, fastqcdir, threads=10)
    runmultiqc(fastqcdir, multiqcdir)
    print("reports are in fastqcdir/ and results/.")


if __name__ == "__main__":
    main()
