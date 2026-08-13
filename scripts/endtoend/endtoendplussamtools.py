import os
import subprocess
import __hello__
import shutil
import glob

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


def sample_name_from_r1(r1path):
    name = os.path.basename(r1path)
    for suffix in fastqvariations:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    for marker in ("_R1_", "_R1"):
        if marker in name:
            name = name.split(marker, 1)[0]
            break
    return name


def runtrimgalore(r1path, r2path, trimmeddir, threads=10):
    if not os.path.isdir(trimmeddir):
        os.makedirs(trimmeddir)

    cmd = ["trim_galore","--paired","--illumina", "--cores", str(threads), "--output_dir", trimmeddir, r1path, r2path]
   
    print(
        "trimming {} and {}".format(os.path.basename(r1path), os.path.basename(r2path))
    )
    subprocess.check_call(cmd)


def runfastqc(OGdir, fastqcdir, threads=10):
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
    
    
def runhisat2(r1, r2, alignmentresults, threads=16):
    sample = sample_name_from_r1(r1)
    sam_path = os.path.join(alignmentresults, sample + ".sam")

    print("running HISAT2 alignment now")
    
    if not os.path.isdir(alignmentresults):
        os.makedirs(alignmentresults)
    
    cmd = [
        "hisat2",
        "-p",
        str(threads),
        "-x",
        "genome_index",
        "-1",
        r1,
        "-2",
        r2,
        "-S",
        sam_path,
    ]
    
    subprocess.check_call(cmd)

def runSAMtoBAM(alignmentresults, BAMfiles, r1path, threads=16):
    print("running samtools to convert SAM to BAM")

    sample = sample_name_from_r1(r1path)
    sam_path = os.path.join(alignmentresults, sample + ".sam")
    bam_path = os.path.join(BAMfiles, sample + ".bam")

    if not os.path.isdir(BAMfiles):
        os.makedirs(BAMfiles)
    
    cmd = [
        "samtools",
        "view",
        "-@",
        str(threads),
        "-b",
        "-o",
        bam_path,
        sam_path,
    ]
    subprocess.check_call(cmd)

def SortAndIndex(BAMfiles, r1path, threads=10):
    print("running samtools to convert SAM to BAM")

    sample = sample_name_from_r1(r1path)
    bam_path = os.path.join(BAMfiles, sample + ".bam")
    sorted_bam_path = os.path.join(BAMfiles, sample + ".sorted.bam")

    cmd1 = [
        "samtools",
        "sort",
        "-@",
        str(threads),
        "-o",
        sorted_bam_path,
        bam_path,
    ]
    subprocess.check_call(cmd1)

    cmd2 = [
        "samtools",
        "index",
        "-@",
        str(threads),
        sorted_bam_path,
    ]
    subprocess.check_call(cmd2)


def main():
    inputdir = "FASTQ"
    trimmeddir = "trimmeddir"
    fastqcdir = "fastqcdir"
    multiqcdir = "results"
    alignmentresults = "alignmentresults"
    BAMfiles = "BAMfiles"


    if not os.path.isdir(trimmeddir):
        os.makedirs(trimmeddir)
    if not os.path.isdir(fastqcdir):
        os.makedirs(fastqcdir)
    if not os.path.isdir(multiqcdir):
        os.makedirs(multiqcdir)

    for r1path, r2path in findpairs(inputdir):
        runtrimgalore(r1path, r2path, trimmeddir, threads=10)
        runhisat2(r1path, r2path, alignmentresults, threads = 16)
    
    runfastqc(trimmeddir, fastqcdir, threads=10)
    runmultiqc(fastqcdir, multiqcdir)
    for r1path, _ in findpairs(inputdir):
        runSAMtoBAM(alignmentresults, BAMfiles, r1path, threads=16)
        SortAndIndex(BAMfiles, r1path, threads=10)
    

if __name__ == "__main__":
    main()
