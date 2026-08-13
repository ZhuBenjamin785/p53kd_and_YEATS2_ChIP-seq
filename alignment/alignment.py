import argparse
import os
import shutil
import subprocess

fastqvariations = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
script_dir = os.path.dirname(os.path.abspath(__file__))


def validate_bowtie2_index(index_prefix):
    """Fail early unless index_prefix names a complete Bowtie 2 index."""
    suffixes = (".1", ".2", ".3", ".4", ".rev.1", ".rev.2")
    for extension in (".bt2", ".bt2l"):
        expected = [index_prefix + suffix + extension for suffix in suffixes]
        if all(os.path.isfile(path) for path in expected):
            return

    raise SystemExit(
        "Error: {!r} is not a complete Bowtie 2 index prefix. "
        "Expected six .bt2 or .bt2l files at that location.".format(index_prefix)
    )


def validate_executables():
    """Fail before doing any work if a required command is unavailable."""
    missing = [command for command in ("bowtie2", "samtools") if not shutil.which(command)]
    if missing:
        raise SystemExit(
            "Error: required command(s) not found in PATH: {}.".format(
                ", ".join(missing)
            )
        )


def validate_input_files(pairs):
    """Reject empty or unreadable FASTQ inputs before starting alignment."""
    invalid = []
    for pair in pairs:
        for path in pair:
            if not os.path.isfile(path) or not os.access(path, os.R_OK) or os.path.getsize(path) == 0:
                invalid.append(path)
    if invalid:
        raise SystemExit(
            "Error: input FASTQ files are missing, unreadable, or empty: {}".format(
                ", ".join(sorted(invalid))
            )
        )


def findpairs(inputdir):
    if not os.path.isdir(inputdir):
        raise SystemExit("Error: FASTQ input directory not found: {}.".format(inputdir))

    filenames = set(
        name for name in os.listdir(inputdir) if name.endswith(fastqvariations)
    )
    pairs = []
    for r1name in sorted(
        name for name in filenames if "_R1_" in name and "_val_1" in name
    ):
        r2name = r1name.replace("_R1_", "_R2_", 1).replace(
            "_val_1", "_val_2", 1
        )
        if r2name not in filenames:
            raise SystemExit("Error: could not find reverse read for {}.".format(r1name))
        pairs.append((os.path.join(inputdir, r1name), os.path.join(inputdir, r2name)))
    if not pairs:
        raise SystemExit(
            "Error: no finalized trimmed pairs named with _R1_..._val_1 and "
            "_R2_..._val_2 found in {}.".format(inputdir)
        )
    return pairs


def sample_name_from_r1(r1path):
    name = os.path.basename(r1path)
    for suffix in fastqvariations:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    for marker in ("_R1_", "_R1", "_val_1", "_val_1_001"):
        if marker in name:
            name = name.split(marker, 1)[0]
            break
    return name


def align_sort_index(
    r1,
    r2,
    BAMfiles,
    genome_index,
    genome_name,
    alignment_threads,
    sort_threads,
):
    """Stream Bowtie 2 SAM output into a coordinate-sorted, indexed BAM."""
    sample = sample_name_from_r1(r1)
    genome_dir = os.path.join(BAMfiles, genome_name)
    sorted_bam_path = os.path.join(genome_dir, sample + ".sorted.bam")

    if not os.path.isdir(genome_dir):
        os.makedirs(genome_dir)
    print(
        "Aligning {} against {} and writing {}".format(
            sample, genome_name, sorted_bam_path
        ),
        flush=True,
    )

    bowtie_cmd = [
        "bowtie2",
        "-p", str(alignment_threads),
        "--very-sensitive",
        "-X", "2000",
        "--no-mixed",
        "--no-discordant",
        "-x", genome_index,
        "-1", r1,
        "-2", r2,
    ]
    sort_cmd = [
        "samtools",
        "sort",
        "-@",
        str(sort_threads),
        "-o",
        sorted_bam_path,
        "-",
    ]

    bowtie_process = subprocess.Popen(bowtie_cmd, stdout=subprocess.PIPE)
    try:
        sort_process = subprocess.Popen(sort_cmd, stdin=bowtie_process.stdout)
        bowtie_process.stdout.close()
        sort_returncode = sort_process.wait()
        bowtie_returncode = bowtie_process.wait()
    except BaseException:
        bowtie_process.kill()
        bowtie_process.wait()
        raise

    if bowtie_returncode != 0:
        raise subprocess.CalledProcessError(bowtie_returncode, bowtie_cmd)
    if sort_returncode != 0:
        raise subprocess.CalledProcessError(sort_returncode, sort_cmd)

    subprocess.check_call([
        "samtools",
        "index",
        "-@",
        str(sort_threads),
        sorted_bam_path,
    ])


def parse_args():
    parser = argparse.ArgumentParser(
        description="Align trimmed p53kd paired-end reads to human and dm6."
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="validate tools, inputs, sample pairing, and indexes without running alignment",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    inputdir = "p53kdtrimmeddir"
    BAMfiles = "p53kdbamfiles"
    genomes = [
        (
            os.environ.get(
                "HUMAN_BOWTIE2_INDEX",
                os.path.join(script_dir, "index", "human", "hgenome_index"),
            ),
            "human",
        ),
        (
            os.environ.get(
                "DM6_BOWTIE2_INDEX",
                os.path.join(script_dir, "index", "dm6", "dm6genome_index"),
            ),
            "dm6",
        ),
    ]

    validate_executables()
    for genome_index, _ in genomes:
        validate_bowtie2_index(genome_index)

    pairs = findpairs(inputdir)
    validate_input_files(pairs)
    samples = [sample_name_from_r1(r1path) for r1path, _ in pairs]
    if len(samples) != len(set(samples)):
        raise SystemExit(
            "Error: multiple read pairs resolve to the same sample name; output "
            "files would overwrite one another."
        )

    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "16"))
    if threads < 2:
        raise SystemExit("Error: SLURM_CPUS_PER_TASK must be at least 2.")
                                                                   
    sort_threads = min(4, max(0, (threads - 2) // 4))
    alignment_threads = threads - sort_threads - 1

    if args.check_only:
        print(
            "Preflight OK: {} paired samples, {} genomes, {} Bowtie 2 threads, "
            "{} additional samtools threads.".format(
                len(pairs), len(genomes), alignment_threads, sort_threads
            )
        )
        return

    for r1path, r2path in pairs:
        for genome_index, genome_name in genomes:
            align_sort_index(
                r1path,
                r2path,
                BAMfiles,
                genome_index,
                genome_name,
                alignment_threads,
                sort_threads,
            )
    

if __name__ == "__main__":
    main()
