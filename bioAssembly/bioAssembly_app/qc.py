import os, subprocess


def run_fastqc(sample, args):
    # CREATE AN OUTPUT DIRECTORY
    """Run Fastp trimming for PAIR-END samples."""
    print(f"\n--- Running Fastp for sample {sample} ---")

    outdir = os.path.join(args.fastqc_dir, f"{sample}_1.fastq.gz")

    # run fastqc command
    cmd = [
        'fastqc',
        os.path.join(args.trimming_dir,'*'),
        '-t',
        str(args.threads),
        '--outdir',
        outdir
    ]

    print(f"Trimming {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    return outdir


def run_multiqc(args):
    """Run Multiqc to evaluate QC analysis."""
    print(f"\n--- Running MultiQC for sample {sample} ---")

    outdir = os.path.join(args.fastqc_dir, f"{sample}_1.fastq.gz")

    # run MultiQC    command
    cmd = [
        'fastqc',
        os.path.join(args.trimming_dir,'*'),
        '-t',
        str(args.threads),
        '-o',
        outdir
    ]

    print(f"Trimming {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
