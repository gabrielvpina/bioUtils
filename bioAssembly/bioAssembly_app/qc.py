import os, subprocess


def run_fastqc(args):
    # CREATE AN OUTPUT DIRECTORY
    """Run FastQC in samples."""
    print(f"\n--- Running Fastqc in {args.inputdir} files. ---")

    os.makedirs(os.path.join(args.outdir, 'fastqc'))

    outdir = os.path.join(args.outdir, 'fastqc')

    # run fastqc command
    cmd = [
        'fastqc',
        os.path.join(args.inputdir,'*'),
        '-t',
        str(args.threads),
        '--outdir',
        outdir
    ]

    print(f"FastQC Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    return outdir


def run_multiqc(args):
    """Run Multiqc to evaluate QC analysis."""
    print(f"\n--- Running MultiQC fastqc directory. ---")

    os.makedirs(os.path.join(args.outdir, 'multiqc'))

    outdir = os.path.join(args.outdir, 'multiqc')

    # run MultiQC command
    cmd = [
        'multiqc',
        os.path.join(args.outdir,'fastqc'),
        '-o',
        outdir
    ]

    print(f"MultiQC Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
