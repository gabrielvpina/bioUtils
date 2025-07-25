import os, sys, subprocess
"""
Execute mapping tools to remove host sequences.
"""

# run STAR for pair-end data
def run_star_PE(sample, file1, file2, args):

    """Run STAR alignment for a pair-end sample."""

    # check if star index exist
    if not os.path.exists(args.star_genome):
        print(f"The {args.star_genome} directory does not exist.")
        sys.exit(1)
    else:
        print(f"\n--- Running STAR in pair-end mode for sample {sample} ---")
    
    # temporary uncompressed files
    temp_file1 = os.path.join(args.fastp_dir, f"{sample}_1.fastq")
    temp_file2 = os.path.join(args.fastp_dir, f"{sample}_2.fastq")
    
    # decompress files
    with open(temp_file1, 'wb') as outfile:
        subprocess.run(['gunzip', '-c', file1], stdout=outfile, check=True)
    
    with open(temp_file2, 'wb') as outfile:
        subprocess.run(['gunzip', '-c', file2], stdout=outfile, check=True)
    
    # STAR command
    star_prefix = os.path.join(args.star_dir, f"{sample}_")
    cmd = [
        'STAR',
        '--runThreadN', str(args.threads),
        '--genomeDir', args.star_genome,
        '--readFilesIn', temp_file1, temp_file2,
        '--outFileNamePrefix', star_prefix,
        '--outReadsUnmapped', 'Fastx'
    ]
    
    print(f"Running STAR for {sample}")
    subprocess.run(cmd, check=True)
    
    # clean temporary uncompressed files
    os.remove(temp_file1)
    os.remove(temp_file2)
    
    # paths to unmapped reads
    unmapped1 = os.path.join(args.star_dir, f"{sample}_Unmapped.out.mate1")
    unmapped2 = os.path.join(args.star_dir, f"{sample}_Unmapped.out.mate2")
    
    # Compress unmapped reads
    unmapped1_gz = f"{unmapped1}.gz"
    unmapped2_gz = f"{unmapped2}.gz"
    
    with open(unmapped1, 'rb') as infile, open(unmapped1_gz, 'wb') as outfile:
        subprocess.run(['gzip', '-c'], stdin=infile, stdout=outfile, check=True)
    
    with open(unmapped2, 'rb') as infile, open(unmapped2_gz, 'wb') as outfile:
        subprocess.run(['gzip', '-c'], stdin=infile, stdout=outfile, check=True)
    
    # remove uncompressedd unmapped files 
    os.remove(unmapped1)
    os.remove(unmapped2)
    
    return unmapped1_gz, unmapped2_gz




def run_star_SE(sample, file1, args):

    """Run STAR alignment for a single-end sample."""

    # check if star index exist
    if not os.path.exists(args.star_genome):
        print(f"The {args.star_genome} directory does not exist.")
        sys.exit(1)
    else:
        print(f"\n--- Running STAR in pair-end mode for sample {sample} ---")
    
    # temporary uncompressed files
    temp_file1 = os.path.join(args.fastp_dir, f"{sample}_1.fastq")
    
    # decompress files
    with open(temp_file1, 'wb') as outfile:
        subprocess.run(['gunzip', '-c', file1], stdout=outfile, check=True)
    
    
    # STAR command
    star_prefix = os.path.join(args.star_dir, f"{sample}_")
    cmd = [
        'STAR',
        '--runThreadN', str(args.threads),
        '--genomeDir', args.star_genome,
        '--readFilesIn', temp_file1,
        '--outFileNamePrefix', star_prefix,
        '--outReadsUnmapped', 'Fastx'
    ]
    
    print(f"Running STAR for {sample}")
    subprocess.run(cmd, check=True)
    
    # clean temporary uncompressed files
    os.remove(temp_file1)
    
    # paths to unmapped reads
    unmapped1 = os.path.join(args.star_dir, f"{sample}_Unmapped.out.mate1")
    
    # Compress unmapped reads
    unmapped1_gz = f"{unmapped1}.gz"
    
    with open(unmapped1, 'rb') as infile, open(unmapped1_gz, 'wb') as outfile:
        subprocess.run(['gzip', '-c'], stdin=infile, stdout=outfile, check=True)
    
    # remove uncompressedd unmapped files 
    os.remove(unmapped1)
    
    return unmapped1_gz







# Quantification routine - a test
def run_star_quantification(sample, file1, file2, args):
    """
    Run STAR alignment and gene quantification for a pair-end sample.

    Args:
        sample (str): The name of the sample.
        file1 (str): Path to the first input FASTQ file (e.g., R1.fastq.gz).
        file2 (str): Path to the second input FASTQ file (e.g., R2.fastq.gz).
        args (object): An object containing command-line arguments or configurations.
                       Expected attributes:
                       - args.star_genome (str): Path to the STAR genome index directory.
                       - args.gtf_file (str): Path to the GTF/GFF annotation file used for quantification.
                       - args.star_quant_dir (str): Directory where quantification output files will be saved.
                       - args.threads (int): Number of threads to use for STAR.

    Returns:
        tuple: A tuple containing paths to the generated gene counts file and sorted BAM file.
               (gene_counts_file_path, sorted_bam_file_path)
    """

    # Check if STAR genome index directory exists
    if not os.path.exists(args.star_genome):
        print(f"Error: The STAR genome directory '{args.star_genome}' does not exist.")
        sys.exit(1)

    # Check if GTF file exists
    if not os.path.exists(args.gtf_file):
        print(f"Error: The GTF/GFF annotation file '{args.gtf_file}' does not exist.")
        sys.exit(1)

    # Create the output directory if it doesn't exist
    os.makedirs(args.star_quant_dir, exist_ok=True)

    print(f"\n--- Running STAR in quantification mode for sample {sample} ---")

    # Define the prefix for output files
    star_prefix = os.path.join(args.star_quant_dir, f"{sample}_")

    # STAR command for alignment and quantification
    # --runThreadN: Number of threads
    # --genomeDir: Path to the STAR genome index
    # --readFilesIn: Input FASTQ files (can be gzipped)
    # --outFileNamePrefix: Prefix for all output files
    # --outSAMtype BAM SortedByCoordinate: Output sorted BAM file
    # --quantMode GeneCounts: Perform gene quantification and output ReadsPerGene.out.tab
    # --sjdbGTFfile: Path to the GTF/GFF annotation file for quantification
    cmd = [
        'STAR',
        '--runThreadN', str(args.threads),
        '--genomeDir', args.star_genome,
        '--readFilesIn', file1, file2,
        '--outFileNamePrHow integrate multiq  args.gtf_file
    ]

    print(f"Running STAR for quantification of {sample}")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"STAR quantification for {sample} completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running STAR for {sample}:")
        print(f"Command: {' '.join(e.cmd)}")
        print(f"Return Code: {e.returncode}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        sys.exit(1)

    # Paths to the main output files
    gene_counts_file = os.path.join(args.star_quant_dir, f"{sample}_ReadsPerGene.out.tab")
    sorted_bam_file = os.path.join(args.star_quant_dir, f"{sample}_Aligned.sortedByCoordinate.out.bam")
    log_final_file = os.path.join(args.star_quant_dir, f"{sample}_Log.final.out")

    # Check if expected output files exist
    if not os.path.exists(gene_counts_file):
        print(f"Warning: Gene counts file '{gene_counts_file}' not found after STAR run.")
    if not os.path.exists(sorted_bam_file):
        print(f"Warning: Sorted BAM file '{sorted_bam_file}' not found after STAR run.")
    if not os.path.exists(log_final_file):
        print(f"Warning: Log file '{log_final_file}' not found after STAR run.")

    print(f"Output gene counts file: {gene_counts_file}")
    print(f"Output sorted BAM file: {sorted_bam_file}")

    return gene_counts_file, sorted_bam_file