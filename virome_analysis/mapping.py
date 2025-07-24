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