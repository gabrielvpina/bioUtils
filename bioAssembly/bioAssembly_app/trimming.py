import os, sys, subprocess
"""
Data processing and trimming.
"""

# Fastp QC routine 
# ========================================================================
def run_fastp_PE(sample, file1, file2, args):
    """Run Fastp trimming for PAIR-END samples."""
    print(f"\n--- Running Fastp for sample {sample} ---")

    if file1.endswith("1.fastq.gz"):
        pass
        if file2.endswith("2.fastq.gz"):
            pass
    else:
        print("ERROR: The fastq files need to be with '1.fastq.gz' AND '2.fastq.gz' suffixes")
        
    
    # Prepare output file paths
    out1 = os.path.join(args.trimming_dir, f"{sample}_1.fastq.gz")
    out2 = os.path.join(args.trimming_dir, f"{sample}_2.fastq.gz")
    
    # run fastp command
    cmd = [
        'fastp', 
        '--trim_poly_g',
        '--in1', file1,
        '--in2', file2,
        '--out1', out1,
        '--out2', out2
    ]
    
    print(f"Trimming {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return out1, out2


def run_fastp_SE(sample, file1, args):
    """Run Fastp trimming for SINGLE-END samples."""
    print(f"\n--- Running Fastp for sample {sample} ---")
    
    # Prepare output file paths
    out1 = os.path.join(args.trimming_dir, f"{sample}.fastq.gz")
    
    # run fastp command
    cmd = [
        'fastp', 
        '--trim_poly_g',
        '--in1', file1,
        '--out1', out1
    ]
    
    print(f"Trimming {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return out1
# ========================================================================





# Trimmomatic QC routine 
# ========================================================================
def run_trimmomatic_PE(sample, file1, file2, args):
    """Run Trimmomatic trimming for PAIR-END samples."""

    if file1.endswith("1.fastq.gz"):
        pass
        if file2.endswith("2.fastq.gz"):
            pass
    else:
        print("ERROR: The fastq files need to be with '1.fastq.gz' AND '2.fastq.gz' suffixes")


    print(f"\n--- Running Trimmomatic for sample {sample} ---")
    
    # Prepare output file paths
    out1 = os.path.join(args.trimming_dir, f"{sample}_1.fastq.gz")
    out2 = os.path.join(args.trimming_dir, f"{sample}_2.fastq.gz")
    
    # run fastp command
    cmd = [
        'Trimmomatic', 'PE',
        '-phred33',
        file1,
        file2,
        out1,
        out2,
        'ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10',
        'LEADING:3',
        'TRAILING:3',
        'SLIDINGWINDOW:4:15',
        'MINLEN:36'
    ]
    
    print(f"Trimming {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return out1, out2


def run_trimmomatic_SE(sample, file1, args):
    """Run Trimmomatic trimming for SINGLE-END samples."""
    print(f"\n--- Running Trimmomatic for sample {sample} ---")
    
    # Prepare output file paths
    out1 = os.path.join(args.trimming_dir, f"{sample}_1.fastq.gz")    
    # run fastp command
    cmd = [
        'Trimmomatic', 'SE',
        '-phred33',
        file1,
        out1,
        'ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10',
        'LEADING:3',
        'TRAILING:3',
        'SLIDINGWINDOW:4:15',
        'MINLEN:36'
    ]
    
    print(f"Trimming {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return out1
# ========================================================================
