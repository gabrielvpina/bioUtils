import os, sys, subprocess
"""
Data processing and trimming.
"""


# check if files exist

def ensure_directories(dirs):
    """Create directories if they don't exist."""
    for directory in dirs:
        os.makedirs(directory, exist_ok=True)
        print(f"Ensured directory exists: {directory}")


def ensure_libs(args):
    for libs in os.listdir(args.fastq_dir):
        if libs.endswith(".fastq.gz"):
            pass
        else:
            print("ERROR: The fastq files need to be in *fastq.gz* format.")
            sys.exit(1)




# Fastp QC routine 
# ========================================================================
def check_fastp_is_installed():
    try:
        subprocess.run(['fastp', '--version'], check=True)
        return True
    except subprocess.CalledProcessError:
        print("ERROR: Fastp is not installed or not in your PATH.")
        sys.exit(1)


def run_fastp_PE(sample, file1, file2, args):
    """Run Fastp trimming for PAIR-END samples."""
    print(f"\n--- Running Fastp for sample {sample} ---")

    if file1.endswith("_1.fastq.gz"):
        pass
        if file2.endswith("_2.fastq.gz"):
            pass
    else:
        print("ERROR: The fastq files need to be with *_1.fastq.gz* AND *_2.fastq.gz* suffixes")
        
    
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
def check_trimmomatic_is_installed():
    try:
        subprocess.run(['trimmomatic', '--version'], check=True)
        return True
    except subprocess.CalledProcessError:
        print("ERROR: Trimmomatic is not installed or not in your PATH.")
        sys.exit(1)


def run_trimmomatic_PE(sample, file1, file2, args):
    """Run Trimmomatic trimming for PAIR-END samples."""
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
