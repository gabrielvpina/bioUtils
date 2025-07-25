#!/usr/bin/env python3
"""
Downstream analysis of assembled contigs.
Trim, map and assemble.
"""

import argparse


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Bioinformatics Pipeline')
    
    # Required arguments
    parser.add_argument('--input_dir', '-i', required=True, help='Directory containing "fastq.gz" files')
    
    # Optional arguments 
    parser.add_argument('--fastp_dir', default='fastp', help='Output directory for Fastp results')
    parser.add_argument('--star_dir', default='STAR_unmapped', help='Output directory for STAR results')
    parser.add_argument('--megahit_dir', default='megahit', help='Output directory for megahit results')
    parser.add_argument('--contigs_dir', default='assembled_contigs', help='Output directory for assembled contigs')
    parser.add_argument('--star_genome', default='star_genome', help='STAR genome directory')
    parser.add_argument('--megahit_path', default='megahit.py', help='Path to megahit executable')
    parser.add_argument('--threads', type=int, default=7, help='Number of threads to use')
    parser.add_argument('--memory', type=int, default=14, help='Memory limit for megahit in GB')
    parser.add_argument('--skip_fastp', action='store_true', help='Skip Fastp trimming step')
    parser.add_argument('--skip_star', action='store_true', help='Skip STAR alignment step')
    parser.add_argument('--skip_megahit', action='store_true', help='Skip megahit assembly step')
    parser.add_argument('--skip_contigs', action='store_true', help='Skip copying contigs step')
    parser.add_argument('--keep_temp_files', action='store_true', help='Keep temporary files after processing')
    
    return parser.parse_args()


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