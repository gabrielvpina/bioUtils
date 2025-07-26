#!/usr/bin/env python3
"""
Downstream analysis of assembled contigs.
Trim, map and assemble.
"""

import argparse, os, sys

def parse_arguments():
    """Parse command line arguments."""

    parser = argparse.ArgumentParser(
    add_help=False,  
    formatter_class=argparse.RawDescriptionHelpFormatter)

    # Required arguments
    parser.add_argument('--inputdir', '-i', required=True)
    parser.add_argument('--outdir', '-o', required=True)
    
    # Modules arguments 
    parser.add_argument('--skip-qc', action='store_true')
    parser.add_argument('--skip-trim', action='store_true')
    parser.add_argument('--trim', choices=['fastp', 'trimmomatic'], required=True)
    parser.add_argument('--skip-align', action='store_true')
    parser.add_argument('--assembly', choices=['megahit', 'spades', 'rnaspades', 'trinity'], required=True)
    parser.add_argument('--integrative-assembly', default='MT') 
    # the initial letter of each assembly tool to be combined: 'MT', SRT', 'MSRT', etc.

    # Temporary files - By default all the temp files are saved after processing
    parser.add_argument('--remove-all', action='store_true')
    parser.add_argument('--remove-trimmed', action='store_true')
    parser.add_argument('--remove-aligned', action='store_true')
    parser.add_argument('--remove-assembled', action='store_true')

    # Config args
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--memory', type=int, default=8)
    
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