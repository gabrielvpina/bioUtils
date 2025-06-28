#!/usr/bin/env python3
"""
Contig assembly for pair-end libraries.
Suffix to read files: "_1.fastq.gz" and "_2.fastq.gz"

This script performs the pipeline:
1. Fastp trimming
2. STAR alignment with unmapped read extraction
3. megahit assembly
4. Extraction of contigs

Usage:
    python bioinformatics_pipeline.py --fastq_dir FASTQ_DIR [options]

Date: May 14, 2025
"""

import os
import sys
import subprocess
import argparse
import glob
import shutil
from pathlib import Path


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Bioinformatics Pipeline')
    
    # Required arguments
    parser.add_argument('--fastq_dir', required=True, help='Directory containing FASTQ files')
    
    # Optional arguments with defaults
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


def ensure_directories(dirs):
    """Create directories if they don't exist."""
    for directory in dirs:
        os.makedirs(directory, exist_ok=True)
        print(f"Ensured directory exists: {directory}")


def run_fastp(sample, file1, file2, args):
    """Run Fastp trimming for a single sample."""
    print(f"\n--- Running Fastp for sample {sample} ---")
    
    # Prepare output file paths
    out1 = os.path.join(args.fastp_dir, f"{sample}_1.fastq.gz")
    out2 = os.path.join(args.fastp_dir, f"{sample}_2.fastq.gz")
    
    # Construct and run fastp command
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


def run_star(sample, file1, file2, args):
    """Run STAR alignment for a single sample."""
    print(f"\n--- Running STAR for sample {sample} ---")
    
    # Create temporary uncompressed files
    temp_file1 = os.path.join(args.fastp_dir, f"{sample}_1.fastq")
    temp_file2 = os.path.join(args.fastp_dir, f"{sample}_2.fastq")
    
    # Decompress files
    with open(temp_file1, 'wb') as outfile:
        subprocess.run(['gunzip', '-c', file1], stdout=outfile, check=True)
    
    with open(temp_file2, 'wb') as outfile:
        subprocess.run(['gunzip', '-c', file2], stdout=outfile, check=True)
    
    # Construct and run STAR command
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
    
    # Clean up temporary uncompressed files
    os.remove(temp_file1)
    os.remove(temp_file2)
    
    # Return paths to unmapped reads
    unmapped1 = os.path.join(args.star_dir, f"{sample}_Unmapped.out.mate1")
    unmapped2 = os.path.join(args.star_dir, f"{sample}_Unmapped.out.mate2")
    
    # Compress unmapped reads
    unmapped1_gz = f"{unmapped1}.gz"
    unmapped2_gz = f"{unmapped2}.gz"
    
    with open(unmapped1, 'rb') as infile, open(unmapped1_gz, 'wb') as outfile:
        subprocess.run(['gzip', '-c'], stdin=infile, stdout=outfile, check=True)
    
    with open(unmapped2, 'rb') as infile, open(unmapped2_gz, 'wb') as outfile:
        subprocess.run(['gzip', '-c'], stdin=infile, stdout=outfile, check=True)
    
    # Remove uncompressed unmapped files
    os.remove(unmapped1)
    os.remove(unmapped2)
    
    return unmapped1_gz, unmapped2_gz


def run_megahit(sample, file1, file2, args):
    """Run megahit assembly for a single sample."""
    print(f"\n--- Running megahit for sample {sample} ---")
    
    # Prepare output directory
    output_dir = os.path.join(args.megahit_dir, f"{sample}_assembly")
    
    # Remove output directory if it exists (megahit requires this)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    # Construct and run megahit command
    cmd = [
        'megahit',
        '-t', str(args.threads),
        '-m', str(args.memory),
        '-1', file1,
        '-2', file2,
        '-o', output_dir
    ]
    
    print(f"Assembling {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return output_dir


def extract_contigs(sample, assembly_dir, args):
    """Extract contigs for a single sample."""
    print(f"\n--- Extracting contigs for sample {sample} ---")
    
    # Path to contigs file
    contigs_file = os.path.join(assembly_dir, 'final.contigs.fa')
    
    # Check if contigs file exists
    if not os.path.isfile(contigs_file):
        print(f"Warning: No contigs file found for {sample}, skipping...")
        return None
    
    # Copy contigs file to output directory
    output_file = os.path.join(args.contigs_dir, f"{sample}.fasta")
    shutil.copy2(contigs_file, output_file)
    print(f"Copied contigs for {sample} to {output_file}")
    
    return output_file


def clean_temp_files(sample, args):
    """Clean up temporary files for a sample."""
    print(f"\n--- Cleaning temporary files for sample {sample} ---")
    
    # Remove fastp output files
    fastp_files = [
        os.path.join(args.fastp_dir, f"{sample}_1.fastq.gz"),
        os.path.join(args.fastp_dir, f"{sample}_2.fastq.gz")
    ]
    
    # Remove STAR output files
    star_pattern = os.path.join(args.star_dir, f"{sample}_*")
    star_files = glob.glob(star_pattern)
    
    # Remove megahit output directory
    megahit_dir = os.path.join(args.megahit_dir, f"{sample}_assembly")
    
    # Delete fastp files
    for file in fastp_files:
        if os.path.exists(file):
            os.remove(file)
            print(f"Removed {file}")
    
    # Delete STAR files
    for file in star_files:
        if os.path.exists(file):
            os.remove(file)
            print(f"Removed {file}")
    
    # Delete megahit directory
    if os.path.exists(megahit_dir):
        shutil.rmtree(megahit_dir)
        print(f"Removed directory {megahit_dir}")


def process_sample(sample, file1, file2, args):
    """Process a single sample through the entire pipeline."""
    print(f"\n=== Processing sample {sample} ===")
    
    fastp_out1 = None
    fastp_out2 = None
    star_out1 = None
    star_out2 = None
    megahit_dir = None
    contig_file = None
    
    try:
        # Step 1: Fastp trimming
        if not args.skip_fastp:
            fastp_out1, fastp_out2 = run_fastp(sample, file1, file2, args)
        else:
            print(f"Skipping Fastp trimming for {sample}")
            fastp_out1, fastp_out2 = file1, file2
        
        # Step 2: STAR alignment
        if not args.skip_star:
            star_out1, star_out2 = run_star(sample, fastp_out1, fastp_out2, args)
            input_to_megahit1, input_to_megahit2 = star_out1, star_out2
        else:
            print(f"Skipping STAR alignment for {sample}")
            input_to_megahit1, input_to_megahit2 = fastp_out1, fastp_out2
        
        # Step 3: megahit assembly
        if not args.skip_megahit:
            megahit_dir = run_megahit(sample, input_to_megahit1, input_to_megahit2, args)
        else:
            print(f"Skipping megahit assembly for {sample}")
        
        # Step 4: Extract contigs
        if not args.skip_contigs and megahit_dir is not None:
            contig_file = extract_contigs(sample, megahit_dir, args)
        else:
            print(f"Skipping contig extraction for {sample}")
        
        # Clean up temporary files if requested
        if not args.keep_temp_files:
            clean_temp_files(sample, args)
        
        print(f"=== Completed processing sample {sample} ===")
        return True
    
    except Exception as e:
        print(f"Error processing sample {sample}: {str(e)}")
        return False


def main():
    """Main function to run the bioinformatics pipeline."""
    args = parse_arguments()
    
    # Ensure all required directories exist
    ensure_directories([
        args.fastq_dir,
        args.fastp_dir,
        args.star_dir,
        args.megahit_dir,
        args.contigs_dir,
    ])
    
    # Find all _1.fastq.gz files in the input directory
    pattern = os.path.join(args.fastq_dir, '*_1.fastq.gz')
    fastq_files = glob.glob(pattern)
    
    if not fastq_files:
        print(f"No FASTQ files found in {args.fastq_dir}")
        sys.exit(1)
    
    print(f"Found {len(fastq_files)} samples to process")
    
    successful_samples = 0
    failed_samples = 0
    
    # Process each sample
    for file1 in fastq_files:
        # Extract sample name without the _1.fastq.gz suffix
        basename = os.path.basename(file1)
        sample = basename.replace('_1.fastq.gz', '')
        
        # Check if corresponding _2.fastq.gz file exists
        file2 = os.path.join(args.fastq_dir, f"{sample}_2.fastq.gz")
        if not os.path.isfile(file2):
            print(f"Warning: No matching pair found for {file1}, skipping...")
            continue
        
        # Process this sample
        success = process_sample(sample, file1, file2, args)
        
        if success:
            successful_samples += 1
        else:
            failed_samples += 1
    
    print(f"\n=== Pipeline Summary ===")
    print(f"Total samples: {len(fastq_files)}")
    print(f"Successfully processed: {successful_samples}")
    print(f"Failed: {failed_samples}")
    print("\n=== Bioinformatics Pipeline Completed ===")


if __name__ == "__main__":
    main()
