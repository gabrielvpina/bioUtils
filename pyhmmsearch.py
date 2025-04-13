#!/usr/bin/env python3
"""
PyHMMER Search Tool - Search protein sequences against HMM profiles

This script performs HMMER searches using pyhmmer to identify domain matches
in protein sequence files based on HMM profiles.
"""

import os
import sys
import argparse
import collections
import pyhmmer
from pyhmmer.easel import SequenceFile, Alphabet
from pyhmmer.plan7 import HMMFile


# Define a structure to store results
Result = collections.namedtuple("Result", ["query", "target", "bitscore", "start", "end"])


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Search protein sequences against HMM profiles using pyhmmer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "-m", "--hmm", 
        required=True,
        help="Path to the HMM profile file"
    )
    
    parser.add_argument(
        "-i", "--input", 
        required=True,
        help="Path to the input FASTA file containing protein sequences"
    )
    
    parser.add_argument(
        "-o", "--output", 
        default="hmmer_results.tsv",
        help="Path to the output TSV file"
    )
    
    parser.add_argument(
        "-c", "--cpus", 
        type=int, 
        default=1,
        help="Number of CPU cores to use"
    )
    
    parser.add_argument(
        "-e", "--evalue", 
        type=float, 
        default=0.01,
        help="E-value threshold for filtering results"
    )
    
    parser.add_argument(
        "--verbose", 
        action="store_true",
        help="Enable verbose output"
    )
    
    return parser.parse_args()


def run_hmmsearch(hmm_file_path, sequence_file_path, cpus=1, e_value=0.01, verbose=False):
    """
    Run HMMER search on input sequences using the specified HMM profile.
    
    Args:
        hmm_file_path: Path to the HMM profile file
        sequence_file_path: Path to the FASTA file with protein sequences
        cpus: Number of CPU cores to use
        e_value: E-value threshold for filtering results
        verbose: Enable verbose output
        
    Returns:
        List of Result namedtuples
    """
    results = []
    
    # Create a protein alphabet
    alphabet = Alphabet.amino()
    
    if verbose:
        print(f"Loading sequences from {sequence_file_path}")
    
    try:
        # Load protein sequences from FASTA file
        with SequenceFile(sequence_file_path, digital=True, alphabet=alphabet) as seq_file:
            sequences = list(seq_file)
            
        if verbose:
            print(f"Loaded {len(sequences)} sequences")
            print(f"Running HMMER search with {cpus} CPUs")
        
        # Run HMMER search
        with HMMFile(hmm_file_path) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, sequences, cpus=cpus, E=e_value):
                hmm_name = hits.query.name.decode()  # Name of the HMM
                
                for hit in hits:
                    if hit.included:
                        seq_name = hit.name.decode()
                        
                        for domain in hit.domains:
                            results.append(Result(
                                seq_name,  # Query sequence name
                                hmm_name,  # Target HMM name
                                hit.score,  # Bitscore
                                domain.env_from,  # Start position
                                domain.env_to  # End position
                            ))
        
        if verbose:
            print(f"Found {len(results)} matches")
            
        return results
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def write_results(results, output_path, verbose=False):
    """
    Write search results to a TSV file.
    
    Args:
        results: List of Result namedtuples
        output_path: Path to the output TSV file
        verbose: Enable verbose output
    """
    try:
        with open(output_path, "w") as out:
            # Write header
            out.write("query\ttarget\tbitscore\tstart\tend\n")
            
            # Write results
            for result in results:
                out.write(f"{result.query}\t{result.target}\t{result.bitscore}\t{result.start}\t{result.end}\n")
        
        if verbose:
            print(f"Results saved to {output_path}")
            
    except IOError as e:
        print(f"Error writing to output file: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    """Main function to run the HMMER search pipeline."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Check if input files exist
    if not os.path.isfile(args.hmm):
        print(f"Error: HMM file not found: {args.hmm}", file=sys.stderr)
        sys.exit(1)
        
    if not os.path.isfile(args.input):
        print(f"Error: Input sequence file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    if args.verbose:
        print("Starting PyHMMER search...")
        print(f"HMM file: {args.hmm}")
        print(f"Input sequence file: {args.input}")
        print(f"Output file: {args.output}")
        print(f"Using {args.cpus} CPU cores")
        print(f"E-value threshold: {args.evalue}")
    
    # Run HMMER search
    results = run_hmmsearch(
        args.hmm, 
        args.input, 
        cpus=args.cpus, 
        e_value=args.evalue,
        verbose=args.verbose
    )
    
    # Write results to file
    write_results(results, args.output, verbose=args.verbose)
    
    if args.verbose:
        print("Search completed successfully")


if __name__ == "__main__":
    main()
