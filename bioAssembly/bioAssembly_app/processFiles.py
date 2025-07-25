import os, shutil, glob


def extract_contigs_megahit(sample, assembly_dir, args):
    """Extract contigs for a single sample assembled with MEGAHIT."""
    print(f"\n--- Extracting contigs for sample {sample} (MEGAHIT) ---")

    # MEGAHIT's final contigs file is typically named 'final.contigs.fa'
    contigs_file = os.path.join(assembly_dir, 'final.contigs.fa')

    # Check if contigs file exists
    if not os.path.isfile(contigs_file):
        print(f"Warning: No MEGAHIT contigs file found for {sample} at {contigs_file}, skipping...")
        return None

    # Copy contigs file to output dir
    output_file = os.path.join(args.contigs_dir, f"{sample}_megahit.fasta")
    shutil.copy2(contigs_file, output_file)
    print(f"Copied MEGAHIT contigs for {sample} to {output_file}")

    return output_file


######################################################
def extract_contigs_spades(sample, assembly_dir, args):
    """Extract contigs for a single sample assembled with SPAdes (including metaSPAdes)."""
    print(f"\n--- Extracting contigs for sample {sample} (SPAdes) ---")

    # SPAdes (including metaSPAdes and regular spades.py) usually outputs contigs.fasta
    # or contigs.fasta in the 'scaffolds.fasta' file (which often includes contigs too)
    # The primary output for contigs is usually 'contigs.fasta'
    contigs_file = os.path.join(assembly_dir, 'contigs.fasta')

    # Check if contigs file exists
    if not os.path.isfile(contigs_file):
        print(f"Warning: No SPAdes contigs file found for {sample} at {contigs_file}, skipping...")
        return None

    # Copy contigs file to output dir
    output_file = os.path.join(args.contigs_dir, f"{sample}_spades.fasta")
    shutil.copy2(contigs_file, output_file)
    print(f"Copied SPAdes contigs for {sample} to {output_file}")

    return output_file


######################################################
def extract_contigs_rnaspades(sample, assembly_dir, args):
    """Extract contigs for a single sample assembled with RNAspades."""
    print(f"\n--- Extracting contigs for sample {sample} (RNAspades) ---")

    # RNAspades also typically outputs 'contigs.fasta' similar to regular SPAdes
    contigs_file = os.path.join(assembly_dir, 'contigs.fasta')

    # Check if contigs file exists
    if not os.path.isfile(contigs_file):
        print(f"Warning: No RNAspades contigs file found for {sample} at {contigs_file}, skipping...")
        return None

    # Copy contigs file to output dir
    output_file = os.path.join(args.contigs_dir, f"{sample}_rnaspades.fasta")
    shutil.copy2(contigs_file, output_file)
    print(f"Copied RNAspades contigs for {sample} to {output_file}")

    return output_file


######################################################
def extract_contigs_trinity(sample, assembly_dir, args):
    """Extract contigs (transcripts) for a single sample assembled with Trinity."""
    print(f"\n--- Extracting contigs for sample {sample} (Trinity) ---")

    # Trinity's main output file containing assembled transcripts (contigs)
    # is usually named 'Trinity.fasta' and located directly in the assembly_dir
    contigs_file = os.path.join(assembly_dir, 'Trinity.fasta')

    # Check if contigs file exists
    if not os.path.isfile(contigs_file):
        print(f"Warning: No Trinity transcripts file found for {sample} at {contigs_file}, skipping...")
        return None

    # Copy contigs file to output dir
    output_file = os.path.join(args.contigs_dir, f"{sample}_trinity.fasta")
    shutil.copy2(contigs_file, output_file)
    print(f"Copied Trinity transcripts for {sample} to {output_file}")

    return output_file







def clean_temp_files(sample, args):
    """Clean up temporary files for a sample."""
    print(f"\n--- Cleaning temporary files for sample {sample} ---")


# Set trimmed files
#######################################################
    # Remove fastp output files
    fastp_files = [
        os.path.join(args.fastp_dir, f"{sample}_1.fastq.gz"),
        os.path.join(args.fastp_dir, f"{sample}_2.fastq.gz")
    ]
    
    # Remove trimmomatic output files
    trimmomatic_files = [
        os.path.join(args.trimmomatic_dir, f"{sample}_1.fastq.gz"),
        os.path.join(args.trimmomatic_dir, f"{sample}_2.fastq.gz")
    ]
#######################################################



    # Remove STAR output files
    star_pattern = os.path.join(args.star_dir, f"{sample}_*")
    star_files = glob.glob(star_pattern)
    
    # Remove assembly output directory
    megahit_dir = os.path.join(args.megahit_dir, f"{sample}_assembly")
    spades_dir = os.path.join(args.spades_dir, f"{sample}_assembly")
    rnaspades_dir = os.path.join(args.rnaspades_dir, f"{sample}_assembly")
    trinity_dir = os.path.join(args.trinity_dir, f"{sample}_assembly")
    

# Delete trimmed files
#######################################################
    # Delete fastp files
    for file in fastp_files:
        if os.path.exists(file):
            os.remove(file)
            print(f"Removed {file}")

    # Delete trimmomatic files
    for file in trimmomatic_files:
        if os.path.exists(file):
            os.remove(file)
            print(f"Removed {file}")
#######################################################

    
    # Delete STAR files
    for file in star_files:
        if os.path.exists(file):
            os.remove(file)
            print(f"Removed {file}")
    
    # Delete megahit directory
    if os.path.exists(megahit_dir):
        shutil.rmtree(megahit_dir)
        print(f"Removed directory {megahit_dir}")