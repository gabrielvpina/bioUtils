import os, subprocess


# MEGAHIT Assembly
# =============================================================================

# PAIR-END
def run_megahit_PE(sample, file1, file2, args):
    """Run megahit assembly for a PAIR-END sample."""
    print(f"\n--- Running megahit for sample {sample} ---")
    
    if not os.path.exists(args.assembly):
        os.makedirs(os.path.join(args.outdir, args.assembly))

    # Prepare output directory
    output_dir = os.path.join(args.outdir, args.assembly, f"{sample}_assembly")
    
    # Remove output directory if it exists (megahit requires this)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    # run megahit command
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



# Single-end
def run_megahit_SE(sample, file1, args):
    """Run Megahit assembly for a SINGLE-END sample."""
    print(f"\n--- Running megahit for sample {sample} ---")

    if not os.path.exists(args.assembly):
        os.makedirs(os.path.join(args.outdir, args.assembly))
    
    # Prepare output directory
    output_dir = os.path.join(args.outdir, args.assembly, f"{sample}_assembly")
    
    # Remove output directory if it exists (megahit requires this)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    # run megahit command
    cmd = [
        'megahit',
        '-t', str(args.threads),
        '-m', str(args.memory),
        '-r', file1,
        '-o', output_dir
    ]
    
    print(f"Assembling {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return output_dir

# =============================================================================





# SPAdes Assembly
# =============================================================================

# PAIR-END
def run_spades_PE(sample, file1, file2, args):
    """Run SPAdes assembly for a PAIR-END sample."""
    print(f"\n--- Running megahit for PAIR-END samples {sample} ---")

    if not os.path.exists(args.assembly):
        os.makedirs(os.path.join(args.outdir, args.assembly))
    
    # Prepare output directory
    output_dir = os.path.join(args.outdir, args.assembly, f"{sample}_assembly")
    
    # Remove output directory if it exists
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    # run spades command
    cmd = [
        'spades.py',
        '-t', str(args.threads),
        '-m', str(args.memory),
        '-1', file1,
        '-2', file2,
        '-o', output_dir
    ]
    
    print(f"Assembling {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return output_dir


# SINGLE-END
def run_spades_SE(sample, file1, args):
    """Run SPAdes assembly for a SINGLE-END sample."""
    print(f"\n--- Running megahit for SINGLE-END samples {sample} ---")

    if not os.path.exists(args.assembly):
        os.makedirs(os.path.join(args.outdir, args.assembly))
    
    # Prepare output directory
    output_dir = os.path.join(args.outdir, args.assembly, f"{sample}_assembly")
    
    # Remove output directory if it exists (megahit requires this)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    # run spades command
    cmd = [
        'spades.py',
        '-t', str(args.threads),
        '-m', str(args.memory),
        '-s', file1,
        '-o', output_dir
    ]
    
    print(f"Assembling {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return output_dir

# =============================================================================




# rnaspades Assembly
# =============================================================================

# PAIR-END
def run_rnaspades_PE(sample, file1, file2, args):
    """Run rnaspades assembly for a PAIR-END sample."""
    print(f"\n--- Running rnaspades.py for PAIR-END samples {sample} ---")

    if not os.path.exists(args.assembly):
        os.makedirs(os.path.join(args.outdir, args.assembly))
    
    # Prepare output directory
    output_dir = os.path.join(args.outdir, args.assembly, f"{sample}_assembly")
    
    # Remove output directory if it exists
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    # run rnaspades command
    cmd = [
        'rnaspades.py',
        '-t', str(args.threads),
        '-m', str(args.memory),
        '-1', file1,
        '-2', file2,
        '-o', output_dir
    ]
    
    print(f"Assembling {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return output_dir


# SINGLE-END
def run_rnaspades_PE(sample, file1, args):
    """Run rnaspades assembly for a SINGLE-END sample."""
    print(f"\n--- Running rnaspades.py for SINGLE-END samples {sample} ---")

    if not os.path.exists(args.assembly):
        os.makedirs(os.path.join(args.outdir, args.assembly))
    
    # Prepare output directory
    output_dir = os.path.join(args.outdir, args.assembly, f"{sample}_assembly")
    
    # Remove output directory if it exists
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    # run rnaspades command
    cmd = [
        'rnaspades.py',
        '-t', str(args.threads),
        '-m', str(args.memory),
        '-s', file1,
        '-o', output_dir
    ]
    
    print(f"Assembling {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return output_dir

# =============================================================================



# Trinity Assembly
# =============================================================================

# PAIR-END
def run_trinity_PE(sample, file1, file2, args):
    """Run trinity assembly for a PAIR-END sample."""
    print(f"\n--- Running trinity for PAIR-END samples {sample} ---")

    if not os.path.exists(args.assembly):
        os.makedirs(os.path.join(args.outdir, args.assembly))
    
    # Prepare output directory
    output_dir = os.path.join(args.outdir, args.assembly, f"{sample}_assembly")
    
    # Remove output directory if it exists
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    # run trinity command
    cmd = [
        'Trinity',
        '--CPU', str(args.threads),
        '--max_memory', str(args.memory),
        '--left', file1,
        '--right', file2,
        '--output', output_dir
    ]
    
    print(f"Assembling {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return output_dir


# SINGLE-END
def run_trinity_PE(sample, file1, args):
    """Run trinity assembly for a SINGLE-END sample."""
    print(f"\n--- Running trinity for SINGLE-END samples {sample} ---")

    if not os.path.exists(args.assembly):
        os.makedirs(os.path.join(args.outdir, args.assembly))
    
    # Prepare output directory
    output_dir = os.path.join(args.outdir, args.assembly, f"{sample}_assembly")
    
    # Remove output directory if it exists
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    # run trinity command
    cmd = [
        'Trinity',
        '--CPU', str(args.threads),
        '--max_memory', str(args.memory),
        '--single', file1,
        '--output', output_dir
    ]
    
    print(f"Assembling {sample}. Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return output_dir