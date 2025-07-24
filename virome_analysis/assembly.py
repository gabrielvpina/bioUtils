import os, subprocess


def run_megahit_PE(sample, file1, file2, args):
    """Run megahit assembly for a PAIR-END sample."""
    print(f"\n--- Running megahit for sample {sample} ---")
    
    # Prepare output directory
    output_dir = os.path.join(args.megahit_dir, f"{sample}_assembly")
    
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