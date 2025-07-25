import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dna_features_viewer import GraphicFeature, GraphicRecord
from BCBio import GFF
import subprocess
import argparse
import matplotlib.pyplot as plt

# parse arguments
parser = argparse.ArgumentParser(description='This script use genomic FASTA, GFF and BED files to create a genomic context of a region.')
parser.add_argument('--genome', '-g', required=True, help='Genomic FASTA file')
parser.add_argument('--annot', '-a', required=False, help='GFF file of the species.')
parser.add_argument('--bed', '-b', required=True, help='BED file with the regions of interest: chrom|start|end|name|0|strand')
parser.add_argument('--context', '-c', type=int, default=10000, help='Integer. Number in nt of flanking regions.')
parser.add_argument('--outdir', '-o', required=True, help='Output directory name.')
parser.add_argument('--merge_distance', '-m', type=int, default=None, 
                    help='Max distance between regions to merge them into a single context. If None, defaults to 2x context size.')

args = parser.parse_args()

# If merge_distance is not specified, set it to 2x the context size
if args.merge_distance is None:
    args.merge_distance = 2 * args.context

# input path
working_dir = "./"

# create output dir
new_directory_path = f'./{args.outdir}'
if not os.path.exists(new_directory_path):
    os.makedirs(new_directory_path)
else:
    print(f"Directory '{new_directory_path}' already exists.")

# Load genome sequences into memory for CDS extraction
print("Loading genome sequences...")
genome_sequences = {}
with open(args.genome, 'r') as genome_file:
    for record in SeqIO.parse(genome_file, 'fasta'):
        genome_sequences[record.id] = record.seq
print(f"Loaded {len(genome_sequences)} sequences from genome")

# open BED file and parse regions
regions_by_scaffold = {}
with open(os.path.join(working_dir, args.bed), "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        pos = line.split("\t")
        if len(pos) < 6:
            print(f"Warning: Skipping invalid BED line: {line}")
            continue
        
        scaffold = pos[0]
        start_coord = int(pos[1])
        end_coord = int(pos[2])
        name = pos[3]
        strand = pos[5]
        
        if scaffold not in regions_by_scaffold:
            regions_by_scaffold[scaffold] = []
        
        regions_by_scaffold[scaffold].append((start_coord, end_coord, name, strand))

# Function to merge overlapping or close-by regions
def merge_proximal_regions(regions, merge_distance):
    if not regions:
        return []
    
    # Sort regions by start position
    sorted_regions = sorted(regions, key=lambda x: x[0])
    
    merged_regions = []
    current_group = [sorted_regions[0]]
    
    for region in sorted_regions[1:]:
        current_start, current_end, _, _ = current_group[-1]
        next_start, next_end, _, _ = region
        
        # If this region is close enough to the previous one, add to current group
        if next_start - current_end <= merge_distance:
            current_group.append(region)
        else:
            # Process the current group and start a new one
            merged_regions.append(current_group)
            current_group = [region]
    
    # Don't forget the last group
    if current_group:
        merged_regions.append(current_group)
    
    return merged_regions

def extract_flanking_regions(gff3_file, scaffold, start, end, region_features, output_gff):
    with open(gff3_file, "r") as gff, open(output_gff, "w") as out_gff:
        for line in gff:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            feature_scaffold = fields[0]
            feature_start = int(fields[3])
            feature_end = int(fields[4])
            # verify the current flanking region
            if feature_scaffold == scaffold and ((feature_start >= start and feature_start <= end) or (feature_end >= start and feature_end <= end)):
                out_gff.write(line)

        # Mark all regions of interest in the output GFF
        for region_start, region_end, name, strand in region_features:
            strand_symbol = "+" if strand == "+" else "-"
            out_gff.write(
                f"{scaffold}\tContextScript\tRegionOfInterest\t{region_start}\t{region_end}\t.\t{strand_symbol}\t.\tID={name}_{region_start}_{region_end};Note=Region_of_Interest\n"
            )

def extract_cds_sequences(gff_file, scaffold, genome_seq, output_fasta):
    """Extract CDS sequences from the genome based on GFF annotations"""
    cds_sequences = []
    
    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith("#") or line.strip() == "":
                continue
            
            fields = line.split("\t")
            if len(fields) < 9:
                continue
                
            feature_scaffold = fields[0]
            feature_type = fields[2]
            feature_start = int(fields[3]) - 1  # Convert to 0-based indexing
            feature_end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            
            # Only process CDS features
            if feature_type.lower() == 'cds' and feature_scaffold == scaffold:
                # Extract sequence
                if scaffold in genome_sequences:
                    cds_seq = genome_sequences[scaffold][feature_start:feature_end]
                    
                    # Reverse complement if on negative strand
                    if strand == '-':
                        cds_seq = cds_seq.reverse_complement()
                    
                    # Parse attributes to get ID
                    attr_dict = {}
                    for attr in attributes.split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attr_dict[key] = value
                    
                    # Create sequence ID
                    seq_id = attr_dict.get('ID', f"CDS_{feature_start+1}_{feature_end}")
                    description = f"{scaffold}:{feature_start+1}-{feature_end}({strand})"
                    
                    # Create SeqRecord
                    seq_record = SeqRecord(
                        cds_seq,
                        id=seq_id,
                        description=description
                    )
                    cds_sequences.append(seq_record)
    
    # Write CDS sequences to FASTA file
    if cds_sequences:
        with open(output_fasta, 'w') as output_handle:
            SeqIO.write(cds_sequences, output_handle, 'fasta')
        print(f"Extracted {len(cds_sequences)} CDS sequences to {output_fasta}")
    else:
        print(f"No CDS sequences found in {gff_file}")

def parse_gff3_original(file_path):
    """read and parse GFF file"""
    features = []
    with open(file_path, 'r') as gff_file:
        for record in GFF.parse(gff_file):
            for feature in record.features:
                start = feature.location.start
                end = feature.location.end
                strand = feature.location.strand
                label = feature.qualifiers.get("ID", [feature.type])[0]
                # Add feature in the format expected by the library
                features.append(
                    GraphicFeature(start=start, end=end, strand=strand, label=label)
                )
    return features

def plot_gff3_files(input_directory, output_directory):
    """process gff files for each plot"""
    os.makedirs(output_directory, exist_ok=True)
    
    for file_name in os.listdir(input_directory):
        if file_name.endswith("context.gff"):
            output_path = os.path.join(output_directory, file_name.replace("_context.gff", "_context.pdf"))
            print(f"Processing {file_name}...")
            
            # Extract features from GFF3 file
            features = parse_gff3_original(os.path.join(input_directory, file_name))
            
            # Check if there are features to process
            if not features:
                print(f"Error: No features found in {file_name}.")
                continue
            
            # Find the minimum start position and adjust coordinates
            min_start = min(feature.start for feature in features)
            max_end = max(feature.end for feature in features)
            
            # Adjust feature coordinates to start from 0
            adjusted_features = []
            for feature in features:
                adjusted_feature = GraphicFeature(
                    start=feature.start - min_start,
                    end=feature.end - min_start,
                    strand=feature.strand,
                    label=feature.label
                )
                adjusted_features.append(adjusted_feature)
            
            # Create GraphicRecord with adjusted sequence length
            sequence_length = max_end - min_start
            graphic_record = GraphicRecord(sequence_length=sequence_length, features=adjusted_features)
            
            # Plot and save the image
            ax, *_ = graphic_record.plot(figure_width=10)
            
            # Add original coordinates to the plot title or labels if needed
            ax.set_title(f"Sequence region: {min_start}-{max_end}")
            
            ax.figure.savefig(output_path, format="pdf")
            print(f"Graph saved at {output_path}")
            
            # Close the figure to free memory
            plt.close(ax.figure)

def create_bed_for_region(region_group, scaffold, output_bed):
    """Create a BED file for a group of regions"""
    with open(output_bed, 'w') as f:
        for start, end, name, strand in region_group:
            f.write(f"{scaffold}\t{start}\t{end}\t{name}\t0\t{strand}\n")

# Process each scaffold separately
for scaffold, regions in regions_by_scaffold.items():
    # Group proximal regions
    merged_region_groups = merge_proximal_regions(regions, args.merge_distance)
    
    print(f"Processing {scaffold}: {len(regions)} regions grouped into {len(merged_region_groups)} contexts")
    
    # Process each group of regions
    for i, region_group in enumerate(merged_region_groups):
        # Create a meaningful group name
        if len(region_group) == 1:
            group_name = region_group[0][2]  # Use the single region's name
        else:
            # Use first and last region names for the group
            group_name = f"{region_group[0][2]}-to-{region_group[-1][2]}"
        
        # Get the overall start and end coordinates for the group
        group_start = min(r[0] for r in region_group)
        group_end = max(r[1] for r in region_group)
        
        # Create individual directory for this group
        ind_path = os.path.join(new_directory_path, group_name)
        os.makedirs(ind_path, exist_ok=True)
        
        # Define flanking regions for the entire group
        contextNumber = args.context
        flanking_start = max(0, group_start - contextNumber)  # avoid negative coordinates
        flanking_end = group_end + contextNumber
        
        # Create output files
        output_gff_file = os.path.join(ind_path, f"{group_name}_context.gff")
        
        if args.annot:
            # Extract flanking regions for all regions in the group
            extract_flanking_regions(args.annot, scaffold, flanking_start, flanking_end, region_group, output_gff_file)

            # Extract CDS sequences from the context GFF
            cds_fasta_file = os.path.join(ind_path, f"{group_name}_CDS_sequences.fasta")
            extract_cds_sequences(output_gff_file, scaffold, genome_sequences[scaffold] if scaffold in genome_sequences else None, cds_fasta_file)

            # Create plot
            plot_gff3_files(ind_path, ind_path)
   
        # Create temporary BED file for this group
        temp_bed_file = os.path.join(ind_path, f"{group_name}_regions.bed")
        create_bed_for_region(region_group, scaffold, temp_bed_file)
        
        # Extract sequence for the entire group context
        flanking_fasta = os.path.join(ind_path, f"{group_name}_flanking{str(contextNumber)}_{scaffold}.fasta")
        
        # Create a BED file with the full region including context
        full_region_bed = os.path.join(ind_path, f"{group_name}_fullregion.bed")
        with open(full_region_bed, 'w') as f:
            f.write(f"{scaffold}\t{flanking_start}\t{flanking_end}\t{group_name}\t0\t+\n")
        
        # Extract the flanking sequence using bedtools
        bedtools_command = [
            "bedtools", "getfasta",
            "-fi", args.genome,
            "-bed", full_region_bed,
            "-fo", flanking_fasta,
            "-name"
        ]
        
        try:
            subprocess.run(bedtools_command, check=True)
            print(f"Successfully extracted sequence for {group_name}")
        except Exception as e:
            print(f"Error in bedtools command for {group_name}: {e}")
            continue
        
        # Extract individual region sequences
        regions_fasta = os.path.join(ind_path, f"{group_name}_individual_regions_{scaffold}.fasta")
        
        bedtools_command = [
            "bedtools", "getfasta",
            "-fi", args.genome,
            "-bed", temp_bed_file,
            "-fo", regions_fasta,
            "-name"
        ]
        
        try:
            subprocess.run(bedtools_command, check=True)
        except Exception as e:
            print(f"Error extracting individual regions for {group_name}: {e}")

print("Analysis completed!")
