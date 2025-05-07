import os
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
from BCBio import GFF
import subprocess
import argparse

# parse arguments
parser = argparse.ArgumentParser(description='This script use genomic FASTA, GFF and BED files to create a genomic context of a region.')
parser.add_argument('--genome', '-g', required=True, help='Genomic FASTA file')
parser.add_argument('--annot', '-a', required=True, help='GFF file of the species.')
parser.add_argument('--bed', '-b', required=True, help='BED file with the regions of interest: chrom|start|end|name|0|strand')
parser.add_argument('--context', '-c', type=int, default=10000, help='Integer. Number in nt of flanking regions.')
parser.add_argument('--outdir', '-o', required=True, help='Output directory name.')

args = parser.parse_args()

# input path
working_dir = "./"

# create output dir
new_directory_path = f'./{args.outdir}'
if not os.path.exists(new_directory_path):
    os.makedirs(new_directory_path)
else:
    print(f"Directory '{new_directory_path}' already exists.")


# open BED file
with open(os.path.join(working_dir, args.bed), "r") as f:
    regions = f.readlines()

def extract_flanking_regions(gff3_file, scaffold, start, end, region_start, region_end, output_gff):
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

        out_gff.write(
            f"{scaffold}\tContextScript\tRegionOfInterest\t{region_start}\t{region_end}\t.\t+\t.\tID={name}_{region_start}_{region_end};Note=Region_of_Interest\n"
        )


# process BED files
for region in regions:

    pos = region.strip().split("\t")
    if len(pos) < 5:
        break
    scaffold = pos[0]
    name = pos[3]
    strand = pos[5]
    start_coord = int(pos[1])
    end_coord = int(pos[2])


    # individual dir
    ind_path = os.path.join(new_directory_path, name)
    os.makedirs(ind_path)

    contextNumber = args.context
    # define flanking regions
    flanking_start = max(0, start_coord - contextNumber)  # avoid negative coordinates
    flanking_end = end_coord + contextNumber

    # out file 
    output_gff_file = os.path.join(ind_path, f"{name}_context.gff")

    output_annot = os.path.join(ind_path, f"{name}_flanking{str(contextNumber)}_{scaffold}.gff")
    extract_flanking_regions(args.annot, scaffold, flanking_start, flanking_end, start_coord, end_coord, output_gff_file)

    ############## ADJUST GFF ########################################
    def adjust_gff_coordinates(input_file, output_file):
        """adjust the gff file to start from 0"""
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            min_start = None

            # find minor value
            for line in infile:
                if line.startswith("#") or line.strip() == "":
                    continue
                parts = line.split("\t")
                start = int(parts[3])
                if min_start is None or start < min_start:
                    min_start = start

            infile.seek(0)  # back to top

            # adjust coordinates based on minor value
            for line in infile:
                if line.startswith("#") or line.strip() == "":
                    outfile.write(line)
                    continue

                parts = line.split("\t")
                parts[3] = str(int(parts[3]) - min_start)  # adjust start
                parts[4] = str(int(parts[4]) - min_start)  # adjust end
                outfile.write("\t".join(parts))  # write adjusted line

    inputFile = output_gff_file
    outputFile = os.path.join(ind_path, f"{name}_contextNorm.gff")

    adjust_gff_coordinates(inputFile, outputFile)
    ########################################################################


    ################ CREATE PLOT #####################################
    def parse_gff3(file_path):
        """read and parse GFF file"""
        features = []
        with open(file_path, 'r') as gff_file:
            for record in GFF.parse(gff_file):
                for feature in record.features:
                    start = feature.location.start
                    end = feature.location.end
                    strand = feature.location.strand
                    label = feature.qualifiers.get("ID", [feature.type])[0]

                    # Add featrures
                    features.append(
                        GraphicFeature(start=start, end=end, strand=strand, label=label)
                    )
        return features

    def plot_gff3_files(input_directory, output_directory):
        """process gff files for each plot"""
        os.makedirs(output_directory, exist_ok=True)

        for file_name in os.listdir(input_directory):
            if file_name.endswith("Norm.gff"):
                input_path = os.path.join(input_directory, file_name)
                output_path = os.path.join(output_directory, f"{file_name}.pdf")
                
                print(f"Processing {file_name}...")
                
                # extract features from gff file
                features = parse_gff3(outputFile)
                
                
                if not features:
                    print(f"Error: None feature found in {file_name}.")
                    continue
                
                # create GraphicRecord
                sequence_length = max(feature.end for feature in features)
                graphic_record = GraphicRecord(sequence_length=sequence_length, features=features)
                
                # plot and save
                ax, _ = graphic_record.plot(figure_width=10)
                ax.figure.savefig(output_path, format="pdf")
                print(f"Gráfico salvo em {output_path}")

    plot_gff3_files(ind_path, ind_path)

    ##################################################################

    # genomic FASTA sequence
    genome_fasta = args.genome
    flanking_fasta = os.path.join(ind_path, f"{name}_flanking{str(contextNumber)}_{scaffold}.fasta")

    # bedtools getfasta
    bedtools_command = [
        "bedtools", "getfasta",
        "-fi", genome_fasta,
        "-bed", args.bed,
        "-fo", flanking_fasta,
        "-name"
    ]

    try:
        subprocess.run(bedtools_command, check=True)
    except Exception as e:
        print(f"Error in bedtools command: {e}")
        continue

    # extract seq from flanking interval
    complete_flanking_sequence = ""
    for record in SeqIO.parse(genome_fasta, "fasta"):
        if record.id == scaffold:
            complete_flanking_sequence = str(record.seq[flanking_start:flanking_end])
            break

    # combine completes seqs from flanking interval
    combined_sequence = (
        complete_flanking_sequence[:contextNumber]  # Região inicial                
        + complete_flanking_sequence[-contextNumber:]  # Região final
    )

    # Criar arquivo final com a sequência combinada
    combined_fasta = os.path.join(ind_path, f"{name}_ContextRegions_{scaffold}.fasta")
    with open(combined_fasta, "w") as combined_out:
        combined_out.write(f">Combined_{name}_flanking_{scaffold}\n{combined_sequence}\n")

print("Analysis completed!")







