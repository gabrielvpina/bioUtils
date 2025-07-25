import os
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord
from BCBio import GFF
import argparse

# parse arguments
parser = argparse.ArgumentParser(description='Display and save genomic regions of GFF files')
parser.add_argument('--gff', '-g', required=True, help='gff file')

args = parser.parse_args()

def parse_gff3(file_path):
    """Reads a GFF3 file and extracts features with focus on exons."""
    features = []
    
    with open(file_path, 'r') as gff_file:
        for record in GFF.parse(gff_file):
            for feature in record.features:
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.strand
                
                # feature type and ID
                feature_type = feature.type
                feature_id = feature.qualifiers.get("ID", [feature_type])[0]
                
                # color based on feature type
                color = "#cccccc" 
                

                if feature_type.lower() == "exon":
                    color = "#ff9900"  # Orange
                    label = f"Exon:{feature_id}"
                elif feature_type.lower() == "cds":
                    color = "#00cc66"  # Green 
                    label = f"CDS:{feature_id}"
                elif feature_type.lower() == "gene":
                    color = "#6600cc"  # Purple 
                    label = f"Gene:{feature_id}"
                else:
                    label = f"{feature_type}:{feature_id}"
                
                # Add the feature in the format expected by the library
                features.append(
                    GraphicFeature(
                        start=start, 
                        end=end, 
                        strand=strand, 
                        label=label,
                        color=color
                    )
                )
    
    return features

def plot_gff3_file(input_file, output_file=None):
    """Plots features from a GFF3 file with a focus on exons."""
    # features from GFF3 file
    features = parse_gff3(input_file)
    
    if not features:
        print(f"No features found in {input_file}")
        return
    
    # GraphicRecord
    sequence_length = max(feature.end for feature in features)
    graphic_record = GraphicRecord(sequence_length=sequence_length, features=features)
    
    # figure with appropriate size
    plt.figure(figsize=(12, 5))
    
    # graphic record
    ax, _ = graphic_record.plot(figure_width=10)
    

    plt.title(f"Features visualization from {os.path.basename(input_file)}")
    
    # save the image if output file is provided
    if output_file:
        output_format = output_file.split('.')[-1] if '.' in output_file else 'pdf'
        plt.savefig(output_file, format=output_format, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        output_file = os.path.splitext(input_file)[0] + '_visualization.pdf'
        plt.savefig(output_file, format='pdf', bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    
    # show  plot
    plt.show()

if __name__ == "__main__":
  
    input_file = args.gff

    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found")
    else:
        # output filename
        output_file = os.path.splitext(input_file)[0] + "_visualization.pdf"
        
        # process and plot
        plot_gff3_file(input_file, output_file)
