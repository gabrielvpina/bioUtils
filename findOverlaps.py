from Bio import SeqIO

"""
Find alignments in 3' and 5' of sequences in a multifasta file.
"""

def check_overlaps(fasta_file, min_overlap=10):
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    overlaps = []

    for id1, seq1 in sequences.items():
        for id2, seq2 in sequences.items():
            if id1 != id2:
                # Check overlap from end of seq1 to start of seq2
                for overlap_len in range(min_overlap, min(len(seq1), len(seq2)) + 1):
                    if seq1[-overlap_len:] == seq2[:overlap_len]:
                        overlaps.append((id1, "end", id2, "start", overlap_len, seq1[-overlap_len:]))

                # Check overlap from start of seq1 to end of seq2
                for overlap_len in range(min_overlap, min(len(seq1), len(seq2)) + 1):
                    if seq1[:overlap_len] == seq2[-overlap_len:]:
                        overlaps.append((id1, "start", id2, "end", overlap_len, seq1[:overlap_len]))

    return overlaps

# Example usage
fasta_file = "allViralContigs_CAP3.fasta"
overlap_results = check_overlaps(fasta_file)

for overlap in overlap_results:
    print(f"Overlap found: {overlap[0]} ({overlap[1]}) aligns with {overlap[2]} ({overlap[3]}) | "
          f"Length: {overlap[4]} | Sequence: {overlap[5]}")
