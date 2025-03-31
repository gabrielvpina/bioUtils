import pyhmmer
import collections

'''
make a multiprocessing hmmsearch using pyHMMER
'''

Result = collections.namedtuple("Result", ["query", "subjID", "bitscore", "start", "end"])

hmm_file = "path/to/model.hmm"  
file = "path/to/seq.fasta"  
output_tsv = "results.tsv"  

alphabet = pyhmmer.easel.Alphabet.amino()

with pyhmmer.easel.SequenceFile(file, digital=True, alphabet=alphabet) as seq_file:
    sequences = list(seq_file)

results = []
with pyhmmer.plan7.HMMFile(hmm_file) as hmm_file:
    for hits in pyhmmer.hmmsearch(hmm_file, sequences, cpus=6):
        hmm_name = hits.query.name.decode() 

        for hit in hits:
            if hit.included:
                for domain in hit.domains:
                    results.append(Result(hit.name.decode(), hmm_name, hit.score, domain.env_from, domain.env_to))

# result in tsv
with open(output_tsv, "w") as out:
    out.write("query\ttarget\tbitscore\tstart\tend\n") 
    for result in results:
        out.write(f"{result.query}\t{result.subjID}\t{result.bitscore}\t{result.start}\t{result.end}\n")

print(f"Resultados salvos em {output_tsv}")
