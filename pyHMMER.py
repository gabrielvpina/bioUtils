import pyhmmer
import collections

# Definir a estrutura para armazenar os resultados
Result = collections.namedtuple("Result", ["query", "subjID", "bitscore", "start", "end"])

# Caminho para os arquivos
hmm_file = "path/to/model.hmm"  # Modelo HMM
file = "path/to/seq.fasta"  # Arquivo FASTA de proteínas
output_tsv = "results.tsv"  # Arquivo de saída

# Criar um alfabeto para proteínas
alphabet = pyhmmer.easel.Alphabet.amino()

# Carregar sequências do arquivo FASTA corretamente
with pyhmmer.easel.SequenceFile(file, digital=True, alphabet=alphabet) as seq_file:
    sequences = list(seq_file)

results = []
with pyhmmer.plan7.HMMFile(hmm_file) as hmm_file:
    for hits in pyhmmer.hmmsearch(hmm_file, sequences, cpus=6):
        hmm_name = hits.query.name.decode()  # Nome do HMM

        for hit in hits:
            if hit.included:
                for domain in hit.domains:
                    results.append(Result(hit.name.decode(), hmm_name, hit.score, domain.env_from, domain.env_to))

# Escrevendo os resultados em TSV
with open(output_tsv, "w") as out:
    out.write("query\ttarget\tbitscore\tstart\tend\n")  # Cabeçalho
    for result in results:
        out.write(f"{result.query}\t{result.subjID}\t{result.bitscore}\t{result.start}\t{result.end}\n")

print(f"Resultados salvos em {output_tsv}")
