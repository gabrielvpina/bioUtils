import pandas as pd
import sys
import json

if len(sys.argv) != 3:
    print("Uso: python converter.py arquivo.tsv arquivo.json")
    sys.exit(1)

tsv_file = sys.argv[1]
json_file = sys.argv[2]

df = pd.read_csv(tsv_file, sep="\t")

with open(json_file, 'w', encoding='utf-8') as f:
    df.to_json(f, orient="records", indent=4, force_ascii=False)

print(f"Arquivo {json_file} gerado com sucesso!")
