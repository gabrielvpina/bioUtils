import pandas as pd
import sys

# Verifica se os argumentos foram passados corretamente
if len(sys.argv) < 3:
    print("Uso: python script.py arquivo.csv arquivo.json")
    sys.exit(1)

# Nome dos arquivos
csv_file = sys.argv[1]
json_file = sys.argv[2]

# Ler o CSV
df = pd.read_csv(csv_file)

# Converter para JSON e salvar
df.to_json(json_file, orient='records', indent=4)

print(f"Conversão concluída! JSON salvo em {json_file}")
