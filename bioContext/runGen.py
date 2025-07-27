import os
import subprocess

# --- Configurações ---
main_bed_file = "dumyat_context.bed"
individual_beds_dir = "beds"
genomes_dir = "genomes"
gencon_script = "genCon.py" # Nome (ou caminho) do seu script genCon.py
context_value = "30000"

# --- 1. Criar diretório para os BEDs individuais, se não existir ---
os.makedirs(individual_beds_dir, exist_ok=True)
print(f"Diretório '{individual_beds_dir}' assegurado.")

# --- 2. Dividir o arquivo BED principal e coletar nomes dos arquivos gerados ---
generated_bed_files = []
try:
    with open(main_bed_file, "r") as infile:
        for line in infile:
            line = line.strip()
            if not line:  # Pular linhas vazias
                continue

            fields = line.split("\t")
            if len(fields) < 4:
                print(f"Aviso: Linha ignorada por ter menos de 4 campos: '{line}'")
                continue

            scaffold_name = fields[3]
            output_bed_filename = f"{scaffold_name}.bed"
            output_bed_filepath = os.path.join(individual_beds_dir, output_bed_filename)

            with open(output_bed_filepath, "w") as outfile_bed:
                outfile_bed.write(line + "\n") # Adiciona a quebra de linha original de volta
            
            generated_bed_files.append(output_bed_filepath)
            print(f"Gerado arquivo BED individual: '{output_bed_filepath}'")

except FileNotFoundError:
    print(f"Erro: Arquivo BED principal '{main_bed_file}' não encontrado.")
    exit()
except Exception as e:
    print(f"Erro ao processar o arquivo BED principal: {e}")
    exit()

if not generated_bed_files:
    print("Nenhum arquivo BED individual foi gerado. Verifique o arquivo de entrada.")
    exit()

print(f"\n--- {len(generated_bed_files)} arquivos BED individuais gerados em '{individual_beds_dir}' ---")

# --- 3. Listar arquivos de genoma ---
try:
    genome_files_all = [f for f in os.listdir(genomes_dir) if f.endswith(".fna")]
    if not genome_files_all:
        print(f"Nenhum arquivo .fna encontrado no diretório de genomas: '{genomes_dir}'")
        exit()
except FileNotFoundError:
    print(f"Erro: Diretório de genomas '{genomes_dir}' não encontrado.")
    exit()
except Exception as e:
    print(f"Erro ao listar arquivos de genoma: {e}")
    exit()

print(f"\n--- Encontrados {len(genome_files_all)} arquivos de genoma em '{genomes_dir}' ---")

# --- 4. Encontrar correspondências e executar genCon.py ---
print("\n--- Iniciando busca por correspondências e execução do genCon.py ---")
for bed_filepath in generated_bed_files:
    bed_filename_base = os.path.basename(bed_filepath) # Ex: FLAVIFRONS_EVE6.bed
    
    # Extrai o prefixo do nome do arquivo BED (antes do primeiro '_')
    # Se não houver '_', usa o nome do arquivo sem a extensão .bed
    bed_name_parts = os.path.splitext(bed_filename_base)[0].split("_", 1)
    bed_prefix = bed_name_parts[0]

    found_genome_match = False
    for genome_filename in genome_files_all: # Ex: AFFINIS_GCF_024516045.1.fna
        # Extrai o prefixo do nome do arquivo do genoma (antes do primeiro '_')
        # Se não houver '_', usa o nome do arquivo sem a extensão .fna
        genome_name_parts = os.path.splitext(genome_filename)[0].split("_", 1)
        genome_prefix = genome_name_parts[0]

        if bed_prefix == genome_prefix:
            print(f"\nCorrespondência encontrada: BED '{bed_filename_base}' (prefixo '{bed_prefix}') com Genoma '{genome_filename}' (prefixo '{genome_prefix}')")
            
            genome_filepath = os.path.join(genomes_dir, genome_filename)
            
            # Definir o nome do arquivo de saída para genCon.py
            output_gencon_name = f"{os.path.splitext(bed_filename_base)[0]}_30k"
            # Não incluir extensão no nome do arquivo de saída, genCon.py pode adicionar a sua própria

            # Montar o comando
            command = [
                "python", gencon_script,
                "-g", genome_filepath,
                "-b", bed_filepath,
                "--context", context_value,
                "-o", output_gencon_name
            ]

            print(f"Executando comando: {' '.join(command)}")

            try:
                # Executar o comando
                result = subprocess.run(command, check=True, capture_output=True, text=True)
                print(f"Saída do {gencon_script} para {bed_filename_base}:")
                if result.stdout:
                    print(f"STDOUT:\n{result.stdout}")
                if result.stderr: # genCon.py pode usar stderr para mensagens de progresso
                    print(f"STDERR:\n{result.stderr}")
                print(f"Comando executado com sucesso para '{bed_filename_base}'. Arquivo de saída esperado: '{output_gencon_name}' (mais a extensão que genCon.py adicionar).")
            
            except FileNotFoundError:
                print(f"Erro: Script '{gencon_script}' não encontrado. Verifique o nome e o caminho.")
                # Poderia decidir parar tudo aqui com exit() ou continuar para os próximos.
            except subprocess.CalledProcessError as e:
                print(f"Erro ao executar {gencon_script} para '{bed_filename_base}':")
                print(f"Comando: {' '.join(e.cmd)}")
                print(f"Código de retorno: {e.returncode}")
                print(f"STDOUT:\n{e.stdout}")
                print(f"STDERR:\n{e.stderr}")
            except Exception as e:
                print(f"Um erro inesperado ocorreu ao executar o subprocesso para {bed_filename_base}: {e}")

            found_genome_match = True
            break # Para de procurar por genomas para este arquivo BED, pois já encontrou um

    if not found_genome_match:
        print(f"\nAviso: Nenhum genoma correspondente encontrado para '{bed_filename_base}' (prefixo '{bed_prefix}')")

print("\n--- Processamento concluído ---")