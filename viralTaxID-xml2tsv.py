import xml.etree.ElementTree as ET
import pandas as pd

def parse_taxonomy(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    data = []
    
    for taxon in root.findall("Taxon"):
        tax_id = taxon.findtext("TaxId")
        sci_name = taxon.findtext("ScientificName")
        
        taxonomy = {"TaxId": tax_id, "ScientificName": sci_name}
        
        for lineage in taxon.findall("LineageEx/Taxon"):
            rank = lineage.findtext("Rank")
            name = lineage.findtext("ScientificName")
            
            if rank and name and rank != "superkingdom":
                taxonomy[rank] = name
        
        data.append(taxonomy)
    
    df = pd.DataFrame(data)
    df.fillna("", inplace=True)  # Preencher valores ausentes com string vazia
    return df

# Exemplo de uso:
df = parse_taxonomy("viral_taxonomy.xml")
df.to_csv("tax.tsv", sep="\t", index=False)
