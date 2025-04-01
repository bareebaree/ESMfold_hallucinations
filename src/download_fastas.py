import requests
import os
import pandas as pd

protein_family = "protease"

# Import structures.tsv and extract accession numbers into a list
df = pd.read_csv(f"./data/initial_proteins/{protein_family}/{protein_family}_structures.tsv", sep="\t")


# Extract the PDB IDs from the 'Accession' column
pdb_ids = df['Accession'].dropna().astype(str).str.lower().unique()

# Output folder
os.makedirs(f"{protein_family}_fastas", exist_ok=True)

for pdb_id in pdb_ids:
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id.upper()}"
    response = requests.get(url)
    if response.status_code == 200:
        with open(f"{protein_family}_fastas/{pdb_id}.fasta", "w") as f:
            f.write(response.text)
        print(f"Downloaded {pdb_id}.fasta")
    else:
        print(f"Failed to download {pdb_id}")
