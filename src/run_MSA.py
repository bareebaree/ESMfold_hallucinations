"""
This script runs a multiple sequence alignment using a locally installed MAFFT. The input is a series of defined protein families that are desired
to be embdedded and clustered according to similiarity. 

"""

import os
import subprocess
from Bio import SeqIO

# Define protein families as list

protein_families = ["protease", "kinase", "dehydrogenase"]

# For loop to open protein fasta files, clean them, then run mafft on them all.

for protein in protein_families:
    print(f"\n🔍 Processing: {protein}")

    input_fasta = f"./data/initial_proteins/{protein}/combined_{protein}s.fasta"
    cleaned_fasta = f"./data/initial_proteins/{protein}/cleaned_{protein}s.fasta"
    aligned_fasta = f"./data/initial_proteins/{protein}/aligned_{protein}.fasta"

    # Step 1: Remove sequences containing 'U'
    with open(cleaned_fasta, "w") as out_handle:
        kept = 0
        skipped = 0
        for record in SeqIO.parse(input_fasta, "fasta"):
            if "U" not in record.seq:
                SeqIO.write(record, out_handle, "fasta")
                kept += 1
            else:
                skipped += 1
        print(f"🧼 Removed {skipped} sequences with 'U'; kept {kept}")

    # Step 2: Run MAFFT on cleaned sequences
    with open(aligned_fasta, "w") as out_f:
        subprocess.run(["mafft", "--auto", cleaned_fasta], stdout=out_f)
    print(f"✅ MAFFT alignment complete → {aligned_fasta}")