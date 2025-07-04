from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import pandas as pd

def analyze_sequence(seq):
    analysis = ProteinAnalysis(str(seq))
    aa_count = analysis.count_amino_acids()
    total_len = len(seq)

    num_pos_residues = aa_count.get('K', 0) + aa_count.get('R', 0) + aa_count.get('H', 0)
    num_neg_residues = aa_count.get('D', 0) + aa_count.get('E', 0)
    sec_struc = analysis.secondary_structure_fraction()

    x_ala = aa_count.get('A', 0) / total_len
    x_val = aa_count.get('V', 0) / total_len
    x_ile = aa_count.get('I', 0) / total_len
    x_leu = aa_count.get('L', 0) / total_len
    aliphatic_index = 100 * (x_ala + 2.9 * x_val + 3.9 * (x_ile + x_leu))

    return {
        'aromaticity': analysis.aromaticity(),
        'instability_index': analysis.instability_index(),
        'isoelectric_point': analysis.isoelectric_point(),
        'sec_struc_helix': sec_struc[0],
        'sec_struc_turn': sec_struc[1],
        'sec_struc_sheet': sec_struc[2],
        'flexibility_mean': sum(analysis.flexibility()) / len(analysis.flexibility()),
        'gravy': analysis.gravy(),
        'num_pos_residues': num_pos_residues,
        'num_neg_residues': num_neg_residues,
        'aliphatic_index': aliphatic_index
    }

def analyze_sequences_file(file_path, filetype):
    results = []

    if filetype.lower() == "fasta":
        for record in SeqIO.parse(file_path, "fasta"):
            pdb_id = record.id[:4]
            seq_metrics = analyze_sequence(record.seq)
            seq_metrics['pdb_id'] = pdb_id
            results.append(seq_metrics)

    elif filetype.lower() == "csv":
        df = pd.read_csv(file_path)
        for _, row in df.iterrows():
            pdb_id = row['pdb_id']
            sequence = row['sequence']
            seq_metrics = analyze_sequence(sequence)
            seq_metrics['pdb_id'] = pdb_id
            results.append(seq_metrics)
    else:
        raise ValueError("filetype must be 'fasta' or 'csv'")

    if results:
        df_out = pd.DataFrame(results)
        df_out.to_csv("output.csv", index=False)
        print("Analysis complete. Results saved to 'output.csv'.")
        return df_out
    else:
        print("No results found.")
        return pd.DataFrame()

df_results = analyze_sequences_file("C:\\Users\\james\\Masters_Degree\\Thesis\\protein_language_model_project\\results\\var_wt_sequences.csv", "csv")