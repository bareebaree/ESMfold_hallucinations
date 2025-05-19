import evo_prot_grad
from experts.EMS_expert_modified import EsmExpertModified

# Instantiate the expert
ems_expert = EsmExpertModified(scoring_strategy='pseudolikelihood_ratio', temperature=1.0)

def read_fasta_sequence(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    # Remove header lines and join the sequence lines
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

# ✅ Define the FASTA path
wt_fasta = '/mnt/c/Users/james/Masters_Degree/Thesis/protein_language_model_project/data/initial_proteins/dehydrogenase/dehydrogenase_fastas/most_distant_sequences/1axe.fasta'

# ✅ Initialize the DirectedEvolution object
evo = evo_prot_grad.DirectedEvolution(
    wt_fasta=wt_fasta,
    output='all',
    experts=[ems_expert],
    parallel_chains=4,
    n_steps=300,
    max_mutations=len(read_fasta_sequence(wt_fasta)),
    verbose=True
)

# ✅ Run it
variants, scores = evo()

# ✅ Save results
evo.save_results(
    csv_filename= '/mnt/c/Users/james/Masters_Degree/Thesis/protein_language_model_project/results/evoprotgrad_output.csv',
    variants=variants,
    scores=scores,
    n_seqs_to_keep=100
)


