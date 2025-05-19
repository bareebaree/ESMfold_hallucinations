"""
This script is required to initialise EvoProtGrad and have it start running. 

"""
import evo_prot_grad
from experts.EMS_expert_modified import EsmExpertModified
import time
current_timestamp = time.time()
print("Timestamp:", current_timestamp)

# Instantiate the expert
ems_expert = EsmExpertModified(scoring_strategy='pseudolikelihood_ratio', temperature=1.0)

def read_fasta_sequence(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    # Remove header lines and join the sequence lines
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

# Define the FASTA path. Replace this with the filepath to your input fasta, which is best to have in the ./data directory in project.
# In this method, the best sequence of the previous iteration will be used as the input.
wt_fasta = '/mnt/c/Users/james/Masters_Degree/Thesis/protein_language_model_project/data/iterations/Iteration_4_sequence.fasta'

# Initialise the DirectedEvolution object.
# Proposed methodology in this project is to do a 200 step, max mutations=8 initialisation run, then subsequent iterations run at 100 steps and
# 4 max mutations.
evo = evo_prot_grad.DirectedEvolution(
    wt_fasta=wt_fasta,
    output='all',
    experts=[ems_expert],
    parallel_chains=4,
    n_steps=50,
    max_mutations=4,
    verbose=True
)

# run
variants, scores = evo()

# Save results
evo.save_results(
    csv_filename= f'/mnt/c/Users/james/Masters_Degree/Thesis/protein_language_model_project/results/evoprotgrad_output_5_{current_timestamp}.csv',
    variants=variants,
    scores=scores,
    n_seqs_to_keep=50*4
)


