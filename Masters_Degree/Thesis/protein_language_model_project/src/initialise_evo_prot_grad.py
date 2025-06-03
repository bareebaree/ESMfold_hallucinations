import evo_prot_grad
from experts.EMS_expert_modified import EsmExpertModified
import time
import numpy as np
import os
import pandas as pd
import glob
from collections import deque

# === CONFIGURATION ===
"""
Set initial directories to be used. Script will utilise the data_dir as an input, and output to results_dir. If using running multiple loops in parallel, it is prudent to 
set different results directories. Additionally, alter your protein folder name, this is called "dehydrogenase/dehydrogenase_fastas" in this.
initial_fastas_dir = './data/initial_proteins/dehydrogenase/dehydrogenase_fastas/most_distant_sequences'
data_dir = './data/iterations/'
results_dir = './results/'
os.makedirs(data_dir, exist_ok=True)
os.makedirs(results_dir, exist_ok=True)
"""
# === INITIALISE EXPERT ===
# Define and instantiate your expert from the experts in the ./experts folder. For the purposes of this project, only the ESM2 expert was used, 
# as the intention was to optimise and focus on a single model's predictions. To utilise CPU instead, change device to 'cpu'.
ems_expert = EsmExpertModified(scoring_strategy='pseudolikelihood_ratio', temperature=1.0, device = 'cuda')

def get_latest_iteration_log(iteration_number: int) -> str:
    """
    Obtains filepath for logfile for any given iteration of the MCMC loop. 
    This is used to extract the filepath for the previous iteration, so the highest scoring sequence can be extracted from it.
    """
    pattern = os.path.join(results_dir, f"evoprotgrad_run_log_iter{iteration_number}_*.csv")
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(f"No log file found for iteration {iteration_number}")
    latest = max(matches, key=os.path.getmtime)
    return latest

def extract_best_variant_to_fasta(log_csv_path: str, output_fasta_path: str, iteration: int):
    """
    From the extracted sequence from get_latest_iteration_log(), extract the sequence paired with the highest score from this run
    and write to fasta, using the iteration, chain, and score as the header for the fasta.
    """
    
    df = pd.read_csv(log_csv_path)
    if df.empty:
        raise ValueError(f"Log {log_csv_path} is empty.")
    best_row = df.loc[df['score'].idxmax()]
    best_seq = best_row['sequence']
    best_score = best_row['score']
    best_chain = best_row['chain_idx']
    with open(output_fasta_path, 'w') as f:
        f.write(f">iteration_{iteration}_chain{best_chain}_score{best_score:.4f}\n")
        f.write(best_seq + "\n")
    print(f"✅ Wrote best variant from {log_csv_path} to {output_fasta_path}")

# === MAIN LOOP OVER ALL FASTA FILES ===
fasta_files = glob.glob(os.path.join(initial_fastas_dir, '*.fasta'))

for initial_fasta in fasta_files:
    # Loop to call the DirectedEvolution class on each fasta in the target fasta file directory. pdb_id is sliced off at the :4 position and
    # extracted for logging.
    pdb_id = os.path.basename(initial_fasta)[:4]
    print(f"\n🚀 Starting evolution for {pdb_id}...")

    # === INITIAL ITERATION ===
    # Initialisation MCMC run.
    evo_init = evo_prot_grad.DirectedEvolution(
        wt_fasta=initial_fasta,
        output='all',
        experts=[ems_expert],
        parallel_chains=4,
        n_steps=100,
        max_mutations=8,
        verbose=True,
        results_dir=results_dir
    )
    evo_init.pdb_id = pdb_id  # store for logging
    evo_init(iteration_number=0)

    init_log_path = get_latest_iteration_log(iteration_number=0)
    init_fasta = os.path.join(data_dir, f'{pdb_id}_iteration_0.fasta')
    extract_best_variant_to_fasta(init_log_path, init_fasta, iteration=0)

    # initialise logging for iterations
    iteration = 1

    #define number of scores to take average of last n for. This is used to allow the loop to end when the model converges.
    recent_best_scores = deque(maxlen=10)
    wt_fasta = init_fasta

    # === EVOLUTION LOOP ===
    # DirectedEvolution() call for the main MCMC loop.
    while True:
        print(f"\n🧬 Running evo iteration {iteration} for {pdb_id}...")

        evo = evo_prot_grad.DirectedEvolution(
            wt_fasta=wt_fasta,
            output='all',
            experts=[ems_expert],
            parallel_chains=4,
            n_steps=50,
            max_mutations=4,
            verbose=True,
            results_dir=results_dir
        )
        evo.pdb_id = pdb_id  # for logging
        evo(iteration_number=iteration)

        log_path = get_latest_iteration_log(iteration)
        fasta_out = os.path.join(data_dir, f'{pdb_id}_iteration_{iteration}.fasta')
        extract_best_variant_to_fasta(log_path, fasta_out, iteration=iteration)

        df = pd.read_csv(log_path)
        best_score = df['score'].max()
        recent_best_scores.append(best_score)

        # This block of code ensures the model stops at close to convergance (average best PoE score of last ten runs is ≤ 0.01.
        # 0.01 is used to approximate convergance as run will take far too long to reach 0, if it can at all. 
        if len(recent_best_scores) == 10:
            avg_score = sum(recent_best_scores) / 10
            print(f"📊 Average best score over last 10 iterations: {avg_score:.4f}")
            if avg_score <= 0.01:
                print(f"\n🛑 Stopping: Average best score over last 10 runs is ≤ 0.01")
                break

        wt_fasta = fasta_out
        iteration += 1

    print(f"\n✅ Completed evolution for {pdb_id}.")
