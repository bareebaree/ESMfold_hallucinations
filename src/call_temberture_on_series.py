

import pandas as pd
from temBERTure import TemBERTure

import sys
sys.path.append("/home/james/TemBERTure")  # or wherever your temBERTure folder lives


# Load the CSV of best sequences per iteration (adjusted to WSL path)
input_csv = "/mnt/c/Users/james/Masters_Degree/Thesis/protein_language_model_project/dehydrogenase_results/dehydrogenase_best_scores_per_iteration.csv"
df = pd.read_csv(input_csv)

# Initialize TemBERTureTM replicas
model_replica1 = TemBERTure(
    adapter_path='./temBERTure_TM/replica1/',
    device='cuda',
    batch_size=16,
    task='regression'
)
model_replica2 = TemBERTure(
    adapter_path='./temBERTure_TM/replica2/',
    device='cuda',
    batch_size=16,
    task='regression'
)
model_replica3 = TemBERTure(
    adapter_path='./temBERTure_TM/replica3/',
    device='cuda',
    batch_size=16,
    task='regression'
)

# To store results per sequence
sequence_results = []

# To compute cumulative scores
cumulative_scores = {}

# Process each sequence
for _, row in df.iterrows():
    pdb_id = row['pdb_id']
    iteration = row['iteration']
    sequence = row['sequence']
    
for _, row in df.iterrows():
    pdb_id = row['pdb_id']
    iteration = row['iteration']
    sequence = row['sequence']
    
    # Predict with each replica (extract scalar from 1-element list)
    pred1 = model_replica1.predict(sequence)[0]
    pred2 = model_replica2.predict(sequence)[0]
    pred3 = model_replica3.predict(sequence)[0]
    
    # Compute average prediction
    avg_prediction = (pred1 + pred2 + pred3) / 3

    
    # Append to results
    sequence_results.append({
        'pdb_id': pdb_id,
        'iteration': iteration,
        'sequence': sequence,
        'avg_prediction': avg_prediction
    })
    
    # Update cumulative score
    if pdb_id not in cumulative_scores:
        cumulative_scores[pdb_id] = 0.0
    cumulative_scores[pdb_id] += avg_prediction

# Save detailed predictions to CSV (using full WSL file path)
sequence_results_csv = "/mnt/c/Users/james/Masters_Degree/Thesis/protein_language_model_project/dehydrogenase_results/temBERTure_sequence_predictions.csv"
sequence_results_df = pd.DataFrame(sequence_results)
sequence_results_df.to_csv(sequence_results_csv, index=False)
print(f"Saved per-sequence predictions to '{sequence_results_csv}'")

# Save cumulative scores to another CSV (using full WSL file path)
cumulative_scores_csv = "/mnt/c/Users/james/Masters_Degree/Thesis/protein_language_model_project/dehydrogenase_results/temBERTure_cumulative_scores.csv"
cumulative_scores_df = pd.DataFrame([
    {"pdb_id": pdb_id, "cumulative_score": cumulative_score}
    for pdb_id, cumulative_score in cumulative_scores.items()
])
cumulative_scores_df.to_csv(cumulative_scores_csv, index=False)
print(f"Saved cumulative scores to '{cumulative_scores_csv}'")

