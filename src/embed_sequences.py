from transformers import EsmModel, EsmTokenizer
import torch
import numpy as np
from Bio import SeqIO

# Initialise protein family
protein_family = "protease"

# Set GPU for more efficient computation
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"🔌 Using device: {device}")

# Load ESM model & tokeniser. Converts amino acid letters into a numerical tokens.
tokeniser = EsmTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")

#load the respective ESM model for the tokeniser
model = EsmModel.from_pretrained("facebook/esm2_t33_650M_UR50D")
model = model.to(device)

# Turn on evaluation mode
model.eval()

def embed_sequence(seq):
    """
    Takes a string of amino acids and converts into a numerical embedding sequence.
    
    Parameters: Sequence (str)
    
    Returns: An embedding: 1D numpy array with a shape equal to hidden size of model (numpy.ndarray)
    """
    tokens = tokeniser(seq, return_tensors="pt", add_special_tokens=True)
    tokens = {k: v.to(device) for k, v in tokens.items()}  # Move to GPU
    with torch.no_grad():
        output = model(**tokens)
    embedding = output.last_hidden_state.mean(dim=1).squeeze().cpu().numpy()
    return embedding

# Embed all sequences
embeddings = []
names = []

fasta_path = f"./data/initial_proteins/{protein_family}/combined_{protein_family}s.fasta"
for record in SeqIO.parse(fasta_path, "fasta"):
    emb = embed_sequence(str(record.seq))
    embeddings.append(emb)
    names.append(record.id)

# Optional: save
np.savez_compressed(f"./data/initial_proteins/{protein_family}/{protein_family}_esm2_embeddings.npz",
                    names=names, embeddings=np.vstack(embeddings))
