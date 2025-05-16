#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 23:03:17 2025

@author: james
"""

from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
import matplotlib
from sklearn.metrics import pairwise_distances
import itertools
import numpy as np
import os
import shutil

# Allow matplotlib to use command line
matplotlib.use("Agg")  

# set protein family
protein_family = "kinase"

# Load data and embeddings
data = np.load(f"./data/initial_proteins/{protein_family}/{protein_family}_esm2_embeddings.npz", allow_pickle=True)
names = data["names"]
embeddings = data["embeddings"]

print(f"Loaded {len(embeddings)} embeddings")

# Initialise kmeans clustering

k = 10  # or however many clusters you used
kmeans = KMeans(n_clusters=k, random_state=42)
labels = kmeans.fit_predict(embeddings)
centroids = kmeans.cluster_centers_

# Map labels to names
for name, label in zip(names, labels):
    print(f"{name} → Cluster {label}")
    

# Save clusters as CSV for each PDB ID
df = pd.DataFrame({"name": names, "cluster": labels})
df.to_csv(f"./data/initial_proteins/{protein_family}/{protein_family}_cluster_assignments.csv", index=False)


"""
The below section serves to determine the maximally distant clusters within this protein family.
It computes a distance matrix so their euclidean distance can be computed, and have them sorted 
by distance.
"""

# Compute all pairwise distances between centroids
dist_matrix = pairwise_distances(centroids)

# Get upper triangle (without diagonal) as flattened index pairs
n_clusters = centroids.shape[0]
pairs = list(itertools.combinations(range(n_clusters), 2))

# Extract distances and sort
distances = [(i, j, dist_matrix[i, j]) for i, j in pairs]
distances.sort(key=lambda x: x[2], reverse=True)  # Sort by distance descending

# Top 10 maximally distant cluster pairs
top_10 = distances[:10]
print("Top 10 maximally distant cluster pairs:")
for i, j, d in top_10:
    print(f"Cluster {i} ↔ Cluster {j} → Distance: {d:.4f}")
    
    
# Plot tSNE plot
plt.figure(figsize=(12, 8))

# Reduce both embeddings and centroids to 2D. Uses tSNE to visualise clusters.
tsne = TSNE(n_components=2, random_state=42)
reduced_all = tsne.fit_transform(np.vstack([embeddings, centroids]))
reduced_embeddings = reduced_all[:-len(centroids)]
reduced_centroids = reduced_all[-len(centroids):]

# Plot individual points
sns.scatterplot(x=reduced_embeddings[:, 0], y=reduced_embeddings[:, 1],
                hue=labels, palette="tab10", s=40, legend='full', alpha=0.8)

# Plot centroids
plt.scatter(reduced_centroids[:, 0], reduced_centroids[:, 1],
            c='black', s=200, marker='X', label='Centroids')

# Draw lines between top-10 most distant centroids
for i, j, d in top_10:
    x_coords = [reduced_centroids[i, 0], reduced_centroids[j, 0]]
    y_coords = [reduced_centroids[i, 1], reduced_centroids[j, 1]]
    plt.plot(x_coords, y_coords, color='black', linestyle='--', alpha=0.6)
    mid_x = (x_coords[0] + x_coords[1]) / 2
    mid_y = (y_coords[0] + y_coords[1]) / 2
    plt.text(mid_x, mid_y, f"{d:.1f}", fontsize=8, color="black")

plt.title("t-SNE of ESM2 Embeddings\nWith Top 10 Maximally Distant Cluster Connections")
plt.xlabel("t-SNE 1")
plt.ylabel("t-SNE 2")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(f"./data/initial_proteins/{protein_family}/tsne_top10_distant_clusters.png", dpi=300)
print("✅ Plot saved as 'tsne_top10_distant_clusters.png'")

# Find the closest sequence to each centroid

closest_sequences = []

for cluster_id in range(k):
    # Get indices of points in this cluster
    indices = np.where(labels == cluster_id)[0]
    
    # Get embeddings of those points
    cluster_embeddings = embeddings[indices]
    
    # Compute distances to the centroid
    dists = np.linalg.norm(cluster_embeddings - centroids[cluster_id], axis=1)
    
    # Find the index of the closest embedding
    closest_idx_in_cluster = indices[np.argmin(dists)]
    
    # Clean PDB_ID by taking only the first 4 characters
    raw_name = names[closest_idx_in_cluster]
    cleaned_name = raw_name[:4]  # First 4 characters only
    
    # Store cleaned name and cluster
    closest_sequences.append({
        "PDB_ID": cleaned_name,
        "Cluster": cluster_id
    })

# Create and save DataFrame
closest_df = pd.DataFrame(closest_sequences)
closest_df.to_csv(f"./data/initial_proteins/{protein_family}/closest_to_centroid.csv", index=False)
print("✅ Saved closest sequences to centroids as 'closest_to_centroid.csv'")

# Copy matching fasta files
fasta_src_dir = f"./data/initial_proteins/{protein_family}/{protein_family}_fastas"
fasta_dst_dir = os.path.join(fasta_src_dir, "most_distant_sequences")
os.makedirs(fasta_dst_dir, exist_ok=True)

for entry in closest_sequences:
    pdb_id_lower = entry["PDB_ID"].lower()
    src_fasta = os.path.join(fasta_src_dir, f"{pdb_id_lower}.fasta")
    dst_fasta = os.path.join(fasta_dst_dir, f"{pdb_id_lower}.fasta")
    
    if os.path.exists(src_fasta):
        shutil.copy(src_fasta, dst_fasta)
        print(f"📁 Copied {pdb_id_lower}.fasta to most_distant_sequences/")
    else:
        print(f"⚠️ FASTA file not found for {pdb_id_lower}")

print("✅ FASTA copying complete.")