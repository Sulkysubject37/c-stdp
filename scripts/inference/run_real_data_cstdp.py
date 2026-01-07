import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy.stats import zscore

# Add project root to path
sys.path.append(os.getcwd())

def vectorized_stdp(spike_trains, A_pos, A_neg, tau_pos, tau_neg, w_max):
    n_genes = len(spike_trains)
    weights = np.zeros((n_genes, n_genes))
    
    print(f"Vectorizing STDP for {n_genes} genes...")
    for i in range(n_genes):
        if i % 20 == 0:
            print(f"  Progress: {i}/{n_genes} genes")
        for j in range(n_genes):
            if i == j:
                continue
            
            sp_i = spike_trains[i]
            sp_j = spike_trains[j]
            
            if len(sp_i) == 0 or len(sp_j) == 0:
                continue
            
            # Matrix of all pairs: dt = t_post (j) - t_pre (i)
            dt = sp_j[np.newaxis, :] - sp_i[:, np.newaxis]
            
            # Causal (dt > 0)
            pos_mask = dt > 0
            # Anti-causal (dt < 0)
            neg_mask = dt < 0
            
            w_sum = 0.0
            if np.any(pos_mask):
                w_sum += np.sum(A_pos * np.exp(-dt[pos_mask] / tau_pos))
            if np.any(neg_mask):
                w_sum -= np.sum(A_neg * np.exp(dt[neg_mask] / tau_neg))
                
            weights[i, j] = w_sum
            
    # Normalize/Clip
    return np.clip(weights, 0, w_max)

def run_cstdp_inference(input_file: str, dataset_name: str):
    """
    Runs C-STDP on real data using Vectorized runner.
    """
    print(f"\n--- Inferring GRN for {dataset_name} (Vectorized) ---")
    
    # 1. Load
    df = pd.read_csv(input_file, index_col=0)
    data = np.nan_to_num(df.values, nan=0.0)
    genes = df.index.values
    n_genes = data.shape[0]
    
    # 2. Spikes
    derivatives = np.gradient(data, axis=1)
    norm_derivatives = np.nan_to_num(zscore(derivatives, axis=1), nan=0.0)
    
    SIGMA = 1.5
    spike_mask = norm_derivatives > SIGMA
    spike_counts = np.sum(spike_mask, axis=1)
    
    # Subset to 200
    SUBSET_SIZE = 200
    top_indices = np.argsort(spike_counts)[-SUBSET_SIZE:]
    genes_sub = genes[top_indices]
    spike_mask_sub = spike_mask[top_indices]
    
    time_points = np.arange(data.shape[1])
    spike_trains = [time_points[spike_mask_sub[i, :]] for i in range(len(genes_sub))]
    
    # 3. Learning (Vectorized)
    # A_pos, A_neg = 0.01, tau = 10 as per previous
    weights = vectorized_stdp(spike_trains, 0.01, 0.01, 10.0, 10.0, 1.0)
    
    # 4. Analysis
    output_dir = f"visuals/real_data/{dataset_name}"
    os.makedirs(output_dir, exist_ok=True)
    
    adj_df = pd.DataFrame(weights, index=genes_sub, columns=genes_sub)
    adj_df.to_csv(f"{output_dir}/inferred_grn_adj.csv")
    
    sparsity = (weights == 0).sum() / weights.size
    print(f"Adjacency Sparsity: {sparsity:.4f}")
    
    out_degree = adj_df.sum(axis=1)
    print("\nTop 10 Regulators (Out-Degree):")
    print(out_degree.sort_values(ascending=False).head(10))
    
    plt.figure()
    plt.hist(out_degree, bins=30)
    plt.title("Out-Degree Distribution (Vectorized)")
    plt.savefig(f"{output_dir}/out_degree_dist.png")
    plt.close()
    print(f"Saved to {output_dir}")

if __name__ == "__main__":
    run_cstdp_inference("data/processed/GSE215865_immune_subset.csv", "GSE215865_Immune")
