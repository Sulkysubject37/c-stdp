import sys
import os
import pandas as pd
import numpy as np
from scipy.stats import zscore

# Add project root to path
sys.path.append(os.getcwd())

from src.visualization.temporal_cascade import plot_temporal_cascade

def generate_cascade(dataset_name, input_file, adj_file, output_file):
    print(f"\n--- Generating Cascade Plot for {dataset_name} ---")
    
    if not os.path.exists(input_file):
        print(f"Input file missing: {input_file}")
        return
    if not os.path.exists(adj_file):
        print(f"Adjacency file missing: {adj_file}")
        return

    # 1. Re-calculate Spike Times
    df = pd.read_csv(input_file, index_col=0)
    data = np.nan_to_num(df.values, nan=0.0)
    genes = df.index.values
    time_points = np.arange(data.shape[1])
    
    derivatives = np.gradient(data, axis=1)
    norm_derivatives = np.nan_to_num(zscore(derivatives, axis=1), nan=0.0)
    
    SIGMA = 1.5
    spike_mask = norm_derivatives > SIGMA
    spike_counts = np.sum(spike_mask, axis=1)
    
    # Identify Top 200 Subset (must match inference)
    SUBSET_SIZE = 200
    top_indices = np.argsort(spike_counts)[-SUBSET_SIZE:]
    
    genes_sub = genes[top_indices]
    spike_mask_sub = spike_mask[top_indices]
    
    # Extract Onset Times (First Spike)
    onset_times = {}
    for i, gene in enumerate(genes_sub):
        spikes = time_points[spike_mask_sub[i, :]]
        if len(spikes) > 0:
            onset_times[gene] = float(spikes[0])
        else:
            # Should not happen for top spiking genes, but technically possible if count > 0
            onset_times[gene] = float(time_points[-1]) # Fallback
            
    # 2. Load Causal Links from Adjacency
    adj_df = pd.read_csv(adj_file, index_col=0)
    
    causal_links = []
    # Only iterate over the subset genes present in adj_df
    # (adj_df should be 200x200)
    for u in adj_df.index:
        for v in adj_df.columns:
            w = adj_df.loc[u, v]
            if w > 0:
                causal_links.append((u, v, w))
                
    # 3. Plot
    plot_temporal_cascade(onset_times, causal_links, output_file, top_n_links=30)

if __name__ == "__main__":
    # Primary
    generate_cascade("GSE215865", 
                     "data/processed/GSE215865_immune_subset.csv",
                     "visuals/real_data/GSE215865_Immune/inferred_grn_adj.csv",
                     "visuals/real_data/GSE215865_Immune/temporal_cascade.png")
    
    # Supplementary
    generate_cascade("GSE157859", 
                     "data/processed/GSE157859_immune_subset.csv",
                     "visuals/real_data/GSE157859_Immune/inferred_grn_adj.csv",
                     "visuals/real_data/GSE157859_Immune/temporal_cascade.png")
