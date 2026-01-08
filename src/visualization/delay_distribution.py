import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy.stats import zscore

def generate_delay_distribution(dataset_name, input_file, adj_file, output_dir):
    print(f"Generating Delay Distribution for {dataset_name}...")
    
    # 1. Load Data
    df = pd.read_csv(input_file, index_col=0)
    data = np.nan_to_num(df.values, nan=0.0)
    genes = df.index.values
    time_points = np.arange(data.shape[1])
    
    derivatives = np.gradient(data, axis=1)
    norm_derivatives = np.nan_to_num(zscore(derivatives, axis=1), nan=0.0)
    spike_mask = norm_derivatives > 1.5
    spike_counts = np.sum(spike_mask, axis=1)
    
    SUBSET_SIZE = 200
    top_indices = np.argsort(spike_counts)[-SUBSET_SIZE:]
    genes_sub = genes[top_indices]
    spike_mask_sub = spike_mask[top_indices]
    
    onset_times = {}
    for i, gene in enumerate(genes_sub):
        spikes = time_points[spike_mask_sub[i, :]]
        if len(spikes) > 0:
            onset_times[gene] = float(spikes[0])
            
    # 2. Load Adjacency
    adj_df = pd.read_csv(adj_file, index_col=0)
    
    delays = []
    
    for u in adj_df.index:
        for v in adj_df.columns:
            w = adj_df.loc[u, v]
            # Only consider links with some weight and where both nodes have onset
            if w > 0.01 and u in onset_times and v in onset_times:
                dt = onset_times[v] - onset_times[u]
                # Filter for causal (positive) delays roughly matching STDP window
                # STDP window is ~10-20. 
                # We plot all to show distribution.
                delays.append(dt)
                
    # 3. Plot
    plt.figure(figsize=(8, 6))
    plt.hist(delays, bins=50, color='teal', alpha=0.7, edgecolor='black')
    
    median_delay = np.median(delays)
    plt.axvline(median_delay, color='red', linestyle='--', label=f'Median: {median_delay:.1f}')
    
    plt.xlabel("Temporal Delay (Pseudo-time units)")
    plt.ylabel("Frequency of Interactions")
    plt.title(f"Inferred Causal Delays: {dataset_name}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    output_path = os.path.join(output_dir, "delay_distribution.png")
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved delay plot to {output_path}")

if __name__ == "__main__":
    generate_delay_distribution(
        "GSE215865",
        "data/processed/GSE215865_immune_subset.csv",
        "visuals/real_data/GSE215865_Immune/inferred_grn_adj.csv",
        "visuals/symposium/GSE215865"
    )
    generate_delay_distribution(
        "GSE157859",
        "data/processed/GSE157859_immune_subset.csv",
        "visuals/real_data/GSE157859_Immune/inferred_grn_adj.csv",
        "visuals/symposium/GSE157859"
    )
