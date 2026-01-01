import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

# Add project root to path
sys.path.append(os.getcwd())

from src.cstdp.core import CausalSTDP
from src.cstdp.utils.spike_encoding import calculate_adaptive_thresholds

def run_cstdp_inference(input_file: str, dataset_name: str):
    """
    Runs C-STDP on real data to infer GRN.
    """
    print(f"\n--- Inferring GRN for {dataset_name} ---")
    
    # 1. Load
    if not os.path.exists(input_file):
        print(f"Skipping {dataset_name}: File not found.")
        return
        
    df = pd.read_csv(input_file, index_col=0)
    data = df.values
    genes = df.index.values
    n_genes, n_timepoints = data.shape
    
    # 2. Spikes
    thresholds = calculate_adaptive_thresholds(data, sigma=1.5)
    cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
    
    spike_trains = cstdp.compute_spike_times(data, np.arange(n_timepoints), thresholds)
    
    # 3. Learning
    print("Running STDP...")
    weights = cstdp.run_cstdp(spike_trains, n_genes)
    
    # 4. Analysis
    # Save Adjacency
    output_dir = f"visuals/real_data/{dataset_name}"
    os.makedirs(output_dir, exist_ok=True)
    
    adj_df = pd.DataFrame(weights, index=genes, columns=genes)
    adj_df.to_csv(f"{output_dir}/inferred_grn_adj.csv")
    
    # Top Regulators (Out-Degree sum)
    out_degree = adj_df.sum(axis=1).sort_values(ascending=False)
    print("Top 5 Regulators:")
    print(out_degree.head(5))
    
    # Visualization: Heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(adj_df, cmap="viridis")
    plt.title(f"Inferred GRN Weights: {dataset_name}")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/grn_heatmap.png")
    plt.close()
    
    # Visualization: Network Graph (Thresholded)
    # Threshold at 95th percentile of non-zero weights
    flat_w = weights.flatten()
    flat_w = flat_w[flat_w > 0]
    if len(flat_w) > 0:
        thresh = np.percentile(flat_w, 90)
    else:
        thresh = 0.0
        
    G = nx.DiGraph()
    for i in range(n_genes):
        for j in range(n_genes):
            if weights[i, j] > thresh:
                G.add_edge(genes[i], genes[j], weight=weights[i, j])
    
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, k=0.5)
    nx.draw_networkx_nodes(G, pos, node_size=500, node_color='lightblue')
    nx.draw_networkx_edges(G, pos, alpha=0.5, arrows=True)
    nx.draw_networkx_labels(G, pos, font_size=8)
    plt.title(f"Top 10% Strongest Interactions ({dataset_name})")
    plt.savefig(f"{output_dir}/grn_graph.png")
    plt.close()

if __name__ == "__main__":
    run_cstdp_inference("data/processed/GSE215865_subset.csv", "GSE215865")
    run_cstdp_inference("data/processed/GSE157859_subset.csv", "GSE157859")
