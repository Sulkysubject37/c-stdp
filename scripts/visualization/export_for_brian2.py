import sys
import os
import pandas as pd
import numpy as np

sys.path.append(os.getcwd())

from src.cstdp.utils.spike_encoding import calculate_adaptive_thresholds
from src.cstdp.core import CausalSTDP

def export_for_brian2():
    print("--- Exporting Data for Brian2 Simulation ---")
    
    # 1. Load Data and Subgraph
    data_path = "data/processed/GSE215865_subset.csv"
    subgraph_path = "results/brian2/mini_adj.csv"
    
    if not os.path.exists(data_path) or not os.path.exists(subgraph_path):
        print("Required files not found.")
        return
        
    df = pd.read_csv(data_path, index_col=0)
    adj_df = pd.read_csv(subgraph_path, index_col=0)
    selected_genes = adj_df.index.values
    
    # 2. Extract Subgraph Data
    sub_data = df.loc[selected_genes].values
    n_genes, n_samples = sub_data.shape
    
    # 3. Compute Spike Times (Reproducing Inference)
    print("Computing spike times for selected genes...")
    thresholds = calculate_adaptive_thresholds(sub_data, sigma=1.5)
    cstdp = CausalSTDP()
    spike_trains = cstdp.compute_spike_times(sub_data, np.arange(n_samples), thresholds)
    
    # 4. Save
    # We save spike indices (times) and the weight matrix
    # Weights are from the inferred adjacency (fixed)
    weights = adj_df.values
    
    # Convert list of arrays to object array for saving
    spikes_dict = {gene: spikes for gene, spikes in zip(selected_genes, spike_trains)}
    
    output_file = "results/brian2/simulation_data.npz"
    np.savez(output_file, 
             weights=weights, 
             genes=selected_genes, 
             **spikes_dict)
             
    print(f"Exported simulation data to {output_file}")
    print(f"Contains {len(selected_genes)} genes and fixed weights.")

if __name__ == "__main__":
    export_for_brian2()
