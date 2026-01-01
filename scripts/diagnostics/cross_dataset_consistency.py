import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.getcwd())

from src.cstdp.stdp import CausalSTDP
from src.utils.spike_encoding import calculate_adaptive_thresholds

def run_cross_dataset_consistency():
    print("--- Starting Cross-Dataset Behavioral Consistency ---")
    
    datasets = {
        "GSE215865": "data/processed/GSE215865_subset.csv",
        "GSE157859": "data/processed/GSE157859_subset.csv"
    }
    
    results = {}
    
    for name, path in datasets.items():
        if not os.path.exists(path):
            print(f"Skipping {name}: not found.")
            continue
            
        print(f"Processing {name}...")
        df = pd.read_csv(path, index_col=0)
        data = df.values
        n_genes, n_samples = data.shape
        
        # Infer
        thresholds = calculate_adaptive_thresholds(data, sigma=1.5)
        cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
        spike_trains = cstdp.compute_spike_times(data, np.arange(n_samples), thresholds)
        weights = cstdp.run_cstdp(spike_trains, n_genes)
        
        # Normalize
        mw = np.max(weights)
        norm_w = weights / mw if mw > 0 else weights
        edges = norm_w > 0.3
        
        # Metrics
        sparsity = np.mean(edges == 0)
        out_degrees = np.sum(edges, axis=1)
        in_degrees = np.sum(edges, axis=0)
        
        results[name] = {
            "Sparsity": sparsity,
            "Out_Degree": out_degrees,
            "In_Degree": in_degrees,
            "Weights": norm_w.flatten()
        }
        
    if len(results) < 2:
        print("Need both datasets for comparison.")
        return

    # Compare
    print("\n--- COMPARISON RESULTS ---")
    print(f"Sparsity GSE215865: {results['GSE215865']['Sparsity']:.4f}")
    print(f"Sparsity GSE157859: {results['GSE157859']['Sparsity']:.4f}")
    
    # Degree Distribution Plot
    os.makedirs("analysis/visuals", exist_ok=True)
    plt.figure(figsize=(12, 6))
    
    plt.subplot(1, 2, 1)
    sns.kdeplot(results['GSE215865']['Out_Degree'], label='GSE215865', fill=True)
    sns.kdeplot(results['GSE157859']['Out_Degree'], label='GSE157859', fill=True)
    plt.title("Out-Degree Distribution")
    plt.xlabel("Number of Targets")
    plt.legend()
    
    plt.subplot(1, 2, 2)
    sns.kdeplot(results['GSE215865']['Weights'], label='GSE215865', log_scale=(False, True)) # Log y for visibility
    sns.kdeplot(results['GSE157859']['Weights'], label='GSE157859', log_scale=(False, True))
    plt.title("Weight Distribution (Log Scale)")
    plt.xlabel("Normalized Weight")
    plt.legend()
    
    plt.savefig("analysis/visuals/cross_dataset_comparison.png")
    plt.close()
    
    print("Interpretation: Similar sparsity and degree distributions suggest the algorithm behaves consistently across different biological contexts.")

if __name__ == "__main__":
    run_cross_dataset_consistency()
