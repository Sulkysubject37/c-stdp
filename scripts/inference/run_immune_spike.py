import sys
import os
import pandas as pd
import numpy as np
from scipy.stats import zscore

# Add project root to path
sys.path.append(os.getcwd())

from src.cstdp.core import CausalSTDP

def run_immune_spike():
    input_file = "data/processed/GSE215865_immune_subset.csv"
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        sys.exit(1)

    print(f"Loading {input_file}...")
    df = pd.read_csv(input_file, index_col=0)
    data = df.values
    genes = df.index.values
    n_genes, n_timepoints = data.shape
    
    print(f"Data Shape: {data.shape}")
    
    # Handle NaNs in input data
    if np.isnan(data).any():
        print("WARNING: Input data contains NaNs. Filling with 0.")
        data = np.nan_to_num(data, nan=0.0)
    
    # 1. Compute Derivative (first difference)
    print("Computing temporal derivative...")
    derivatives = np.gradient(data, axis=1)
    
    # 2. Gene-wise Normalization (Z-score of derivative)
    print("Applying gene-wise Z-score normalization to derivatives...")
    # zscore computes (x - mean) / std along axis
    norm_derivatives = zscore(derivatives, axis=1)
    
    # Handle NaNs (const genes result in nan zscore)
    norm_derivatives = np.nan_to_num(norm_derivatives, nan=0.0)
    
    # 3. Thresholding
    print(f"Normalized Derivative range: {norm_derivatives.min():.2f} to {norm_derivatives.max():.2f}")
    
    SIGMA = 1.5
    print(f"Applying threshold Z > {SIGMA}...")
    
    spike_mask = norm_derivatives > SIGMA
    spike_counts = np.sum(spike_mask, axis=1)
    
    # Prepare spike trains for C-STDP (list of times)
    time_points = np.arange(n_timepoints)
    spike_trains = []
    for i in range(n_genes):
        times = time_points[spike_mask[i, :]]
        spike_trains.append(times)
        
    # Stats
    avg_spikes = np.mean(spike_counts)
    total_spikes = np.sum(spike_counts)
    genes_with_spikes = np.sum(spike_counts > 0)
    
    print(f"Total Spikes: {total_spikes}")
    print(f"Avg Spikes/Gene: {avg_spikes:.2f}")
    print(f"Genes with Spikes: {genes_with_spikes}/{n_genes}")
    
    if total_spikes == 0:
        print(f"STOP: No bursts detected with Z-score > {SIGMA}.")
        max_z = np.max(norm_derivatives, axis=1)
        print("Top 5 Max Z-scores:")
        print(np.sort(max_z)[-5:])
        sys.exit(1)
    
    # Spike counts per gene (top 10)
    spike_counts_named = sorted(zip(genes, spike_counts), key=lambda x: x[1], reverse=True)
    print("\nTop 10 Spiking Genes:")
    for gene, count in spike_counts_named[:10]:
        print(f"  {gene}: {count}")

    # Onset distribution
    first_spikes = []
    for s in spike_trains:
        if len(s) > 0:
            first_spikes.append(s[0])
            
    print("\nTemporal Onset Distribution (First Spike Time):")
    if len(first_spikes) > 0:
        hist, bins = np.histogram(first_spikes, bins=10)
        for i in range(len(hist)):
            print(f"  T={bins[i]:.1f}-{bins[i+1]:.1f}: {hist[i]} genes")
    else:
        print("  No start times detected.")

if __name__ == "__main__":
    run_immune_spike()
