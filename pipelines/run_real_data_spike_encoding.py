import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Add project root to path
sys.path.append(os.getcwd())

from src.cstdp.stdp import CausalSTDP
from src.utils.spike_encoding import calculate_adaptive_thresholds, plot_raster

def run_spike_encoding(input_file: str, dataset_name: str):
    """
    Runs spike encoding on a real dataset and generates visuals.
    """
    print(f"--- Processing {dataset_name} ---")
    
    # 1. Load Data
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return
        
    df = pd.read_csv(input_file, index_col=0)
    print(f"Loaded {df.shape[0]} genes, {df.shape[1]} samples.")
    
    data = df.values
    genes = df.index.values
    n_genes, n_timepoints = data.shape
    
    # 2. Time Points (Pseudo-time: just indices 0, 1, 2...)
    # In real analysis, we would parse column names for time metadata.
    # Here we assume columns are ordered time-series (or we treat them as such for demonstration).
    time_points = np.arange(n_timepoints)
    
    # 3. Calculate Thresholds
    # Using sigma=1.5 as in synthetic test
    thresholds = calculate_adaptive_thresholds(data, sigma=1.5)
    
    # 4. Compute Spikes
    cstdp = CausalSTDP() # params dont matter for encoding, only for learning
    spike_trains = cstdp.compute_spike_times(data, time_points, thresholds)
    
    # 5. Visuals
    output_dir = f"visuals/real_data/{dataset_name}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Raster Plot
    plt.figure(figsize=(12, 8))
    plot_raster(spike_trains, (0, min(200, n_timepoints)))
    plt.title(f"Spike Raster: {dataset_name} (First 200 pts)")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/raster_plot.png")
    plt.close()
    
    # Activation Order (First spike time per gene)
    first_spikes = []
    for i, sp in enumerate(spike_trains):
        if len(sp) > 0:
            first_spikes.append(sp[0])
        else:
            first_spikes.append(np.nan)
            
    first_spikes = np.array(first_spikes)
    
    # Sort genes by first activation
    sorted_indices = np.argsort(first_spikes)
    sorted_genes = genes[sorted_indices]
    sorted_times = first_spikes[sorted_indices]
    
    # Filter NaNs
    valid_mask = ~np.isnan(sorted_times)
    sorted_genes = sorted_genes[valid_mask]
    sorted_times = sorted_times[valid_mask]
    
    plt.figure(figsize=(10, 12))
    plt.scatter(sorted_times, range(len(sorted_genes)), alpha=0.7)
    plt.yticks(range(len(sorted_genes)), sorted_genes, fontsize=8)
    plt.xlabel("Time of First Spike")
    plt.title(f"Gene Activation Order: {dataset_name}")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/activation_order.png")
    plt.close()
    
    print(f"Visuals saved to {output_dir}")

if __name__ == "__main__":
    # Primary
    run_spike_encoding("data/processed/GSE215865_subset.csv", "GSE215865")
    
    # Secondary
    run_spike_encoding("data/processed/GSE157859_subset.csv", "GSE157859")
