import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import zscore
import os

def plot_onset_timeline(input_file, output_dir):
    print("Generating Absolute Onset Timeline for GSE157859...")
    
    # 1. Load Data & Onset
    df = pd.read_csv(input_file, index_col=0)
    data = np.nan_to_num(df.values, nan=0.0)
    time_points = np.arange(data.shape[1]) # True staged time
    
    derivatives = np.gradient(data, axis=1)
    norm_derivatives = np.nan_to_num(zscore(derivatives, axis=1), nan=0.0)
    spike_mask = norm_derivatives > 1.5
    
    onset_times = []
    genes = []
    
    # Use top 100 spiking genes for clarity
    spike_counts = np.sum(spike_mask, axis=1)
    top_indices = np.argsort(spike_counts)[-100:]
    
    for i in top_indices:
        spikes = time_points[spike_mask[i, :]]
        if len(spikes) > 0:
            onset_times.append(spikes[0])
            genes.append(df.index[i])
            
    # Sort
    sorted_pairs = sorted(zip(genes, onset_times), key=lambda x: x[1])
    sorted_genes = [x[0] for x in sorted_pairs]
    sorted_times = [x[1] for x in sorted_pairs]
    
    # 2. Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    ax.scatter(sorted_times, range(len(sorted_genes)), c='purple', alpha=0.6, s=20)
    
    # Connect to Y-axis
    for i, t in enumerate(sorted_times):
        ax.hlines(i, 0, t, color='lightgray', linewidth=0.5)
        
    ax.set_yticks(range(len(sorted_genes)))
    ax.set_yticklabels(sorted_genes, fontsize=6)
    
    ax.set_xlabel("Time (Infection Stages)")
    ax.set_ylabel("Immune Genes (Ordered by Onset)")
    ax.set_title("Absolute Temporal Onset (GSE157859)")
    
    ax.grid(True, axis='x', linestyle='--', alpha=0.3)
    
    output_path = os.path.join(output_dir, "GSE157859_onset_timeline.png")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved timeline to {output_path}")

if __name__ == "__main__":
    plot_onset_timeline(
        "data/processed/GSE157859_immune_subset.csv",
        "visuals/symposium/GSE157859"
    )
