import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy.stats import zscore

def get_onset_times(input_file):
    df = pd.read_csv(input_file, index_col=0)
    data = np.nan_to_num(df.values, nan=0.0)
    genes = df.index.values
    time_points = np.arange(data.shape[1])
    
    derivatives = np.gradient(data, axis=1)
    norm_derivatives = np.nan_to_num(zscore(derivatives, axis=1), nan=0.0)
    spike_mask = norm_derivatives > 1.5
    spike_counts = np.sum(spike_mask, axis=1)
    
    # Subset
    SUBSET_SIZE = 50 # Larger for comparison
    top_indices = np.argsort(spike_counts)[-SUBSET_SIZE:]
    
    genes_sub = genes[top_indices]
    spike_mask_sub = spike_mask[top_indices]
    
    onset_times = {}
    for i, gene in enumerate(genes_sub):
        spikes = time_points[spike_mask_sub[i, :]]
        if len(spikes) > 0:
            onset_times[gene] = float(spikes[0])
            
    return onset_times

def compare_cascades(file1, name1, file2, name2, output_dir):
    print(f"Comparing Cascades: {name1} vs {name2}...")
    
    onset1 = get_onset_times(file1)
    onset2 = get_onset_times(file2)
    
    # Find shared genes
    shared_genes = set(onset1.keys()).intersection(set(onset2.keys()))
    print(f"Shared genes in top active subset: {len(shared_genes)}")
    
    # Prepare Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 8), sharey=True)
    
    # Plot 1
    sorted_genes1 = sorted(onset1.keys(), key=lambda g: onset1[g])
    times1 = [onset1[g] for g in sorted_genes1]
    ranks1 = range(len(sorted_genes1))
    
    axes[0].scatter(times1, ranks1, c='blue', alpha=0.7)
    axes[0].set_yticks(ranks1)
    axes[0].set_yticklabels(sorted_genes1, fontsize=8)
    axes[0].set_title(f"{name1} Cascade")
    axes[0].set_xlabel("Onset Time")
    axes[0].grid(True, alpha=0.3)
    
    # Highlight shared in Red
    for i, g in enumerate(sorted_genes1):
        if g in shared_genes:
            axes[0].scatter([times1[i]], [ranks1[i]], c='red', s=50, edgecolors='black')

    # Plot 2
    sorted_genes2 = sorted(onset2.keys(), key=lambda g: onset2[g])
    times2 = [onset2[g] for g in sorted_genes2]
    ranks2 = range(len(sorted_genes2))
    
    axes[1].scatter(times2, ranks2, c='green', alpha=0.7)
    axes[1].set_yticks(ranks2)
    axes[1].set_yticklabels(sorted_genes2, fontsize=8)
    axes[1].set_title(f"{name2} Cascade")
    axes[1].set_xlabel("Onset Time")
    axes[1].grid(True, alpha=0.3)
    
    # Highlight shared
    for i, g in enumerate(sorted_genes2):
        if g in shared_genes:
            axes[1].scatter([times2[i]], [ranks2[i]], c='red', s=50, edgecolors='black')
            
    plt.suptitle("Side-by-Side Infection Cascades (Red = Shared Regulators)", fontsize=14)
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, "compare_cascades.png")
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved comparison to {output_path}")

if __name__ == "__main__":
    compare_cascades(
        "data/processed/GSE215865_immune_subset.csv", "GSE215865",
        "data/processed/GSE157859_immune_subset.csv", "GSE157859",
        "visuals/symposium"
    )
