import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
from scipy.stats import zscore
import os

# Functional Keywords for Coloring
FUNCTION_MAP = {
    'IFN': ['IFN', 'ISG', 'IRF', 'STAT', 'OAS', 'MX', 'IFI'],
    'Cytokine': ['IL', 'CXC', 'CCL', 'TNF', 'TGF', 'CSF'],
    'Sensing': ['TLR', 'NLR', 'RIG', 'DDX', 'MAVS'],
    'Stress': ['HSP', 'HMOX', 'EIF2', 'ATF']
}

def get_functional_class(gene):
    gene_upper = str(gene).upper()
    for cat, keywords in FUNCTION_MAP.items():
        for k in keywords:
            if k in gene_upper:
                return cat
    return 'Other'

COLOR_MAP = {
    'IFN': 'blue',
    'Cytokine': 'red',
    'Sensing': 'green',
    'Stress': 'orange',
    'Other': 'gray'
}

def generate_activation_timeline(input_file, dataset_name, output_dir):
    print(f"Generating Activation Timeline for {dataset_name}...")
    
    # 1. Load Data & Detect Onset (Same logic as inference)
    df = pd.read_csv(input_file, index_col=0)
    data = np.nan_to_num(df.values, nan=0.0)
    genes = df.index.values
    time_points = np.arange(data.shape[1])
    
    # Spike Encoding
    derivatives = np.gradient(data, axis=1)
    norm_derivatives = np.nan_to_num(zscore(derivatives, axis=1), nan=0.0)
    
    spike_mask = norm_derivatives > 1.5
    spike_counts = np.sum(spike_mask, axis=1)
    
    # Subset to Top 200 (Consistency with GRN)
    SUBSET_SIZE = 200
    top_indices = np.argsort(spike_counts)[-SUBSET_SIZE:]
    
    genes_sub = genes[top_indices]
    spike_mask_sub = spike_mask[top_indices]
    
    # Get Onset Times
    onset_times = []
    for i, gene in enumerate(genes_sub):
        spikes = time_points[spike_mask_sub[i, :]]
        if len(spikes) > 0:
            onset_times.append((gene, float(spikes[0])))
        else:
            onset_times.append((gene, float(time_points[-1]))) # Fallback
            
    # Sort by Onset
    onset_times.sort(key=lambda x: x[1])
    sorted_genes = [x[0] for x in onset_times]
    sorted_times = [x[1] for x in onset_times]
    
    # Assign Colors
    colors = [COLOR_MAP[get_functional_class(g)] for g in sorted_genes]
    
    # 2. Animation Setup
    fig, ax = plt.subplots(figsize=(10, 8))
    
    ax.set_xlim(0, max(sorted_times) + 5)
    ax.set_ylim(-5, len(sorted_genes) + 5)
    ax.set_xlabel("Pseudo-Time (Onset)")
    ax.set_ylabel("Immune Genes (Sorted by Activation)")
    ax.set_title(f"Immune Activation Timeline: {dataset_name}")
    
    scat = ax.scatter([], [], c=[], s=50, alpha=0.8)
    
    # Add Legend
    markers = [plt.Line2D([0,0],[0,0], color=color, marker='o', linestyle='') for color in COLOR_MAP.values()]
    ax.legend(markers, COLOR_MAP.keys(), loc='upper left', fontsize=8)
    
    def update(frame):
        # Frame represents current time
        # Show all genes activated up to 'frame'
        
        current_time = frame
        
        # Filter points
        x_data = []
        y_data = []
        c_data = []
        
        for i, (g, t) in enumerate(onset_times):
            if t <= current_time:
                x_data.append(t)
                y_data.append(i) # Rank on Y-axis
                c_data.append(colors[i])
                
        if x_data:
            scat.set_offsets(np.c_[x_data, y_data])
            scat.set_color(c_data)
        
        ax.set_title(f"Immune Activation: {dataset_name} (T={current_time:.1f})")
        return scat,

    # Frames: from 0 to max time
    max_t = int(max(sorted_times)) + 2
    ani = animation.FuncAnimation(fig, update, frames=np.linspace(0, max_t, 100), interval=50, blit=False)
    
    # Save Dynamic
    mp4_path = os.path.join(output_dir, "activation_timeline.mp4")
    ani.save(mp4_path, writer='ffmpeg', fps=20)
    print(f"Saved animation to {mp4_path}")
    
    # Save Static (Final Frame)
    update(max_t)
    png_path = os.path.join(output_dir, "activation_timeline_static.png")
    plt.savefig(png_path, dpi=150)
    plt.close()
    print(f"Saved static snapshot to {png_path}")

if __name__ == "__main__":
    generate_activation_timeline(
        "data/processed/GSE215865_immune_subset.csv",
        "GSE215865",
        "visuals/symposium/GSE215865"
    )
    generate_activation_timeline(
        "data/processed/GSE157859_immune_subset.csv",
        "GSE157859",
        "visuals/symposium/GSE157859"
    )
