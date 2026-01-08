import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_activation_layers(rank_file, output_dir):
    print("Generating Activation Layers Plot for GSE215865...")
    
    if not os.path.exists(rank_file):
        print("Rank file not found.")
        return
        
    df = pd.read_csv(rank_file)
    df = df.sort_values('Percentile')
    
    # Categorize
    # Early: 0-20%
    # Inter: 20-60%
    # Late: 60-100%
    
    early = df[df['Percentile'] <= 20]
    inter = df[(df['Percentile'] > 20) & (df['Percentile'] <= 60)]
    late = df[df['Percentile'] > 60]
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # We plot Gene Index (sorted) vs Percentile
    # Or just blocks?
    # "X-axis: activation percentile, Y-axis: genes"
    # Genes are categorical, too many (3488).
    # We plot them as lines or points sorted by percentile.
    
    y_pos = range(len(df))
    
    ax.scatter(early['Percentile'], range(len(early)), c='#1f77b4', s=5, label='Early (Top 20%)', edgecolors='none')
    
    offset = len(early)
    ax.scatter(inter['Percentile'], range(offset, offset+len(inter)), c='#ff7f0e', s=5, label='Intermediate (20-60%)', edgecolors='none')
    
    offset += len(inter)
    ax.scatter(late['Percentile'], range(offset, offset+len(late)), c='#2ca02c', s=5, label='Late (60-100%)', edgecolors='none')
    
    # Labels
    ax.set_xlabel("Relative Activation Order (Percentile)")
    ax.set_ylabel("Immune Genes (Sorted by Rank)")
    ax.set_title("Rank-Based Activation Layers (Pseudo-time Cohort)")
    
    # Boundaries
    ax.axvline(20, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(60, color='gray', linestyle='--', alpha=0.5)
    
    ax.legend(loc='lower right')
    ax.grid(True, linestyle=':', alpha=0.3)
    
    # Add explicit text warning about time units
    plt.figtext(0.5, 0.01, "Note: Axis represents relative rank order, not absolute time.", 
                ha="center", fontsize=9, style='italic', color='gray')

    output_path = os.path.join(output_dir, "GSE215865_activation_layers.png")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved figure to {output_path}")

if __name__ == "__main__":
    plot_activation_layers(
        "visuals/symposium/GSE215865/activation_ranks.csv",
        "visuals/symposium/GSE215865"
    )
