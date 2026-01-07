import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Tuple

def plot_temporal_cascade(spike_times: Dict[str, float], 
                         causal_links: List[Tuple[str, str, float]],
                         output_path: str,
                         top_n_links: int = 20):
    """
    Generates the 'Infection Cascade' visualization for the symposium.
    
    Args:
        spike_times: Dictionary {Gene: Onset_Time}
        causal_links: List of (Source, Target, Weight)
        output_path: File path to save the plot.
        top_n_links: Number of strongest links to visualize.
    """
    # Sort genes by onset time
    sorted_genes = sorted(spike_times.keys(), key=lambda g: spike_times[g])
    gene_rank = {g: i for i, g in enumerate(sorted_genes)}
    
    times = [spike_times[g] for g in sorted_genes]
    ranks = list(range(len(sorted_genes)))
    
    plt.figure(figsize=(12, 8))
    
    # Plot gene onset events
    plt.scatter(times, ranks, c='darkred', s=50, zorder=3, label='Sudden Onset Event')
    
    # Draw causal arrows
    # Filter for top links that exist in our spike_times set
    valid_links = [
        (u, v, w) for u, v, w in causal_links 
        if u in spike_times and v in spike_times
    ]
    # Sort by weight
    valid_links.sort(key=lambda x: x[2], reverse=True)
    
    print(f"Plotting top {top_n_links} temporal links out of {len(valid_links)} candidates...")
    
    for u, v, w in valid_links[:top_n_links]:
        t1, r1 = spike_times[u], gene_rank[u]
        t2, r2 = spike_times[v], gene_rank[v]
        
        # Draw arrow
        plt.arrow(t1, r1, t2-t1, r2-r1, 
                  head_width=0.2, head_length=0.02 * (max(times)-min(times)), 
                  fc='gray', ec='gray', alpha=0.6, length_includes_head=True)
        
    # Formatting
    plt.yticks(ranks, sorted_genes, fontsize=8)
    plt.xlabel("Time (Pseudo-time / Hours Post-Infection)")
    plt.ylabel("Host Factors (Sorted by Onset)")
    plt.title("Temporal Cascade of Host Immune Response")
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Cascade plot saved to {output_path}")
