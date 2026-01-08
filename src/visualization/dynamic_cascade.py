import matplotlib.pyplot as plt
import matplotlib.animation as animation
import networkx as nx
import pandas as pd
import numpy as np
import os
from scipy.stats import zscore

def generate_dynamic_cascade(dataset_name, input_file, adj_file, output_dir):
    print(f"Generating Dynamic Cascade for {dataset_name}...")
    
    # 1. Load Data & Logic (Replicated for consistency)
    df = pd.read_csv(input_file, index_col=0)
    data = np.nan_to_num(df.values, nan=0.0)
    genes = df.index.values
    time_points = np.arange(data.shape[1])
    
    derivatives = np.gradient(data, axis=1)
    norm_derivatives = np.nan_to_num(zscore(derivatives, axis=1), nan=0.0)
    spike_mask = norm_derivatives > 1.5
    spike_counts = np.sum(spike_mask, axis=1)
    
    SUBSET_SIZE = 30 # Reduced for graph clarity (avoid hairball)
    top_indices = np.argsort(spike_counts)[-SUBSET_SIZE:]
    
    genes_sub = genes[top_indices]
    spike_mask_sub = spike_mask[top_indices]
    
    onset_times = {}
    for i, gene in enumerate(genes_sub):
        spikes = time_points[spike_mask_sub[i, :]]
        if len(spikes) > 0:
            onset_times[gene] = float(spikes[0])
        else:
            onset_times[gene] = float(time_points[-1])
            
    # 2. Load Adjacency (Subset to top 30)
    full_adj = pd.read_csv(adj_file, index_col=0)
    # Filter adj to just genes_sub
    # Intersect
    valid_genes = [g for g in genes_sub if g in full_adj.index]
    adj_sub = full_adj.loc[valid_genes, valid_genes]
    
    # Build Graph
    G = nx.DiGraph()
    for g in valid_genes:
        G.add_node(g, onset=onset_times[g])
        
    for u in valid_genes:
        for v in valid_genes:
            w = adj_sub.loc[u, v]
            if w > 0.1: # Visibility threshold
                G.add_edge(u, v, weight=w)
                
    # Layout (Deterministic)
    pos = nx.spring_layout(G, k=0.5, seed=42)
    
    # 3. Animation
    fig, ax = plt.subplots(figsize=(10, 10))
    
    def update(frame):
        ax.clear()
        current_time = frame
        
        # Nodes: Active vs Dormant
        active_nodes = [n for n in G.nodes if onset_times[n] <= current_time]
        dormant_nodes = [n for n in G.nodes if onset_times[n] > current_time]
        
        # Draw Dormant (Gray, Ghost)
        nx.draw_networkx_nodes(G, pos, nodelist=dormant_nodes, node_color='lightgray', alpha=0.3, node_size=100, ax=ax)
        
        # Draw Active (Red, Bold)
        if active_nodes:
            nx.draw_networkx_nodes(G, pos, nodelist=active_nodes, node_color='darkred', alpha=1.0, node_size=300, ax=ax)
            nx.draw_networkx_labels(G, pos, labels={n:n for n in active_nodes}, font_size=10, ax=ax)
            
        # Edges: Only if Source is Active
        active_edges = [(u, v) for u, v in G.edges if onset_times[u] <= current_time]
        
        if active_edges:
            nx.draw_networkx_edges(G, pos, edgelist=active_edges, edge_color='black', alpha=0.5, arrows=True, ax=ax)
            
        ax.set_title(f"Infection Cascade: {dataset_name} (T={current_time:.1f})")
        ax.axis('off')

    max_t = max(onset_times.values()) + 5
    ani = animation.FuncAnimation(fig, update, frames=np.linspace(0, max_t, 100), interval=50)
    
    mp4_path = os.path.join(output_dir, "dynamic_cascade.mp4")
    ani.save(mp4_path, writer='ffmpeg', fps=20)
    print(f"Saved cascade to {mp4_path}")
    
    # Static Snapshot
    update(max_t)
    plt.savefig(os.path.join(output_dir, "dynamic_cascade_static.png"), dpi=150)
    plt.close()

if __name__ == "__main__":
    generate_dynamic_cascade(
        "GSE215865",
        "data/processed/GSE215865_immune_subset.csv",
        "visuals/real_data/GSE215865_Immune/inferred_grn_adj.csv",
        "visuals/symposium/GSE215865"
    )
    generate_dynamic_cascade(
        "GSE157859",
        "data/processed/GSE157859_immune_subset.csv",
        "visuals/real_data/GSE157859_Immune/inferred_grn_adj.csv",
        "visuals/symposium/GSE157859"
    )
