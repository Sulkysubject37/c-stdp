import matplotlib.pyplot as plt
import matplotlib.animation as animation
import networkx as nx
import pandas as pd
import numpy as np
import os
from scipy.stats import zscore

def generate_regulator_replay(dataset_name, regulator_gene_base, input_file, adj_file, output_dir):
    print(f"Generating Replay for Regulator {regulator_gene_base} in {dataset_name}...")
    
    # 1. Load Adjacency
    adj_df = pd.read_csv(adj_file, index_col=0)
    
    # Fuzzy match for regulator ID
    regulator_gene = None
    for gene in adj_df.index:
        if str(gene).startswith(regulator_gene_base):
            regulator_gene = gene
            break
            
    if not regulator_gene:
        print(f"Regulator {regulator_gene_base} not found in adjacency matrix.")
        return
        
    print(f"  Found full ID: {regulator_gene}")

    # Find downstream targets (direct children)
    # We want a subgraph: Regulator -> Targets -> Secondary Targets?
    # Just Depth 1 for clarity in replay
    targets = adj_df.loc[regulator_gene]
    targets = targets[targets > 0.01].index.tolist()
    
    if not targets:
        print(f"No strong targets found for {regulator_gene}.")
        return
        
    nodes = [regulator_gene] + targets
    
    # 2. Get Onset Times for these nodes
    df_expr = pd.read_csv(input_file, index_col=0)
    # subset
    valid_nodes = [n for n in nodes if n in df_expr.index]
    df_sub = df_expr.loc[valid_nodes]
    data = np.nan_to_num(df_sub.values, nan=0.0)
    
    derivatives = np.gradient(data, axis=1)
    norm_derivatives = np.nan_to_num(zscore(derivatives, axis=1), nan=0.0)
    spike_mask = norm_derivatives > 1.5
    
    time_points = np.arange(data.shape[1])
    onset_times = {}
    for i, gene in enumerate(df_sub.index):
        spikes = time_points[spike_mask[i, :]]
        if len(spikes) > 0:
            onset_times[gene] = float(spikes[0])
        else:
            onset_times[gene] = float(time_points[-1])
            
    # 3. Build Graph
    G = nx.DiGraph()
    for n in valid_nodes:
        G.add_node(n)
    
    # Edges from Regulator only
    for t in targets:
        if t in valid_nodes:
            w = adj_df.loc[regulator_gene, t]
            G.add_edge(regulator_gene, t, weight=w)
            
    pos = nx.spring_layout(G, seed=42)
    
    # 4. Animation
    fig, ax = plt.subplots(figsize=(8, 8))
    
    def update(frame):
        ax.clear()
        current_time = frame
        
        # Color nodes based on activation
        node_colors = []
        for n in G.nodes:
            if onset_times.get(n, 9999) <= current_time:
                if n == regulator_gene:
                    node_colors.append('gold') # Regulator active
                else:
                    node_colors.append('red') # Target active
            else:
                node_colors.append('lightgray') # Dormant
                
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=500, ax=ax)
        nx.draw_networkx_labels(G, pos, labels={n:n.split('.')[0] for n in G.nodes}, font_size=8, ax=ax)
        
        # Edges
        edge_colors = []
        for u, v in G.edges:
            if onset_times.get(u, 9999) <= current_time:
                edge_colors.append('black')
            else:
                edge_colors.append((0, 0, 0, 0)) # Invisible (tuple, not string)
                
        nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=2, ax=ax)
        
        ax.set_title(f"Regulator Replay: {regulator_gene.split('.')[0]} (T={current_time:.1f})")
        ax.axis('off')

    min_t = min(onset_times.values())
    max_t = max(onset_times.values()) + 5
    
    ani = animation.FuncAnimation(fig, update, frames=np.linspace(min_t-2, max_t, 60), interval=100)
    
    mp4_path = os.path.join(output_dir, "regulator_replay.mp4")
    ani.save(mp4_path, writer='ffmpeg', fps=10)
    print(f"Saved replay to {mp4_path}")

if __name__ == "__main__":
    generate_regulator_replay(
        "GSE215865",
        "ENSG00000186834", 
        "data/processed/GSE215865_immune_subset.csv",
        "visuals/real_data/GSE215865_Immune/inferred_grn_adj.csv",
        "visuals/symposium/GSE215865"
    )
    
    generate_regulator_replay(
        "GSE157859",
        "ENSG00000100292", 
        "data/processed/GSE157859_immune_subset.csv",
        "visuals/real_data/GSE157859_Immune/inferred_grn_adj.csv",
        "visuals/symposium/GSE157859"
    )
