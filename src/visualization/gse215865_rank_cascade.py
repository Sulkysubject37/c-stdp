import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os

def plot_rank_cascade(rank_file, adj_file, output_dir):
    print("Generating Rank-Based Cascade for GSE215865...")
    
    # 1. Load Ranks & Adj
    ranks = pd.read_csv(rank_file)
    rank_map = pd.Series(ranks.Percentile.values, index=ranks.Gene).to_dict()
    
    adj = pd.read_csv(adj_file, index_col=0)
    
    # 2. Extract Top 20 Strongest Temporal Links
    # Flatten adj
    links = []
    for u in adj.index:
        for v in adj.columns:
            w = adj.loc[u, v]
            if w > 0.0:
                links.append((u, v, w))
                
    links.sort(key=lambda x: x[2], reverse=True)
    top_links = links[:20]
    
    # 3. Build Graph
    G = nx.DiGraph()
    
    layers = {'Early': [], 'Inter': [], 'Late': []}
    
    for u, v, w in top_links:
        # Determine Layer
        def get_layer(gene):
            p = rank_map.get(gene, -1)
            if p == -1: return None
            if p <= 20: return 'Early'
            if p <= 60: return 'Inter'
            return 'Late'
            
        l_u = get_layer(u)
        l_v = get_layer(v)
        
        if l_u and l_v:
            G.add_edge(u, v, weight=w)
            if u not in layers[l_u]: layers[l_u].append(u)
            if v not in layers[l_v]: layers[l_v].append(v)
            
    # 4. Plot with Multipartite Layout
    pos = {}
    
    # Manually position layers (X-axis) and disperse nodes (Y-axis)
    layer_x = {'Early': 0, 'Inter': 1, 'Late': 2}
    
    for layer, nodes in layers.items():
        x = layer_x[layer]
        for i, node in enumerate(nodes):
            # Center Y
            y = (i - len(nodes)/2) 
            pos[node] = (x, y)
            
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Draw Nodes
    for layer, nodes in layers.items():
        color = {'Early': '#1f77b4', 'Inter': '#ff7f0e', 'Late': '#2ca02c'}[layer]
        nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=color, node_size=500, alpha=0.9, ax=ax)
        
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold', ax=ax)
    
    # Draw Edges (Curved?)
    nx.draw_networkx_edges(G, pos, edge_color='gray', width=1.5, arrowsize=20, ax=ax, connectionstyle='arc3,rad=0.1')
    
    # Annotations
    ax.text(0, min([p[1] for p in pos.values()]) - 1, "Early Layer\n(0-20%)", ha='center', fontweight='bold')
    ax.text(1, min([p[1] for p in pos.values()]) - 1, "Intermediate\n(20-60%)", ha='center', fontweight='bold')
    ax.text(2, min([p[1] for p in pos.values()]) - 1, "Late Layer\n(60-100%)", ha='center', fontweight='bold')
    
    ax.set_title("Simplified Directional Cascade (Rank-Ordered)")
    ax.axis('off')
    
    output_path = os.path.join(output_dir, "GSE215865_rank_cascade.png")
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved cascade to {output_path}")

if __name__ == "__main__":
    plot_rank_cascade(
        "visuals/symposium/GSE215865/activation_ranks.csv",
        "visuals/real_data/GSE215865_Immune/inferred_grn_adj.csv",
        "visuals/symposium/GSE215865"
    )
