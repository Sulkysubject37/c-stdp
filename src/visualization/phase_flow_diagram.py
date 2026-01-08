import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import os

def plot_phase_flow(matrix_file, output_dir):
    print("Generating Phase Flow Diagram...")
    
    if not os.path.exists(matrix_file):
        print("Matrix file not found.")
        return
        
    df = pd.read_csv(matrix_file, index_col=0)
    
    # Build Graph
    G = nx.DiGraph()
    phases = ['Phase I', 'Phase II', 'Phase III', 'Phase IV']
    
    for p in phases:
        G.add_node(p)
        
    for u in df.index:
        for v in df.columns:
            w = df.loc[u, v]
            if w > 0.1: # Threshold for visibility
                G.add_edge(u, v, weight=w)
                
    # Layout: Linear or Circular?
    # Linear: I -> II -> III -> IV
    pos = {
        'Phase I': (0, 0),
        'Phase II': (1, 0),
        'Phase III': (2, 0),
        'Phase IV': (3, 0)
    }
    
    fig, ax = plt.subplots(figsize=(10, 4))
    
    # Draw Nodes
    nx.draw_networkx_nodes(G, pos, node_size=2000, node_color='lightgray', edgecolors='black', ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', ax=ax)
    
    # Draw Edges
    weights = [G[u][v]['weight'] for u, v in G.edges]
    # Normalize width
    if weights:
        max_w = max(weights)
        widths = [(w / max_w) * 5 + 1 for w in weights]
    else:
        widths = []
        
    nx.draw_networkx_edges(G, pos, width=widths, arrowsize=20, connectionstyle='arc3,rad=-0.3', ax=ax)
    
    # Annotation
    ax.set_title("Phase-to-Phase Temporal Influence Flow")
    ax.text(1.5, -0.5, "Edges represent aggregated temporal precedence weights.", ha='center', style='italic')
    
    if not G.edges:
        ax.text(1.5, 0.5, "No inter-phase links detected in top 200 regulators.", ha='center', color='red')
    
    ax.axis('off')
    
    output_path = os.path.join(output_dir, "phase_flow.png")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved flow diagram to {output_path}")

if __name__ == "__main__":
    plot_phase_flow(
        "visuals/symposium_final/phase_influence.csv",
        "visuals/symposium_final"
    )
