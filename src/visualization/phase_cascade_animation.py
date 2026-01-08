import matplotlib.pyplot as plt
import matplotlib.animation as animation
import networkx as nx
import pandas as pd
import os

def animate_phase_cascade(matrix_file, output_dir):
    print("Generating Phase Cascade Animation...")
    
    if not os.path.exists(matrix_file):
        print("Matrix file not found.")
        return
        
    df = pd.read_csv(matrix_file, index_col=0)
    
    # Graph Setup
    G = nx.DiGraph()
    phases = ['Phase I', 'Phase II', 'Phase III', 'Phase IV']
    for p in phases:
        G.add_node(p)
    
    for u in df.index:
        for v in df.columns:
            w = df.loc[u, v]
            if w > 0.1:
                G.add_edge(u, v, weight=w)
                
    pos = {
        'Phase I': (0, 0),
        'Phase II': (1, 0),
        'Phase III': (2, 0),
        'Phase IV': (3, 0)
    }
    
    fig, ax = plt.subplots(figsize=(10, 4))
    
    # Timing (Conceptual Percentile)
    activation_thresholds = {
        'Phase I': 0,
        'Phase II': 20,
        'Phase III': 50,
        'Phase IV': 80
    }
    
    colors = {
        'Phase I': '#1f77b4',
        'Phase II': '#ff7f0e',
        'Phase III': '#2ca02c',
        'Phase IV': '#d62728'
    }

    def update(frame):
        ax.clear()
        current_time = frame
        
        # Nodes
        node_colors = []
        for p in phases:
            if current_time >= activation_thresholds[p]:
                node_colors.append(colors[p])
            else:
                node_colors.append('lightgray')
                
        nx.draw_networkx_nodes(G, pos, node_size=2000, node_color=node_colors, edgecolors='black', ax=ax)
        nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', ax=ax)
        
        # Edges (Only show if source and target are active)
        active_edges = []
        widths = []
        for u, v in G.edges:
            if current_time >= activation_thresholds[u] and current_time >= activation_thresholds[v]:
                active_edges.append((u, v))
                widths.append(G[u][v]['weight'] / df.values.max() * 5 + 1)
                
        if active_edges:
             nx.draw_networkx_edges(G, pos, edgelist=active_edges, width=widths, arrowsize=20, connectionstyle='arc3,rad=-0.3', ax=ax)
        
        ax.set_title(f"Phase Activation Cascade (T={current_time}%)")
        ax.axis('off')

    ani = animation.FuncAnimation(fig, update, frames=range(0, 101, 2), interval=100)
    
    mp4_path = os.path.join(output_dir, "phase_cascade.mp4")
    ani.save(mp4_path, writer='ffmpeg', fps=20)
    print(f"Saved animation to {mp4_path}")

if __name__ == "__main__":
    animate_phase_cascade(
        "visuals/symposium_final/phase_influence.csv",
        "visuals/symposium_final"
    )
