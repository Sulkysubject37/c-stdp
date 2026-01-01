import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.animation as animation
from brian2 import *

import matplotlib
matplotlib.use('Agg') # Essential for headless environments

def run_brian2_animation():
    print("--- Starting Brian2 Dynamic Animation ---")
    
    # 1. Load Data
    npz_path = "results/brian2/simulation_data.npz"
    if not os.path.exists(npz_path):
        return
        
    data = np.load(npz_path, allow_pickle=True)
    weights = data['weights']
    genes = data['genes']
    n_genes = len(genes)
    
    # Collect spikes
    spikes_dict = {}
    all_times = []
    
    for i, gene in enumerate(genes):
        spikes = data[gene]
        # Store as simple list of times (ms)
        # Sample index * 10 ms
        times_ms = spikes * 10.0
        spikes_dict[i] = set(times_ms)
        all_times.extend(times_ms)
        
    duration_ms = int(max(all_times)) + 50 if all_times else 100
    
    # 2. Setup Graph Layout
    G_nx = nx.DiGraph()
    for i in range(n_genes):
        G_nx.add_node(i, label=genes[i])
        for j in range(n_genes):
            if weights[i, j] > 0:
                G_nx.add_edge(i, j, weight=weights[i, j])
                
    pos = nx.circular_layout(G_nx)
    
    # 3. Animation Setup
    fig, ax = plt.subplots(figsize=(8, 8))
    
    def update(frame):
        ax.clear()
        current_time = frame * 10.0 # 10ms steps per frame
        
        # Color nodes if they fired in this window (current_time +/- 5ms)
        node_colors = []
        for i in range(n_genes):
            fired = False
            for t in spikes_dict[i]:
                if abs(t - current_time) < 5.0:
                    fired = True
                    break
            node_colors.append('red' if fired else 'lightblue')
            
        # Draw
        nx.draw_networkx_nodes(G_nx, pos, ax=ax, node_color=node_colors, node_size=600)
        nx.draw_networkx_labels(G_nx, pos, ax=ax, labels={i:genes[i] for i in range(n_genes)}, font_size=8)
        
        # Edges (static for now, or could light up)
        nx.draw_networkx_edges(G_nx, pos, ax=ax, arrowstyle='->', arrowsize=15, edge_color='gray', alpha=0.5)
        
        ax.set_title(f"Time: {current_time:.0f} ms")
        ax.axis('off')
        
    # Limit frames to first 200 steps (2 seconds) to keep GIF size manageable
    num_frames = min(200, duration_ms // 10)
    
    print(f"Generating animation ({num_frames} frames)...")
    ani = animation.FuncAnimation(fig, update, frames=num_frames, interval=100)
    
    output_dir = "results/visuals/brian2"
    os.makedirs(output_dir, exist_ok=True)
    ani.save(f"{output_dir}/network_activity.gif", writer='pillow', fps=10)
    print(f"Animation saved to {output_dir}/network_activity.gif")

if __name__ == "__main__":
    run_brian2_animation()
