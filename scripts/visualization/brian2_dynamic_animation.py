import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.animation as animation
import matplotlib.cm as cm
import matplotlib.colors as mcolors

import matplotlib
matplotlib.use('Agg')

def run_brian2_advanced_animations():
    print("--- Starting Advanced Brian2 Dynamic Animations ---")
    
    # 1. Load Data
    npz_path = "results/brian2/simulation_data.npz"
    if not os.path.exists(npz_path):
        return
        
    data = np.load(npz_path, allow_pickle=True)
    weights = data['weights']
    genes = data['genes']
    n_genes = len(genes)
    
    # Collect spikes
    spikes_dict = {} # gene_idx -> list of times
    all_times = []
    
    for i, gene in enumerate(genes):
        spikes = data[gene]
        times_ms = spikes * 10.0 # 10ms steps
        spikes_dict[i] = times_ms
        all_times.extend(times_ms)
        
    duration_ms = int(max(all_times)) + 50 if all_times else 100
    
    # Setup Output
    output_dir = "results/visuals/brian2"
    os.makedirs(output_dir, exist_ok=True)
    
    # --- Animation 1: Network Pulse (Signal Propagation) ---
    print("Generating Animation 1: Network Pulse...")
    
    G_nx = nx.DiGraph()
    for i in range(n_genes):
        G_nx.add_node(i, label=genes[i])
        for j in range(n_genes):
            if weights[i, j] > 0:
                G_nx.add_edge(i, j, weight=weights[i, j])
    
    # Layout: Hierarchical (using dot if available, else spring)
    try:
        pos = nx.nx_agraph.graphviz_layout(G_nx, prog='dot')
    except:
        pos = nx.spring_layout(G_nx, seed=42)
        
    fig, ax = plt.subplots(figsize=(10, 10))
    # Scientific style
    plt.style.use('seaborn-v0_8-paper')
    
    # Node sizes based on Out-Degree
    out_deg = dict(G_nx.out_degree(weight='weight'))
    sizes = [300 + out_deg[i]*500 for i in range(n_genes)]
    
    def update_pulse(frame):
        ax.clear()
        current_time = frame * 10.0
        
        # Node Colors: "Cool" (inactive) -> "Hot" (Spiking)
        # Decay function for visual persistence
        node_colors = []
        for i in range(n_genes):
            # Find most recent spike
            last_spike = -999
            for t in spikes_dict[i]:
                if t <= current_time:
                    last_spike = max(last_spike, t)
            
            dt = current_time - last_spike
            # Flash intensity: exp decay over 50ms
            intensity = np.exp(-dt / 50.0) if dt >= 0 else 0
            # Colormap: Viridis (Dark Blue -> Yellow)
            # Or Plasma (Purple -> Orange) -> Fire-like
            color = cm.plasma(intensity)
            node_colors.append(color)
            
        # Draw Edges first
        nx.draw_networkx_edges(G_nx, pos, ax=ax, arrowstyle='->', arrowsize=15, 
                               edge_color='gray', alpha=0.3, width=1.5)
        
        # Draw Nodes
        nx.draw_networkx_nodes(G_nx, pos, ax=ax, node_color=node_colors, 
                               node_size=sizes, edgecolors='white', linewidths=1.5)
        
        # Labels
        nx.draw_networkx_labels(G_nx, pos, ax=ax, labels={i:genes[i].split('.')[0] for i in range(n_genes)}, 
                                font_size=9, font_weight='bold', font_color='black')
        
        ax.set_title(f"Network Activity (Pulse) | T={current_time:.0f} ms", fontsize=14)
        ax.axis('off')
        
    ani_pulse = animation.FuncAnimation(fig, update_pulse, frames=200, interval=50)
    ani_pulse.save(f"{output_dir}/anim_1_network_pulse.gif", writer='pillow', fps=15)
    plt.close()

    # --- Animation 2: Cumulative Activation (Heatmap Accumulation) ---
    print("Generating Animation 2: Cumulative Activation...")
    
    fig2, ax2 = plt.subplots(figsize=(10, 10))
    
    def update_accum(frame):
        ax2.clear()
        current_time = frame * 10.0
        
        # Node Color based on Total Spikes so far
        node_colors = []
        counts = []
        for i in range(n_genes):
            count = sum(1 for t in spikes_dict[i] if t <= current_time)
            counts.append(count)
            
        # Normalize
        max_c = max(counts) if max(counts) > 0 else 1
        norm = mcolors.Normalize(vmin=0, vmax=max_c)
        cmap = cm.viridis
        
        node_colors = [cmap(norm(c)) for c in counts]
        
        nx.draw_networkx_edges(G_nx, pos, ax=ax2, edge_color='gray', alpha=0.3)
        nx.draw_networkx_nodes(G_nx, pos, ax=ax2, node_color=node_colors, node_size=sizes)
        nx.draw_networkx_labels(G_nx, pos, ax=ax2, labels={i:genes[i].split('.')[0] for i in range(n_genes)}, font_size=9)
        
        # Add colorbar only once? Hard in FuncAnimation. 
        # We'll rely on visual intensity.
        
        ax2.set_title(f"Cumulative Activation | T={current_time:.0f} ms", fontsize=14)
        ax2.axis('off')

    ani_accum = animation.FuncAnimation(fig2, update_accum, frames=200, interval=50)
    ani_accum.save(f"{output_dir}/anim_2_cumulative.gif", writer='pillow', fps=15)
    plt.close()
    
    print("Animations saved.")

if __name__ == "__main__":
    run_brian2_advanced_animations()