import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from brian2 import *

# Ensure we use the non-GUI backend for matplotlib if needed, 
# but Brian2 often handles this. We'll set it just in case.
import matplotlib
matplotlib.use('Agg')

def run_brian2_viz():
    print("--- Starting Brian2 Visualization (No Learning) ---")
    
    # 1. Load Data
    npz_path = "results/brian2/simulation_data.npz"
    if not os.path.exists(npz_path):
        print("Simulation data not found.")
        return
        
    data = np.load(npz_path, allow_pickle=True)
    weights = data['weights']
    genes = data['genes']
    n_genes = len(genes)
    
    # Collect all spike indices and times for SpikeGeneratorGroup
    all_indices = []
    all_times = []
    
    print("Reconstructing spike trains...")
    for i, gene in enumerate(genes):
        spikes = data[gene]
        # Spikes are in "sample index" units. We treat 1 sample = 10 ms
        # to make it biologically visible in Brian2
        time_scale = 10 * ms
        if len(spikes) > 0:
            all_indices.extend([i] * len(spikes))
            all_times.extend(spikes * time_scale)
            
    # Sort for Brian2 efficiency
    sort_idx = np.argsort(all_times)
    all_indices = np.array(all_indices)[sort_idx]
    all_times = np.array(all_times)[sort_idx] * second # Ensure units
    
    duration = (max(all_times) if len(all_times) > 0 else 0*ms) + 50*ms
    
    # 2. Setup Brian2 Model
    start_scope()
    
    # Input group (forcing the spikes)
    G = SpikeGeneratorGroup(n_genes, all_indices, all_times)
    
    # Target Neurons (Leaky Integrate and Fire to visualize propagation)
    # tau_m = 10ms (membrane time constant)
    # v_rest = 0
    # v_thresh = 10 (arbitrary high to avoid self-spiking unless driven, 
    # but here we just want to see potential fluctuations driven by inputs)
    # Actually, we want to see the "Network State".
    # Since inputs are "forced", we can just monitor the "Post-Synaptic Potential"
    # accumulated by a ghost neuron layer connected to the input layer.
    
    eqs = '''
    dv/dt = -v / (10*ms) : 1
    '''
    
    # We simulate a "Ghost" layer that receives input from G
    # to visualize the effect of the weights.
    Ghost = NeuronGroup(n_genes, eqs, method='exact')
    
    # Synapses with FIXED weights
    # Connect Input G -> Ghost
    # The weights matrix is [Source, Target]
    S = Synapses(G, Ghost, model='w : 1', on_pre='v += w')
    
    sources, targets = np.where(weights > 0)
    S.connect(i=sources, j=targets)
    
    # Assign weights (normalize to be visible, e.g. max 1.0 -> 1.0mV boost)
    # Our weights are ~0.3 to 1.0. 
    # v threshold is arbitrary. Let's say v is normalized units.
    # We need to map the numpy matrix to the Synapses object.
    # S.w[i] accesses the weight of the i-th synapse in the flattened array of connections.
    # We iterate manually to be safe.
    
    print("Setting fixed synaptic weights...")
    for k in range(len(sources)):
        src = sources[k]
        tgt = targets[k]
        w_val = weights[src, tgt]
        # Find the specific synapse index
        # Brian2 indexing can be tricky. S.w['i==src and j==tgt'] = w_val is cleaner.
        S.w[f'i=={src} and j=={tgt}'] = w_val
        
    # Monitors
    M = StateMonitor(Ghost, 'v', record=True)
    # We also monitor the input spikes
    SM = SpikeMonitor(G)
    
    # 3. Run
    print(f"Running simulation for {duration}...")
    run(duration)
    
    # 4. Plotting
    output_dir = "results/visuals/brian2"
    os.makedirs(output_dir, exist_ok=True)
    
    # Plot 1: Raster
    plt.figure(figsize=(12, 6))
    plt.plot(SM.t/ms, SM.i, '.k')
    plt.yticks(range(n_genes), genes, fontsize=8)
    plt.xlabel('Time (ms)')
    plt.ylabel('Gene Index')
    plt.title('Replayed Spike Raster (Input)')
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{output_dir}/raster.png")
    plt.close()
    
    # Plot 2: Membrane Potentials (Causal Integration)
    # Show traces for the 3 targets identified in Step 1
    # We need to pick them. We know the indices.
    # Let's pick 3 genes with high in-degree in the mini-adj
    in_degrees = np.sum(weights > 0, axis=0)
    target_indices = np.argsort(in_degrees)[::-1][:3]
    
    plt.figure(figsize=(12, 6))
    for idx in target_indices:
        plt.plot(M.t/ms, M.v[idx], label=f"{genes[idx]} (In-Deg: {in_degrees[idx]})")
    
    plt.xlabel('Time (ms)')
    plt.ylabel('Integrated Potential (v)')
    plt.title('Causal Integration (Virtual Membrane Potential)')
    plt.legend()
    plt.savefig(f"{output_dir}/membrane_traces.png")
    plt.close()
    
    # Plot 3: Static Graph
    G_nx = nx.DiGraph()
    for i in range(n_genes):
        for j in range(n_genes):
            if weights[i, j] > 0:
                G_nx.add_edge(genes[i], genes[j], weight=weights[i, j])
                
    plt.figure(figsize=(8, 8))
    pos = nx.circular_layout(G_nx)
    nx.draw_networkx_nodes(G_nx, pos, node_color='lightblue', node_size=500)
    nx.draw_networkx_labels(G_nx, pos, font_size=8)
    edges = nx.draw_networkx_edges(G_nx, pos, arrowstyle='->', arrowsize=15, edge_color='gray')
    plt.title("Simulated Subgraph Structure")
    plt.savefig(f"{output_dir}/static_graph.png")
    plt.close()
    
    print(f"Visuals saved to {output_dir}")

if __name__ == "__main__":
    run_brian2_viz()
