import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from brian2 import *

import matplotlib
matplotlib.use('Agg')

def run_membrane_animation():
    print("--- Starting Brian2 Membrane Trace Animation ---")
    
    # 1. Load Data
    npz_path = "results/brian2/simulation_data.npz"
    if not os.path.exists(npz_path):
        return
        
    data = np.load(npz_path, allow_pickle=True)
    weights = data['weights']
    genes = data['genes']
    n_genes = len(genes)
    
    # Reconstruct Spikes
    all_indices = []
    all_times = []
    for i, gene in enumerate(genes):
        spikes = data[gene]
        if len(spikes) > 0:
            all_indices.extend([i] * len(spikes))
            all_times.extend(spikes * 10 * ms)
            
    sort_idx = np.argsort(all_times)
    all_indices = np.array(all_indices)[sort_idx]
    all_times = np.array(all_times)[sort_idx] * second
    
    duration = (max(all_times) if len(all_times) > 0 else 0*ms) + 50*ms
    
    # 2. Run Simulation (Re-run to capture state)
    start_scope()
    G = SpikeGeneratorGroup(n_genes, all_indices, all_times)
    Ghost = NeuronGroup(n_genes, 'dv/dt = -v / (10*ms) : 1', method='exact')
    S = Synapses(G, Ghost, model='w : 1', on_pre='v += w')
    sources, targets = np.where(weights > 0)
    S.connect(i=sources, j=targets)
    for k in range(len(sources)):
        S.w[f'i=={sources[k]} and j=={targets[k]}'] = weights[sources[k], targets[k]]
        
    M = StateMonitor(Ghost, 'v', record=True)
    SM = SpikeMonitor(G)
    
    run(duration)
    
    # 3. Animate
    # We will show the top 3 genes with highest in-degree (integrators)
    in_degrees = np.sum(weights > 0, axis=0)
    target_indices = np.argsort(in_degrees)[::-1][:3]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    plt.style.use('dark_background') # Cool look for traces
    
    times = M.t/ms
    window_size = 200 # ms to show at once? No, let's scroll.
    
    # Lines
    lines = []
    colors = ['#FF9AA2', '#FFB7B2', '#FFDAC1', '#E2F0CB', '#B5EAD7', '#C7CEEA']
    for i, idx in enumerate(target_indices):
        line, = ax1.plot([], [], lw=2, color=colors[i % len(colors)], label=genes[idx].split('.')[0])
        lines.append(line)
        
    ax1.set_ylabel("Membrane Potential (a.u.)")
    ax1.legend(loc='upper right')
    ax1.set_ylim(-0.1, np.max(M.v[target_indices]) * 1.1)
    
    # Raster dots
    raster_scatter = ax2.scatter([], [], s=10, c='white', alpha=0.8)
    ax2.set_ylabel("Gene Index")
    ax2.set_xlabel("Time (ms)")
    ax2.set_ylim(-1, n_genes)
    
    def update(frame):
        # Frame corresponds to time index in M.t
        # Downsample for speed: step = 10 indices (assuming 0.1ms dt default? Brian2 default is 0.1ms)
        # Let's say frame is actual time index.
        # Speed up: skip
        current_idx = frame * 10
        if current_idx >= len(times): return
        
        current_t = times[current_idx]
        
        # Show window [0, current_t] or scrolling [current_t - window, current_t]
        # Let's do accumulating trace
        
        # Update traces
        for i, idx in enumerate(target_indices):
            lines[i].set_data(times[:current_idx], M.v[idx, :current_idx])
            
        # Update Raster
        # Find spikes up to current_t
        spike_mask = SM.t/ms <= current_t
        if np.any(spike_mask):
            ax2.collections[0].set_offsets(np.c_[SM.t[spike_mask]/ms, SM.i[spike_mask]])
            
        ax1.set_xlim(0, max(current_t, 10))
        
        return lines + [raster_scatter]
        
    # Frames: Total simulation steps / 10
    total_frames = min(500, len(times) // 10)
    
    print(f"Generating membrane animation ({total_frames} frames)...")
    ani = animation.FuncAnimation(fig, update, frames=total_frames, interval=30, blit=False)
    
    output_dir = "results/visuals/brian2"
    ani.save(f"{output_dir}/anim_3_membrane_traces.gif", writer='pillow', fps=20)
    print(f"Saved to {output_dir}/anim_3_membrane_traces.gif")

if __name__ == "__main__":
    run_membrane_animation()
