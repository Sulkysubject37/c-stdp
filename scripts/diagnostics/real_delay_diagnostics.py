import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.getcwd())

from src.cstdp.utils.spike_encoding import calculate_adaptive_thresholds
from src.cstdp.core import CausalSTDP

def run_delay_diagnostics():
    print("--- Starting Delay Structure Diagnostics (GSE215865) ---")
    
    # 1. Load Data
    input_file = "data/processed/GSE215865_subset.csv"
    df = pd.read_csv(input_file, index_col=0)
    data = df.values
    n_genes, n_timepoints = data.shape
    time_points = np.arange(n_timepoints)
    
    # 2. Run Inference with Trace
    thresholds = calculate_adaptive_thresholds(data, sigma=1.5)
    cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
    spike_trains = cstdp.compute_spike_times(data, time_points, thresholds)
    
    weights, trace = cstdp.run_cstdp(spike_trains, n_genes, return_trace=True)
    
    # 3. Filter for Strong Edges and Compute Weighted Delay
    thresh = 0.3
    strong_edges_indices = np.argwhere(weights > thresh)
    
    weighted_delays = []
    
    print(f"Analyzing {len(strong_edges_indices)} strong edges...")
    
    for src, tgt in strong_edges_indices:
        events = trace.get((src, tgt), [])
        if not events:
            continue
            
        # Events are (time, dw). We need dt. 
        # The trace stores 'time' of the second spike (max(t_i, t_j)).
        # But we don't have t_i stored directly in trace.
        # However, for Potentiation (dw > 0), dt is related to dw:
        # dw = A * exp(-dt/tau) => log(dw/A) = -dt/tau => dt = -tau * log(dw/A)
        
        # Invert the STDP formula to recover dt
        # Valid only for positive dw
        
        total_dt_weight = 0.0
        total_w = 0.0
        
        for t, dw in events:
            if dw > 1e-9: # Positive contribution
                # dt = -tau * ln(dw / A_pos)
                try:
                    dt = -10.0 * np.log(dw / 0.05)
                    if dt > 0:
                        total_dt_weight += dt * dw
                        total_w += dw
                except:
                    pass
                    
        if total_w > 0:
            weighted_delays.append(total_dt_weight / total_w)
            
    inferred_delays = np.array(weighted_delays)
    
    print(f"Number of analyzed edges: {len(inferred_delays)}")
    if len(inferred_delays) > 0:
        print(f"Mean Weighted Delay: {np.mean(inferred_delays):.2f} steps")
        print(f"Median Weighted Delay: {np.median(inferred_delays):.2f} steps")
    
    # 4. Visuals
    os.makedirs("results/visuals", exist_ok=True)
    plt.figure(figsize=(10, 6))
    sns.histplot(inferred_delays, bins=30, kde=True, color='purple')
    plt.title("Inferred Temporal Delays (Delta t) for Strong Edges")
    plt.xlabel("Delay (Pseudo-Time Steps)")
    plt.ylabel("Count")
    plt.axvline(x=10, color='r', linestyle='--', label="Tau_pos (10)")
    plt.legend()
    plt.savefig("results/visuals/real_delay_distribution.png")
    plt.close()
    
    print("Interpretation: Delta t reflects ordering preference under the chosen temporal discretization.")

if __name__ == "__main__":
    run_delay_diagnostics()
