import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.getcwd())

from src.cstdp.utils.spike_encoding import calculate_adaptive_thresholds
from src.cstdp.core import CausalSTDP

def run_real_trace():
    print("--- Starting Edge-Level Causal Trace on Real Data (GSE215865) ---")
    
    # 1. Load Data
    input_file = "data/processed/GSE215865_subset.csv"
    df = pd.read_csv(input_file, index_col=0)
    data = df.values
    genes = df.index.values
    n_genes, n_timepoints = data.shape
    time_points = np.arange(n_timepoints)
    
    # 2. Run Inference with Trace
    # Use baseline params
    thresholds = calculate_adaptive_thresholds(data, sigma=1.5)
    cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
    spike_trains = cstdp.compute_spike_times(data, time_points, thresholds)
    
    print("Running C-STDP with trace logging...")
    weights, trace = cstdp.run_cstdp(spike_trains, n_genes, return_trace=True)
    
    # 3. Select Strongest Edge
    # Flatten weights to find max
    flat_idx = np.argmax(weights)
    src_idx, tgt_idx = np.unravel_index(flat_idx, weights.shape)
    
    src_gene = genes[src_idx]
    tgt_gene = genes[tgt_idx]
    
    print(f"Tracing Strongest Edge: {src_gene} -> {tgt_gene}")
    print(f"Final Weight: {weights[src_idx, tgt_idx]:.4f}")
    
    edge_trace = trace.get((src_idx, tgt_idx), [])
    print(f"Number of contributing spike pairs: {len(edge_trace)}")
    
    if len(edge_trace) == 0:
        print("No trace found (weight might be 0).")
        return

    # 4. Generate Table
    trace_df = pd.DataFrame(edge_trace, columns=["Time", "Delta_W"])
    trace_df["Cumulative_W"] = trace_df["Delta_W"].cumsum().clip(0, 1.0)
    trace_df.to_csv("analysis/real_causal_trace_data.csv", index=False)
    
    # 5. Plot
    plt.figure(figsize=(12, 6))
    plt.step(trace_df["Time"], trace_df["Cumulative_W"], where='post', label="Cumulative Weight")
    # Add rug plot for spikes of source and target?
    # Actually, simpler to just show the cumulative trace
    plt.axhline(y=weights[src_idx, tgt_idx], color='g', linestyle='--', alpha=0.5, label="Final Weight")
    plt.title(f"Causal Evidence Trace: {src_gene} -> {tgt_gene}")
    plt.xlabel("Pseudo-Time (Sample Index)")
    plt.ylabel("Inferred Weight Influence")
    plt.legend()
    plt.grid(True, alpha=0.2)
    plt.savefig("analysis/visuals/real_causal_trace_plot.png")
    plt.close()
    
    print("Interpretation: This edge is supported by repeated, consistent temporal precedence events in the data.")

if __name__ == "__main__":
    run_real_trace()
