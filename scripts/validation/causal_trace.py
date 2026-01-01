import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.getcwd())

from src.cstdp.utils.simulate_grn import generate_synthetic_grn, simulate_expression
from src.cstdp.utils.spike_encoding import calculate_adaptive_thresholds
from src.cstdp.core import CausalSTDP

def run_causal_trace():
    print("--- Starting Edge-Level Causal Trace Extraction ---")
    
    # 1. Setup Data
    np.random.seed(42)
    n_genes = 10
    n_timepoints = 1000
    true_adj, delays = generate_synthetic_grn(n_genes, connection_prob=0.2, seed=42)
    expression = simulate_expression(n_genes, n_timepoints, true_adj, delays, seed=42)
    time_points = np.arange(n_timepoints)
    
    # 2. Run Inference with Trace
    thresholds = calculate_adaptive_thresholds(expression, sigma=1.5)
    cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
    spike_trains = cstdp.compute_spike_times(expression, time_points, thresholds)
    
    weights, trace = cstdp.run_cstdp(spike_trains, n_genes, return_trace=True)
    
    # 3. Select an Edge to Trace
    # Let's find a True Positive (TP) edge
    sources, targets = np.where(true_adj > 0)
    if len(sources) == 0:
        print("No true edges found to trace.")
        return
        
    src, tgt = sources[0], targets[0]
    edge_trace = trace[(src, tgt)]
    
    print(f"Tracing Edge: Gene {src} -> Gene {tgt}")
    print(f"Number of contributing spike pairs: {len(edge_trace)}")
    
    # 4. Generate Table
    trace_df = pd.DataFrame(edge_trace, columns=["Time", "Delta_W"])
    trace_df["Cumulative_W"] = trace_df["Delta_W"].cumsum().clip(0, 1.0)
    trace_df.to_csv("analysis/causal_trace_data.csv", index=False)
    
    # 5. Plot
    plt.figure(figsize=(12, 6))
    plt.step(trace_df["Time"], trace_df["Cumulative_W"], where='post', label="Cumulative Weight")
    plt.scatter(trace_df["Time"], [0]*len(trace_df), marker='|', color='r', alpha=0.5, label="Spike Event")
    plt.axhline(y=weights[src, tgt], color='g', linestyle='--', alpha=0.5, label="Final Weight")
    plt.title(f"Causal Weight Trace: Gene {src} -> Gene {tgt}")
    plt.xlabel("Time (Steps)")
    plt.ylabel("Weight Influence")
    plt.legend()
    plt.grid(True, alpha=0.2)
    plt.savefig("analysis/visuals/causal_trace_plot.png")
    plt.close()
    
    print("Interpretation: Edge strength grows with each causal spike pair (Pre before Post).")

if __name__ == "__main__":
    run_causal_trace()
