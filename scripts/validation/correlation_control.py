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
from src.cstdp.utils.evaluate import calculate_metrics

def run_correlation_control():
    print("--- Starting Competing Explanation Control ---")
    
    # 1. Setup Data
    np.random.seed(42)
    n_genes = 15
    n_timepoints = 1000
    true_adj, delays = generate_synthetic_grn(n_genes, connection_prob=0.15, seed=42)
    expression = simulate_expression(n_genes, n_timepoints, true_adj, delays, seed=42)
    
    # 2. C-STDP Inference
    thresholds = calculate_adaptive_thresholds(expression, sigma=1.5)
    cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
    spike_trains = cstdp.compute_spike_times(expression, np.arange(n_timepoints), thresholds)
    weights_stdp = cstdp.run_cstdp(spike_trains, n_genes)
    
    # Normalize STDP
    max_stdp = np.max(weights_stdp)
    norm_stdp = weights_stdp / max_stdp if max_stdp > 0 else weights_stdp
    
    # 3. Time-Lagged Correlation (Baseline)
    # Correlation between Gene i (t) and Gene j (t+lag)
    # We use lag = 10 (avg delay in simulation)
    lag = 10
    corr_matrix = np.zeros((n_genes, n_genes))
    
    for i in range(n_genes):
        for j in range(n_genes):
            if i == j: continue
            # Calculate correlation with lag
            # Gene i is "cause" (earlier), Gene j is "effect" (later)
            series_i = expression[i, :-lag]
            series_j = expression[j, lag:]
            c = np.corrcoef(series_i, series_j)[0, 1]
            corr_matrix[i, j] = max(0, c) # Only positive correlations
            
    # Normalize Correlation
    max_corr = np.max(corr_matrix)
    norm_corr = corr_matrix / max_corr if max_corr > 0 else corr_matrix
    
    # 4. Compare Edges
    thresh = 0.3
    edges_stdp = norm_stdp > thresh
    edges_corr = norm_corr > thresh
    
    stdp_only = edges_stdp & ~edges_corr
    corr_only = edges_corr & ~edges_stdp
    both = edges_stdp & edges_corr
    
    print("\n--- COMPARISON ---")
    print(f"Edges by STDP:    {np.sum(edges_stdp)}")
    print(f"Edges by Correlation: {np.sum(edges_corr)}")
    print(f"Intersection:     {np.sum(both)}")
    print(f"STDP Only:        {np.sum(stdp_only)}")
    print(f"Correlation Only: {np.sum(corr_only)}")
    
    # Accuracy vs Truth
    m_stdp = calculate_metrics(true_adj, norm_stdp, threshold=thresh)
    m_corr = calculate_metrics(true_adj, norm_corr, threshold=thresh)
    
    print(f"\nSTDP Precision:    {m_stdp['precision']:.4f}")
    print(f"Corr Precision:    {m_corr['precision']:.4f}")
    
    # 5. Visuals
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sns.heatmap(norm_stdp, ax=axes[0], cmap="viridis", vmin=0, vmax=1)
    axes[0].set_title("C-STDP Weights")
    sns.heatmap(norm_corr, ax=axes[1], cmap="viridis", vmin=0, vmax=1)
    axes[1].set_title(f"Time-Lagged Correlation (Lag={lag})")
    plt.savefig("results/visuals/stdp_vs_correlation.png")
    plt.close()
    
    print("\nInterpretation: STDP filters edges by event-precedence, which is stricter than lagged correlation.")

if __name__ == "__main__":
    run_correlation_control()
