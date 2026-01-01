import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.tsa.stattools import grangercausalitytests

sys.path.append(os.getcwd())

from src.cstdp.utils.simulate_grn import generate_synthetic_grn, simulate_expression
from src.cstdp.utils.spike_encoding import calculate_adaptive_thresholds
from src.cstdp.core import CausalSTDP
from src.cstdp.utils.evaluate import calculate_metrics

def run_granger_baseline():
    print("--- Starting Granger Causality Baseline Comparison ---")
    
    # 1. Setup Data
    np.random.seed(42)
    n_genes = 10
    n_timepoints = 1000
    true_adj, delays = generate_synthetic_grn(n_genes, connection_prob=0.2, seed=42)
    expression = simulate_expression(n_genes, n_timepoints, true_adj, delays, seed=42)
    
    # 2. C-STDP Inference
    thresholds = calculate_adaptive_thresholds(expression, sigma=1.5)
    cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
    spike_trains = cstdp.compute_spike_times(expression, np.arange(n_timepoints), thresholds)
    weights_stdp = cstdp.run_cstdp(spike_trains, n_genes)
    
    # Normalize STDP
    max_stdp = np.max(weights_stdp)
    norm_stdp = weights_stdp / max_stdp if max_stdp > 0 else weights_stdp
    
    # 3. Granger Causality Inference
    # We test if i Granger-causes j
    # Use maxlag=10
    maxlag = 10
    granger_matrix = np.zeros((n_genes, n_genes))
    
    print("Running Pairwise Granger tests (this may take a moment)...")
    for i in range(n_genes):
        for j in range(n_genes):
            if i == j: continue
            # Data format: effects (j) in col 0, causes (i) in col 1
            data = expression[[j, i], :].T
            try:
                # verbose=False to keep logs clean
                res = grangercausalitytests(data, maxlag=maxlag, verbose=False)
                # Take minimum p-value across lags for the SSR F-test
                min_p = min([res[lag][0]['ssr_ftest'][1] for lag in range(1, maxlag+1)])
                # Store 1-p as "weight" (higher is more causal)
                granger_matrix[i, j] = 1.0 - min_p
            except:
                granger_matrix[i, j] = 0.0
                
    # Normalize Granger
    max_g = np.max(granger_matrix)
    norm_granger = granger_matrix / max_g if max_g > 0 else granger_matrix
    
    # 4. Metrics
    thresh = 0.3
    m_stdp = calculate_metrics(true_adj, norm_stdp, threshold=thresh)
    m_granger = calculate_metrics(true_adj, norm_granger, threshold=0.95) # Higher threshold for significance (1-p)
    
    print("\n--- BASELINE RESULTS ---")
    print(f"STDP F1 Score:    {m_stdp['f1']:.4f}")
    print(f"Granger F1 Score: {m_granger['f1']:.4f}")
    print(f"STDP Precision:   {m_stdp['precision']:.4f}")
    print(f"Granger Precision: {m_granger['precision']:.4f}")
    
    # 5. Visuals
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sns.heatmap(norm_stdp, ax=axes[0], cmap="viridis", vmin=0, vmax=1)
    axes[0].set_title("C-STDP Weights")
    sns.heatmap(norm_granger, ax=axes[1], cmap="magma", vmin=0, vmax=1)
    axes[1].set_title("Granger Causality (1-p)")
    plt.savefig("results/visuals/stdp_vs_granger.png")
    plt.close()
    
    # Save for summary
    summary = pd.DataFrame([
        {"Method": "STDP", "Precision": m_stdp['precision'], "Recall": m_stdp['recall'], "F1": m_stdp['f1']},
        {"Method": "Granger", "Precision": m_granger['precision'], "Recall": m_granger['recall'], "F1": m_granger['f1']}
    ])
    summary.to_csv("results/baseline_comparison.csv", index=False)

if __name__ == "__main__":
    run_granger_baseline()
