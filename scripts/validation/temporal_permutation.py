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

def run_permutation_control():
    print("--- Starting Temporal Permutation Control ---")
    
    # 1. Setup Data
    np.random.seed(42)
    n_genes = 15
    n_timepoints = 1000
    true_adj, delays = generate_synthetic_grn(n_genes, connection_prob=0.15, seed=42)
    expression_orig = simulate_expression(n_genes, n_timepoints, true_adj, delays, 
                                          noise_level=0.1, burst_prob=0.02, decay=0.2, seed=42)
    time_points = np.arange(n_timepoints)
    
    # 2. Permute Samples
    perm_indices = np.random.permutation(n_timepoints)
    expression_perm = expression_orig[:, perm_indices]
    
    # 3. Run Pipeline on Both
    def get_grn(expr):
        thresholds = calculate_adaptive_thresholds(expr, sigma=1.5)
        cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
        # We still use the original time_points [0, 1, 2...] but data is shuffled
        # This simulates destroying the temporal relationship between neighbors
        spike_trains = cstdp.compute_spike_times(expr, time_points, thresholds)
        weights = cstdp.run_cstdp(spike_trains, n_genes)
        max_w = np.max(weights)
        norm_w = weights / max_w if max_w > 0 else weights
        return norm_w, weights

    print("Running C-STDP on original data...")
    grn_orig, raw_orig = get_grn(expression_orig)
    print("Running C-STDP on permuted data...")
    grn_perm, raw_perm = get_grn(expression_perm)
    
    # 4. Metrics
    eval_thresh = 0.3
    metrics_orig = calculate_metrics(true_adj, grn_orig, threshold=eval_thresh)
    metrics_perm = calculate_metrics(true_adj, grn_perm, threshold=eval_thresh)
    
    # Jaccard Overlap between Inferred Graphs
    edges_orig = grn_orig > eval_thresh
    edges_perm = grn_perm > eval_thresh
    
    intersection = np.sum(edges_orig & edges_perm)
    union = np.sum(edges_orig | edges_perm)
    jaccard = intersection / union if union > 0 else 0.0
    
    print("\n--- PERMUTATION RESULTS ---")
    print(f"Original Data Precision: {metrics_orig['precision']:.4f}")
    print(f"Permuted Data Precision: {metrics_perm['precision']:.4f}")
    print(f"Overlap (Jaccard):       {jaccard:.4f}")
    
    # 5. Visuals
    os.makedirs("results/visuals", exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sns.heatmap(grn_orig, ax=axes[0], cmap="viridis", vmin=0, vmax=1)
    axes[0].set_title("Original Inferred GRN")
    sns.heatmap(grn_perm, ax=axes[1], cmap="viridis", vmin=0, vmax=1)
    axes[1].set_title("Permuted Inferred GRN")
    plt.savefig("results/visuals/permutation_comparison.png")
    plt.close()
    
    # 6. Weight Distribution Shift
    plt.figure(figsize=(10, 6))
    plt.hist(raw_orig.flatten(), bins=50, alpha=0.5, label="Original", density=True)
    plt.hist(raw_perm.flatten(), bins=50, alpha=0.5, label="Permuted", density=True)
    plt.title("Weight Distribution Shift under Permutation")
    plt.legend()
    plt.savefig("results/visuals/permutation_dist_shift.png")
    plt.close()
    
    # REQUIRED CONCLUSION CHECK
    # Expect precision to collapse to chance or near zero
    chance_level = 0.15
    if metrics_perm['precision'] < metrics_orig['precision'] * 0.7 or metrics_perm['precision'] < chance_level + 0.05:
        print("\n✅ SUCCESS: Temporal structure is required for causality.")
    else:
        print("\n❌ FAILURE: Inferred GRN is invariant to permutation. STDP logic compromised.")
        sys.exit(1)

if __name__ == "__main__":
    run_permutation_control()
