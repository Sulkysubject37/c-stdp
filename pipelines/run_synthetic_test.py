import sys
import os
sys.path.append(os.getcwd())

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from src.utils.simulate_grn import generate_synthetic_grn, simulate_expression
from src.utils.spike_encoding import calculate_adaptive_thresholds, plot_raster
from src.cstdp.stdp import CausalSTDP
from src.utils.evaluate import calculate_metrics

def main():
    # 1. Setup
    np.random.seed(42)
    n_genes = 15
    n_timepoints = 1000
    
    print("1. Generating Synthetic GRN...")
    true_adj, delays = generate_synthetic_grn(n_genes, connection_prob=0.15)
    
    print("2. Simulating Expression...")
    expression = simulate_expression(n_genes, n_timepoints, true_adj, delays, 
                                     noise_level=0.1, burst_prob=0.02, decay=0.2)
    
    # 3. Spike Encoding
    print("3. Encoding Spikes...")
    thresholds = calculate_adaptive_thresholds(expression, sigma=1.5)
    
    # Use the compute_spike_times from STDP class (or wrapper)
    # We instantiate STDP first
    # FIX: Increased A_neg to 0.06 (> A_pos) to favor depression and sparsity.
    cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
    
    time_points = np.arange(n_timepoints)
    spike_trains = cstdp.compute_spike_times(expression, time_points, thresholds)
    
    # Plot Raster
    print("   Saving Raster Plot...")
    plot_raster(spike_trains, (0, 200))
    plt.savefig("visuals/synthetic/raster_plot.png")
    plt.close()
    
    # 4. Run C-STDP
    print("4. Running C-STDP...")
    inferred_weights = cstdp.run_cstdp(spike_trains, n_genes)
    
    # Debug: Weight Stats
    print("\n[Weight Stats]")
    print(f"Mean: {np.mean(inferred_weights):.4f}, Std: {np.std(inferred_weights):.4f}")
    print(f"Max: {np.max(inferred_weights):.4f}")
    print(f"25%: {np.percentile(inferred_weights, 25):.4f}, 50%: {np.percentile(inferred_weights, 50):.4f}, 75%: {np.percentile(inferred_weights, 75):.4f}")
    
    # 5. Evaluate
    # Normalize inferred weights for comparison
    max_w = np.max(inferred_weights)
    if max_w > 0:
        inferred_norm = inferred_weights / max_w
    else:
        inferred_norm = inferred_weights
        
    # Threshold for evaluation
    eval_threshold = 0.3
    metrics = calculate_metrics(true_adj, inferred_norm, threshold=eval_threshold)
    
    print("\n--- RESULTS ---")
    print(f"Precision: {metrics['precision']:.4f}")
    print(f"Recall:    {metrics['recall']:.4f}")
    print(f"F1 Score:  {metrics['f1']:.4f}")
    print(f"SHD:       {metrics['shd']}")
    
    # 6. Visualize Matrices
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    sns.heatmap(true_adj, ax=axes[0], cmap="viridis", vmin=0, vmax=1)
    axes[0].set_title("Ground Truth")
    
    sns.heatmap(inferred_norm, ax=axes[1], cmap="viridis", vmin=0, vmax=1)
    axes[1].set_title("Inferred Weights (Normalized)")
    
    # Binarized
    sns.heatmap((inferred_norm > eval_threshold).astype(int), ax=axes[2], cmap="Greys", cbar=False)
    axes[2].set_title(f"Inferred Binary (Thresh={eval_threshold})")
    
    plt.savefig("visuals/synthetic/comparison_matrices.png")
    plt.close()
    
    # Check Success Criteria
    # "Recovered above chance"
    # Chance precision for sparse graph (prob 0.15) is ~0.15
    if metrics['precision'] > 0.2 and metrics['recall'] > 0.1:
        print("\n✅ SUCCESS: Synthetic GRN recovered above chance.")
    else:
        print("\n❌ FAILURE: Recovery poor.")
        sys.exit(1)

if __name__ == "__main__":
    main()
