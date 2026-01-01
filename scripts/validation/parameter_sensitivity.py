import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product

sys.path.append(os.getcwd())

from src.cstdp.utils.simulate_grn import generate_synthetic_grn, simulate_expression
from src.cstdp.utils.spike_encoding import calculate_adaptive_thresholds
from src.cstdp.core import CausalSTDP
from src.cstdp.utils.evaluate import calculate_metrics

def run_sensitivity_analysis():
    print("--- Starting Parameter Sensitivity Analysis ---")
    
    # 1. Fixed Synthetic Data (for fairness)
    np.random.seed(42)
    n_genes = 15
    n_timepoints = 1000
    true_adj, delays = generate_synthetic_grn(n_genes, connection_prob=0.15, seed=42)
    expression = simulate_expression(n_genes, n_timepoints, true_adj, delays, 
                                     noise_level=0.1, burst_prob=0.02, decay=0.2, seed=42)
    time_points = np.arange(n_timepoints)
    
    # 2. Parameter Grid
    # Focused sweep based on previous knowledge (A_neg > A_pos is critical)
    
    a_pos_vals = [0.01, 0.05, 0.10]
    a_neg_multipliers = [0.8, 1.0, 1.1, 1.2, 1.5, 2.0] # Relative to A_pos
    tau_vals = [5, 10, 20, 50]
    sigma_vals = [1.0, 1.5, 2.0, 3.0]
    
    results = []
    
    # We will do a few slices to keep compute reasonable
    # Slice 1: A_pos vs A_neg (fixed tau=10, sigma=1.5)
    print("Running Slice 1: Amplitude Sensitivity...")
    for a_p, mul in product(a_pos_vals, a_neg_multipliers):
        a_n = a_p * mul
        
        # Run Pipeline
        thresholds = calculate_adaptive_thresholds(expression, sigma=1.5)
        cstdp = CausalSTDP(w_max=1.0, A_pos=a_p, A_neg=a_n, tau_pos=10, tau_neg=10)
        spike_trains = cstdp.compute_spike_times(expression, time_points, thresholds)
        weights = cstdp.run_cstdp(spike_trains, n_genes)
        
        # Normalize
        max_w = np.max(weights)
        norm_w = weights / max_w if max_w > 0 else weights
        
        # Evaluate
        metrics = calculate_metrics(true_adj, norm_w, threshold=0.3)
        
        # Sparsity
        sparsity = np.mean(norm_w == 0)
        
        results.append({
            "Experiment": "Amplitude",
            "A_pos": a_p,
            "A_neg": a_n,
            "Ratio_An_Ap": mul,
            "Tau": 10,
            "Sigma": 1.5,
            "Precision": metrics['precision'],
            "Recall": metrics['recall'],
            "F1": metrics['f1'],
            "SHD": metrics['shd'],
            "Sparsity": sparsity
        })

    # Slice 2: Tau vs Sigma (fixed A_pos=0.05, A_neg=0.06)
    print("Running Slice 2: Time/Threshold Sensitivity...")
    for tau, sig in product(tau_vals, sigma_vals):
        # Run Pipeline
        thresholds = calculate_adaptive_thresholds(expression, sigma=sig)
        cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=tau, tau_neg=tau)
        spike_trains = cstdp.compute_spike_times(expression, time_points, thresholds)
        weights = cstdp.run_cstdp(spike_trains, n_genes)
        
        # Normalize
        max_w = np.max(weights)
        norm_w = weights / max_w if max_w > 0 else weights
        
        metrics = calculate_metrics(true_adj, norm_w, threshold=0.3)
        sparsity = np.mean(norm_w == 0)
        
        results.append({
            "Experiment": "TauSigma",
            "A_pos": 0.05,
            "A_neg": 0.06,
            "Ratio_An_Ap": 1.2,
            "Tau": tau,
            "Sigma": sig,
            "Precision": metrics['precision'],
            "Recall": metrics['recall'],
            "F1": metrics['f1'],
            "SHD": metrics['shd'],
            "Sparsity": sparsity
        })
        
    df = pd.DataFrame(results)
    df.to_csv("analysis/parameter_sensitivity.csv", index=False)
    print("Saved results to analysis/parameter_sensitivity.csv")
    
    # 3. Generate Heatmaps
    generate_heatmaps(df)

def generate_heatmaps(df):
    os.makedirs("analysis/visuals", exist_ok=True)
    
    # Heatmap 1: A_pos vs Ratio (Amplitude Stability)
    subset1 = df[df["Experiment"] == "Amplitude"]
    pivot1 = subset1.pivot(index="A_pos", columns="Ratio_An_Ap", values="F1")
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(pivot1, annot=True, cmap="viridis", vmin=0, vmax=1)
    plt.title("F1 Score Stability: Amplitude (Tau=10, Sigma=1.5)")
    plt.ylabel("A_pos")
    plt.xlabel("Ratio (A_neg / A_pos)")
    plt.savefig("analysis/visuals/heatmap_amplitude_f1.png")
    plt.close()
    
    # Heatmap 2: Tau vs Sigma (Time Stability)
    subset2 = df[df["Experiment"] == "TauSigma"]
    pivot2 = subset2.pivot(index="Tau", columns="Sigma", values="F1")
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(pivot2, annot=True, cmap="viridis", vmin=0, vmax=1)
    plt.title("F1 Score Stability: Time Constants (A_pos=0.05, A_neg=0.06)")
    plt.ylabel("Tau")
    plt.xlabel("Sigma (Threshold)")
    plt.savefig("analysis/visuals/heatmap_tau_sigma_f1.png")
    plt.close()
    
    print("Saved heatmaps to analysis/visuals/")

if __name__ == "__main__":
    run_sensitivity_analysis()
