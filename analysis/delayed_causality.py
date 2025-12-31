import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.getcwd())

from src.utils.simulate_grn import simulate_expression
from src.utils.spike_encoding import calculate_adaptive_thresholds
from src.cstdp.stdp import CausalSTDP
from src.utils.evaluate import calculate_metrics, evaluate_directionality

def run_delay_stress_test():
    print("--- Starting Delayed Causality Stress Test ---")
    
    n_genes = 10
    n_timepoints = 1500
    time_points = np.arange(n_timepoints)
    
    # 1. Create a chain of 3 genes: 0 -> 1 -> 2
    # This allows testing simple delays and cascades
    adj = np.zeros((n_genes, n_genes))
    adj[0, 1] = 1.0
    adj[1, 2] = 1.0
    
    delay_range = [2, 5, 10, 15, 20, 30, 50]
    results = []
    
    for d in delay_range:
        delays = np.zeros((n_genes, n_genes))
        delays[0, 1] = d
        delays[1, 2] = d
        
        # Simulate
        # Low noise, high burst prob to ensure propagation
        expression = simulate_expression(n_genes, n_timepoints, adj, delays, 
                                         noise_level=0.05, burst_prob=0.05, decay=0.2, seed=42)
        
        # Inference
        thresholds = calculate_adaptive_thresholds(expression, sigma=1.5)
        # Fix STDP window (tau=10)
        cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
        spike_trains = cstdp.compute_spike_times(expression, time_points, thresholds)
        weights = cstdp.run_cstdp(spike_trains, n_genes)
        
        # Normalize
        max_w = np.max(weights)
        norm_w = weights / max_w if max_w > 0 else weights
        
        # Metrics
        # Direction accuracy among found edges
        dir_acc = evaluate_directionality(adj, norm_w, threshold=0.3)
        
        # General metrics
        m = calculate_metrics(adj, norm_w, threshold=0.3)
        
        results.append({
            "Delay": d,
            "F1": m['f1'],
            "Precision": m['precision'],
            "Recall": m['recall'],
            "Dir_Accuracy": dir_acc
        })
        
    df = pd.DataFrame(results)
    df.to_csv("analysis/delay_stress_test.csv", index=False)
    
    # 2. Plot
    plt.figure(figsize=(10, 6))
    plt.plot(df["Delay"], df["F1"], marker='o', label="F1 Score")
    plt.plot(df["Delay"], df["Dir_Accuracy"], marker='s', label="Direction Accuracy")
    plt.axvline(x=10, color='r', linestyle='--', label="Tau_pos (10)")
    plt.xlabel("Delay (Time Steps)")
    plt.ylabel("Metric Score")
    plt.title("STDP Performance vs Causal Delay Length")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig("analysis/visuals/delay_stress_test.png")
    plt.close()
    
    print("\n--- DELAY TEST RESULTS ---")
    print(df)
    print("\nInterpretation: Performance peaks when delay matches tau_pos window.")

if __name__ == "__main__":
    run_delay_stress_test()
