import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.getcwd())

from src.cstdp.core import CausalSTDP
from src.cstdp.utils.spike_encoding import calculate_adaptive_thresholds

def run_cohort_consistency():
    print("--- Starting Cohort Consistency Analysis (GSE215865) ---")
    
    # 1. Load Processed Data
    input_file = "data/processed/GSE215865_subset.csv"
    if not os.path.exists(input_file):
        print("Processed data not found.")
        return
        
    df = pd.read_csv(input_file, index_col=0)
    genes = df.index.values
    n_genes = len(genes)
    n_samples = df.shape[1]
    
    # 2. Split Cohorts
    np.random.seed(42)
    indices = np.random.permutation(n_samples)
    mid = n_samples // 2
    
    cohort1_idx = indices[:mid]
    cohort2_idx = indices[mid:]
    
    df1 = df.iloc[:, cohort1_idx]
    df2 = df.iloc[:, cohort2_idx]
    
    print(f"Cohort 1: {df1.shape[1]} samples, Cohort 2: {df2.shape[1]} samples.")
    
    # 3. Infer GRNs
    def infer(data_df):
        data = data_df.values
        thresholds = calculate_adaptive_thresholds(data, sigma=1.5)
        cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
        # Treated as pseudo-time in each cohort
        time_points = np.arange(data.shape[1])
        spike_trains = cstdp.compute_spike_times(data, time_points, thresholds)
        weights = cstdp.run_cstdp(spike_trains, n_genes)
        
        # Normalize
        mw = np.max(weights)
        norm_w = weights / mw if mw > 0 else weights
        return norm_w

    print("Inferring GRN for Cohort 1...")
    grn1 = infer(df1)
    print("Inferring GRN for Cohort 2...")
    grn2 = infer(df2)
    
    # 4. Measure Overlap
    thresh = 0.3
    edges1 = grn1 > thresh
    edges2 = grn2 > thresh
    
    intersection = np.sum(edges1 & edges2)
    union = np.sum(edges1 | edges2)
    jaccard = intersection / union if union > 0 else 0.0
    
    # Top Regulators Correlation
    out1 = np.sum(grn1, axis=1)
    out2 = np.sum(grn2, axis=1)
    # Correlation between out-degree vectors
    out_corr = np.corrcoef(out1, out2)[0, 1]
    
    print("\n--- CONSISTENCY RESULTS ---")
    print(f"Edge Overlap (Jaccard): {jaccard:.4f}")
    print(f"Regulator Correlation:  {out_corr:.4f}")
    
    # 5. Visuals
    os.makedirs("results/visuals", exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sns.heatmap(grn1, ax=axes[0], cmap="viridis", vmin=0, vmax=1)
    axes[0].set_title("Cohort 1 GRN")
    sns.heatmap(grn2, ax=axes[1], cmap="viridis", vmin=0, vmax=1)
    axes[1].set_title("Cohort 2 GRN")
    plt.savefig("results/visuals/cohort_consistency_heatmaps.png")
    plt.close()
    
    plt.figure(figsize=(8, 8))
    plt.scatter(out1, out2, alpha=0.6)
    plt.xlabel("Cohort 1 Out-Degree")
    plt.ylabel("Cohort 2 Out-Degree")
    plt.title("Stability of Gene Out-Degrees")
    plt.savefig("results/visuals/cohort_regulator_stability.png")
    plt.close()

if __name__ == "__main__":
    run_cohort_consistency()
