import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product
from scipy.stats import spearmanr

sys.path.append(os.getcwd())

from src.cstdp.utils.spike_encoding import calculate_adaptive_thresholds
from src.cstdp.core import CausalSTDP

def run_real_sensitivity_analysis():
    print("--- Starting Parameter Sensitivity Analysis on Real Data (GSE215865) ---")
    
    # 1. Load Real Data
    input_file = "data/processed/GSE215865_subset.csv"
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return
        
    df = pd.read_csv(input_file, index_col=0)
    data = df.values
    genes = df.index.values
    n_genes, n_timepoints = data.shape
    time_points = np.arange(n_timepoints)
    
    # 2. Define Baseline (Reference) Parameters
    # From synthetic analysis: A_pos=0.05, A_neg=0.06 (Ratio 1.2), Tau=10, Sigma=1.5
    base_params = {
        "A_pos": 0.05,
        "A_neg": 0.06,
        "tau": 10,
        "sigma": 1.5
    }
    
    # Get Baseline GRN for comparison
    print("Computing Baseline GRN...")
    thresh_base = calculate_adaptive_thresholds(data, sigma=base_params["sigma"])
    cstdp_base = CausalSTDP(w_max=1.0, A_pos=base_params["A_pos"], A_neg=base_params["A_neg"], 
                            tau_pos=base_params["tau"], tau_neg=base_params["tau"])
    spikes_base = cstdp_base.compute_spike_times(data, time_points, thresh_base)
    w_base = cstdp_base.run_cstdp(spikes_base, n_genes)
    
    # Normalize Baseline
    mw_b = np.max(w_base)
    norm_w_base = w_base / mw_b if mw_b > 0 else w_base
    edges_base = norm_w_base > 0.3 # Threshold from synthetic
    out_degree_base = np.sum(norm_w_base, axis=1)
    
    # 3. Parameter Grid (Same as synthetic)
    a_pos_vals = [0.01, 0.05, 0.10]
    a_neg_multipliers = [0.8, 1.0, 1.1, 1.2, 1.5, 2.0]
    tau_vals = [5, 10, 20, 50]
    sigma_vals = [1.0, 1.5, 2.0, 3.0]
    
    results = []
    
    # Slice 1: Amplitude Sensitivity
    print("Running Slice 1: Amplitude Sensitivity...")
    for a_p, mul in product(a_pos_vals, a_neg_multipliers):
        a_n = a_p * mul
        
        # Run
        thresh = calculate_adaptive_thresholds(data, sigma=base_params["sigma"])
        cstdp = CausalSTDP(w_max=1.0, A_pos=a_p, A_neg=a_n, 
                           tau_pos=base_params["tau"], tau_neg=base_params["tau"])
        spikes = cstdp.compute_spike_times(data, time_points, thresh)
        w = cstdp.run_cstdp(spikes, n_genes)
        
        # Metrics
        mw = np.max(w)
        norm_w = w / mw if mw > 0 else w
        edges = norm_w > 0.3
        
        # Sparsity
        sparsity = np.mean(edges == 0)
        
        # Jaccard vs Baseline
        intersect = np.sum(edges & edges_base)
        union = np.sum(edges | edges_base)
        jaccard = intersect / union if union > 0 else 0.0
        
        # Rank Correlation of Out-Degrees
        out_degree = np.sum(norm_w, axis=1)
        # Handle constant output (all zeros)
        if np.std(out_degree) == 0 or np.std(out_degree_base) == 0:
            corr = 0.0
        else:
            corr, _ = spearmanr(out_degree, out_degree_base)
            
        results.append({
            "Experiment": "Amplitude",
            "A_pos": a_p,
            "Ratio": mul,
            "Tau": base_params["tau"],
            "Sigma": base_params["sigma"],
            "Sparsity": sparsity,
            "Jaccard_vs_Base": jaccard,
            "Degree_Rank_Corr": corr
        })

    # Slice 2: Time/Threshold Sensitivity
    print("Running Slice 2: Time/Threshold Sensitivity...")
    for tau, sig in product(tau_vals, sigma_vals):
        # Run
        thresh = calculate_adaptive_thresholds(data, sigma=sig)
        cstdp = CausalSTDP(w_max=1.0, A_pos=base_params["A_pos"], A_neg=base_params["A_neg"], 
                           tau_pos=tau, tau_neg=tau)
        spikes = cstdp.compute_spike_times(data, time_points, thresh)
        w = cstdp.run_cstdp(spikes, n_genes)
        
        # Metrics
        mw = np.max(w)
        norm_w = w / mw if mw > 0 else w
        edges = norm_w > 0.3
        
        sparsity = np.mean(edges == 0)
        
        intersect = np.sum(edges & edges_base)
        union = np.sum(edges | edges_base)
        jaccard = intersect / union if union > 0 else 0.0
        
        out_degree = np.sum(norm_w, axis=1)
        if np.std(out_degree) == 0 or np.std(out_degree_base) == 0:
            corr = 0.0
        else:
            corr, _ = spearmanr(out_degree, out_degree_base)
            
        results.append({
            "Experiment": "TauSigma",
            "A_pos": base_params["A_pos"],
            "Ratio": 1.2,
            "Tau": tau,
            "Sigma": sig,
            "Sparsity": sparsity,
            "Jaccard_vs_Base": jaccard,
            "Degree_Rank_Corr": corr
        })
        
    df_res = pd.DataFrame(results)
    df_res.to_csv("results/real_data_sensitivity.csv", index=False)
    
    # Heatmaps
    plot_heatmaps(df_res)

def plot_heatmaps(df):
    os.makedirs("results/visuals", exist_ok=True)
    
    # 1. Amplitude Stability (Degree Correlation)
    sub1 = df[df["Experiment"] == "Amplitude"]
    piv1 = sub1.pivot(index="A_pos", columns="Ratio", values="Degree_Rank_Corr")
    
    plt.figure(figsize=(8,6))
    sns.heatmap(piv1, annot=True, cmap="coolwarm", vmin=0, vmax=1)
    plt.title("Real Data Stability: Regulator Rank Correlation\n(Amplitude Sweep)")
    plt.ylabel("A_pos")
    plt.xlabel("Ratio (A_neg/A_pos)")
    plt.savefig("results/visuals/real_heatmap_amp_corr.png")
    plt.close()
    
    # 2. Tau/Sigma Stability (Degree Correlation)
    sub2 = df[df["Experiment"] == "TauSigma"]
    piv2 = sub2.pivot(index="Tau", columns="Sigma", values="Degree_Rank_Corr")
    
    plt.figure(figsize=(8,6))
    sns.heatmap(piv2, annot=True, cmap="coolwarm", vmin=0, vmax=1)
    plt.title("Real Data Stability: Regulator Rank Correlation\n(Time/Threshold Sweep)")
    plt.ylabel("Tau")
    plt.xlabel("Sigma")
    plt.savefig("results/visuals/real_heatmap_tausigma_corr.png")
    plt.close()
    
    print("Visuals saved to results/visuals/")

if __name__ == "__main__":
    run_real_sensitivity_analysis()
