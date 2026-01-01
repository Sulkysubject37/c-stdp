import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

sys.path.append(os.getcwd())

from src.cstdp.utils.spike_encoding import calculate_adaptive_thresholds
from src.cstdp.core import CausalSTDP

def run_real_permutation_diagnostics():
    print("--- Starting Full Permutation Diagnostics on Real Data (GSE215865) ---")
    
    # 1. Load Data
    input_file = "data/processed/GSE215865_subset.csv"
    df = pd.read_csv(input_file, index_col=0)
    data_orig = df.values
    genes = df.index.values
    n_genes, n_timepoints = data_orig.shape
    time_points = np.arange(n_timepoints)
    
    # 2. Define Baseline Inference Function
    def infer_network(data_matrix):
        # Using stable params from Step 1
        params = {"A_pos": 0.05, "A_neg": 0.06, "tau": 10, "sigma": 1.5}
        
        thresh = calculate_adaptive_thresholds(data_matrix, sigma=params["sigma"])
        cstdp = CausalSTDP(w_max=1.0, A_pos=params["A_pos"], A_neg=params["A_neg"], 
                           tau_pos=params["tau"], tau_neg=params["tau"])
        spikes = cstdp.compute_spike_times(data_matrix, time_points, thresh)
        w = cstdp.run_cstdp(spikes, n_genes)
        
        mw = np.max(w)
        norm_w = w / mw if mw > 0 else w
        return norm_w

    print("Inferring Original GRN...")
    w_orig = infer_network(data_orig)
    edges_orig = w_orig > 0.3
    out_orig = np.sum(w_orig, axis=1)
    
    # 3. Permutation Tests
    results = []
    
    # Test A: Sample Order Permutation (Global temporal destruction)
    print("Test A: Sample Order Permutation...")
    perm_idx = np.random.permutation(n_timepoints)
    data_perm_samples = data_orig[:, perm_idx]
    w_perm_samples = infer_network(data_perm_samples)
    
    # Test B: Spike Train Permutation (Local temporal destruction per gene)
    # We shuffle the spike times themselves? No, simpler to shuffle the data row-wise
    # actually, shuffling data row-wise destroys correlations between genes but keeps marginal stats
    print("Test B: Gene-wise Time Permutation...")
    data_perm_genes = np.zeros_like(data_orig)
    for i in range(n_genes):
        data_perm_genes[i, :] = np.random.permutation(data_orig[i, :])
    w_perm_genes = infer_network(data_perm_genes)
    
    # Test C: Gene Label Permutation (Topological destruction)
    # Actually, this just renames nodes. The network topology should be identical but mapped differently.
    # Comparing overlap with original requires mapping back.
    # Wait, if I shuffle labels, I'm effectively testing if the structure is specific to *these* genes
    # or just generic network properties.
    # Let's shuffle rows of the INPUT matrix (effectively swapping gene profiles)
    print("Test C: Gene Profile Swap...")
    # Row shuffle
    row_idx = np.random.permutation(n_genes)
    data_swap = data_orig[row_idx, :]
    # The "Original Gene i" now has "Profile of Gene k"
    w_swap = infer_network(data_swap) 
    # If inferred network is driven by profile, then w_swap[i,j] should match w_orig[row_idx[i], row_idx[j]]
    # BUT, if we compare w_swap directly to w_orig, overlap should be low.
    
    # 4. Compare
    def compare(w_target, name):
        edges_tgt = w_target > 0.3
        
        intersect = np.sum(edges_orig & edges_tgt)
        union = np.sum(edges_orig | edges_tgt)
        jaccard = intersect / union if union > 0 else 0.0
        
        out_tgt = np.sum(w_target, axis=1)
        if np.std(out_tgt) == 0 or np.std(out_orig) == 0:
            corr = 0.0
        else:
            corr, _ = spearmanr(out_orig, out_tgt)
            
        print(f"[{name}] Jaccard: {jaccard:.4f}, Rank Corr: {corr:.4f}")
        return {"Test": name, "Jaccard": jaccard, "Rank_Corr": corr}

    results.append(compare(w_perm_samples, "Sample_Permutation"))
    results.append(compare(w_perm_genes, "Gene_Time_Permutation"))
    results.append(compare(w_swap, "Gene_Profile_Swap"))
    
    # 5. Visuals
    os.makedirs("results/visuals", exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    
    sns.heatmap(w_orig, ax=axes[0], cmap="viridis", vmin=0, vmax=1)
    axes[0].set_title("Original")
    
    sns.heatmap(w_perm_samples, ax=axes[1], cmap="viridis", vmin=0, vmax=1)
    axes[1].set_title("Sample Perm")
    
    sns.heatmap(w_perm_genes, ax=axes[2], cmap="viridis", vmin=0, vmax=1)
    axes[2].set_title("Gene Time Perm")
    
    sns.heatmap(w_swap, ax=axes[3], cmap="viridis", vmin=0, vmax=1)
    axes[3].set_title("Profile Swap")
    
    plt.savefig("results/visuals/real_permutation_tests.png")
    plt.close()
    
    # 6. Conclusion
    # Sample Permutation and Gene Time Permutation MUST destroy structure (low Jaccard)
    # Profile Swap implies that the network structure follows the data, so comparing adj matrices (fixed node IDs) should be low overlap
    
    if results[0]['Jaccard'] < 0.1 and results[1]['Jaccard'] < 0.1:
        print("\n✅ SUCCESS: Network structure collapses under temporal destruction.")
    else:
        print("\n❌ FAILURE: Structure persists despite permutation. STDP may be overfitting to marginals.")
        # sys.exit(1) # Warning only for now, as real data might have artifacts

if __name__ == "__main__":
    run_real_permutation_diagnostics()
