import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.getcwd())

from src.cstdp.stdp import CausalSTDP
from src.utils.spike_encoding import calculate_adaptive_thresholds

def run_clinical_constraint():
    print("--- Starting Multi-Modal Constraint (Clinical Proxy) ---")
    
    # 1. Load Data
    input_file = "data/processed/GSE215865_subset.csv"
    df = pd.read_csv(input_file, index_col=0)
    data = df.values
    genes = df.index.values
    n_genes, n_samples = data.shape
    
    # 2. Simulate Clinical Event (Severity Score)
    # Assume samples are time-ordered (pseudo-time). 
    # Severity increases linearly with noise.
    np.random.seed(42)
    severity = np.linspace(0, 1, n_samples) + np.random.normal(0, 0.1, n_samples)
    
    # 3. Identify Severity-Associated Genes (Targets)
    corrs = [np.corrcoef(data[i, :], severity)[0, 1] for i in range(n_genes)]
    targets = [genes[i] for i in range(n_genes) if abs(corrs[i]) > 0.1]
    print(f"identified {len(targets)} severity-associated genes (Targets) with |corr| > 0.1.")
    
    # 4. Infer GRN
    thresholds = calculate_adaptive_thresholds(data, sigma=1.5)
    cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
    spike_trains = cstdp.compute_spike_times(data, np.arange(n_samples), thresholds)
    weights = cstdp.run_cstdp(spike_trains, n_genes)
    
    # Normalize
    mw = np.max(weights)
    norm_w = weights / mw if mw > 0 else weights
    
    # 5. Check Consistency
    # Do Regulators of Targets appear "earlier"?
    # We check if edges (Regulator -> Target) are enriched.
    # Regulator: High Out-Degree
    # Target: High Severity Correlation
    
    regulator_scores = np.sum(norm_w, axis=1)
    target_scores = np.array(corrs)
    
    # We expect Regulators (Sources) to NOT necessarily be correlated with Severity (Late outcome)
    # OR they might be negatively correlated (early peaks).
    # Let's check correlation between "Out-Degree" and "Severity Correlation".
    # If network is Feed-Forward: Sources -> Targets (Severity).
    # Then Sources should have LOW correlation with Severity (active early) 
    # or just distinct from Targets.
    
    consistency_corr = np.corrcoef(regulator_scores, target_scores)[0, 1]
    
    print("\n--- RESULTS ---")
    print(f"Correlation between Out-Degree and Severity-Association: {consistency_corr:.4f}")
    
    # 6. Visuals
    os.makedirs("analysis/visuals", exist_ok=True)
    plt.figure(figsize=(8, 6))
    plt.scatter(target_scores, regulator_scores, alpha=0.7)
    plt.xlabel("Correlation with Clinical Severity (Late)")
    plt.ylabel("Inferred Out-Degree (Early/Source)")
    plt.title("Temporal Ordering Validation")
    plt.axvline(0, color='grey', linestyle='--')
    plt.axhline(0, color='grey', linestyle='--')
    plt.savefig("analysis/visuals/clinical_constraint_check.png")
    plt.close()
    
    print("Interpretation: A low or negative correlation suggests regulators are distinct from late-stage effectors, consistent with a causal cascade.")

if __name__ == "__main__":
    run_clinical_constraint()
