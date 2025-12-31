import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.getcwd())

from src.cstdp.stdp import CausalSTDP
from src.utils.spike_encoding import calculate_adaptive_thresholds

def run_negative_controls():
    print("--- Starting Negative Control Genes Analysis ---")
    
    # 1. Load Data
    input_file = "data/processed/GSE215865_subset.csv"
    df = pd.read_csv(input_file, index_col=0)
    real_genes = df.index.values
    n_real = len(real_genes)
    n_samples = df.shape[1]
    
    # 2. Generate Negative Controls (Randomly shuffled real genes)
    np.random.seed(42)
    n_controls = 10
    control_data = []
    control_names = []
    
    for i in range(n_controls):
        # Pick a random real gene and shuffle its temporal order
        gene_to_shuff = np.random.choice(real_genes)
        shuff_profile = np.random.permutation(df.loc[gene_to_shuff].values)
        control_data.append(shuff_profile)
        control_names.append(f"CTRL_{i}_{gene_to_shuff[:10]}")
        
    df_controls = pd.DataFrame(control_data, index=control_names, columns=df.columns)
    
    # Combine
    df_combined = pd.concat([df, df_controls])
    print(f"Combined data: {df_combined.shape[0]} total genes.")
    
    # 3. Infer GRN
    data = df_combined.values
    thresholds = calculate_adaptive_thresholds(data, sigma=1.5)
    cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
    time_points = np.arange(n_samples)
    spike_trains = cstdp.compute_spike_times(data, time_points, thresholds)
    weights, _ = cstdp.run_cstdp(spike_trains, df_combined.shape[0])
    
    # 4. Analyze Out-Degrees
    out_degrees = np.sum(weights, axis=1)
    
    res_real = out_degrees[:n_real]
    res_ctrl = out_degrees[n_real:]
    
    mean_real = np.mean(res_real)
    mean_ctrl = np.mean(res_ctrl)
    
    print("\n--- RESULTS ---")
    print(f"Mean Out-Degree (Real): {mean_real:.4f}")
    print(f"Mean Out-Degree (CTRL): {mean_ctrl:.4f}")
    
    # Check if any CTRL is in Top 10%
    sorted_idx = np.argsort(out_degrees)[::-1]
    top_10_percent_count = int(0.1 * len(out_degrees))
    top_genes = np.array(df_combined.index)[sorted_idx[:top_10_percent_count]]
    
    ctrls_in_top = [g for g in top_genes if g.startswith("CTRL")]
    print(f"Controls in Top 10%: {len(ctrls_in_top)} / {n_controls}")
    
    # 5. Visuals
    os.makedirs("analysis/visuals", exist_ok=True)
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=[res_real, res_ctrl], palette="pastel")
    plt.xticks([0, 1], ["Real Genes", "Negative Controls"])
    plt.ylabel("Out-Degree sum")
    plt.title("Out-Degree Distribution: Real vs Negative Controls")
    plt.savefig("analysis/visuals/negative_control_boxplot.png")
    plt.close()
    
    if len(ctrls_in_top) > n_controls * 0.2: # Allow small noise
        print("\n❌ FAILURE: Negative controls emerged as top regulators.")
        sys.exit(1)
    else:
        print("\n✅ SUCCESS: Negative controls correctly suppressed.")

if __name__ == "__main__":
    run_negative_controls()
