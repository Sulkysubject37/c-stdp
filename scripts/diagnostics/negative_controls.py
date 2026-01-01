import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(os.getcwd())

from src.cstdp.core import CausalSTDP
from src.cstdp.utils.spike_encoding import calculate_adaptive_thresholds

def run_negative_controls():
    print("--- Starting Negative Control Genes Analysis ---")
    
    # 1. Load Data
    input_file = "data/processed/GSE215865_subset.csv"
    df = pd.read_csv(input_file, index_col=0)
    real_genes = df.index.values
    n_real = len(real_genes)
    n_samples = df.shape[1]
    
    # 2. Generate Negative Controls
    np.random.seed(42)
    n_controls = 10
    control_data = []
    control_names = []
    
    # Group A: Shuffled Real Profiles
    for i in range(n_controls):
        gene_to_shuff = np.random.choice(real_genes)
        shuff_profile = np.random.permutation(df.loc[gene_to_shuff].values)
        control_data.append(shuff_profile)
        control_names.append(f"SHUFF_{i}")
        
    # Group B: Pure Poisson Noise (approximated by random normal for continuous input)
    # Since we use adaptive thresholds on continuous data, random normal noise 
    # will generate random spikes based on sigma.
    for i in range(n_controls):
        noise_profile = np.random.normal(0, 1, n_samples)
        control_data.append(noise_profile)
        control_names.append(f"NOISE_{i}")
        
    df_controls = pd.DataFrame(control_data, index=control_names, columns=df.columns)
    
    # Combine
    df_combined = pd.concat([df, df_controls])
    print(f"Combined data: {df_combined.shape[0]} total genes (Real + {n_controls} Shuff + {n_controls} Noise).")
    
    # 3. Infer GRN
    data = df_combined.values
    thresholds = calculate_adaptive_thresholds(data, sigma=1.5)
    cstdp = CausalSTDP(w_max=1.0, A_pos=0.05, A_neg=0.06, tau_pos=10, tau_neg=10)
    time_points = np.arange(n_samples)
    spike_trains = cstdp.compute_spike_times(data, time_points, thresholds)
    weights = cstdp.run_cstdp(spike_trains, df_combined.shape[0])
    
    # 4. Analyze Out-Degrees
    out_degrees = np.sum(weights, axis=1)
    
    res_real = out_degrees[:n_real]
    res_shuff = out_degrees[n_real:n_real+n_controls]
    res_noise = out_degrees[n_real+n_controls:]
    
    print("\n--- RESULTS ---")
    print(f"Mean Out-Degree (Real):  {np.mean(res_real):.4f}")
    print(f"Mean Out-Degree (Shuff): {np.mean(res_shuff):.4f}")
    print(f"Mean Out-Degree (Noise): {np.mean(res_noise):.4f}")
    
    # Check Top 10%
    sorted_idx = np.argsort(out_degrees)[::-1]
    top_k = int(0.1 * len(out_degrees))
    top_genes = np.array(df_combined.index)[sorted_idx[:top_k]]
    
    bad_actors = [g for g in top_genes if g.startswith("SHUFF") or g.startswith("NOISE")]
    print(f"Controls in Top {top_k}: {len(bad_actors)} (Names: {bad_actors})")
    
    # 5. Visuals
    os.makedirs("results/visuals", exist_ok=True)
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=[res_real, res_shuff, res_noise], palette="pastel")
    plt.xticks([0, 1, 2], ["Real", "Shuffled", "Poisson Noise"])
    plt.ylabel("Out-Degree sum")
    plt.title("Out-Degree Distribution: Real vs Controls")
    plt.savefig("results/visuals/negative_control_boxplot.png")
    plt.close()
    
    if len(bad_actors) > (2 * n_controls) * 0.2:
        print("\n❌ FAILURE: Controls dominate top regulators.")
    else:
        print("\n✅ SUCCESS: Negative controls suppressed.")

if __name__ == "__main__":
    run_negative_controls()
