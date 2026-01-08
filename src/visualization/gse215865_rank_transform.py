import pandas as pd
import numpy as np
from scipy.stats import zscore
import os

def calculate_activation_ranks(input_file, output_name):
    print(f"Calculating Activation Ranks for {output_name} (Pseudo-time Cohort)...")
    
    # 1. Load Data
    df = pd.read_csv(input_file, index_col=0)
    data = np.nan_to_num(df.values, nan=0.0)
    genes = df.index.values
    time_points = np.arange(data.shape[1])
    
    # 2. Spike Encoding (Derivative Z-score > 1.5)
    derivatives = np.gradient(data, axis=1)
    norm_derivatives = np.nan_to_num(zscore(derivatives, axis=1), nan=0.0)
    spike_mask = norm_derivatives > 1.5
    
    # 3. Determine Onset Index (First Spike)
    onset_indices = []
    valid_genes = []
    
    for i in range(len(genes)):
        spikes = np.where(spike_mask[i, :])[0]
        if len(spikes) > 0:
            onset_indices.append(spikes[0])
            valid_genes.append(genes[i])
            
    if not onset_indices:
        print("STOP: No onset events detected.")
        return

    onset_indices = np.array(onset_indices)
    
    # Check Variance
    if np.var(onset_indices) < 1e-5:
        print(f"STOP: Onset variance is near zero ({np.var(onset_indices)}). All genes activate simultaneously?")
        # This is a critical check requested by user
        import sys; sys.exit(1)
        
    # 4. Rank Transformation
    # Rank 0 to N-1
    ranks = np.argsort(np.argsort(onset_indices)) # Tied ranks? argsort twice gives rank
    # Better: rankdata
    from scipy.stats import rankdata
    ranks = rankdata(onset_indices, method='min')
    
    # Percentile (0 to 100)
    percentiles = (ranks - 1) / (len(ranks) - 1) * 100
    
    # Summary
    print(f"Total Genes with Onset: {len(valid_genes)}")
    print(f"Onset Index Range: {min(onset_indices)} - {max(onset_indices)}")
    print(f"Rank Range: {min(ranks)} - {max(ranks)}")
    
    # Save Ranks for downstream
    df_ranks = pd.DataFrame({
        'Gene': valid_genes,
        'OnsetIndex': onset_indices,
        'Rank': ranks,
        'Percentile': percentiles
    })
    
    output_path = f"visuals/symposium/{output_name}/activation_ranks.csv"
    df_ranks.to_csv(output_path, index=False)
    print(f"Saved ranks to {output_path}")

if __name__ == "__main__":
    calculate_activation_ranks(
        "data/processed/GSE215865_immune_subset.csv",
        "GSE215865"
    )
