import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
import os

def plot_temporal_phase_progression(input_file, phase_file, output_dir):
    print("Generating Temporal Phase Progression (GSE157859)...")
    
    if not os.path.exists(input_file) or not os.path.exists(phase_file):
        print("Files not found.")
        return
        
    # 1. Load Data
    df = pd.read_csv(input_file, index_col=0)
    data = np.nan_to_num(df.values, nan=0.0)
    time_points = np.arange(data.shape[1])
    
    # 2. Detect Onset in GSE157859
    derivatives = np.gradient(data, axis=1)
    norm_derivatives = np.nan_to_num(zscore(derivatives, axis=1), nan=0.0)
    spike_mask = norm_derivatives > 1.5
    
    onset_times = []
    genes = []
    
    for i in range(len(df)):
        spikes = time_points[spike_mask[i, :]]
        if len(spikes) > 0:
            onset_times.append(float(spikes[0]))
            genes.append(df.index[i])
            
    df_onset = pd.DataFrame({'Gene': genes, 'Onset': onset_times})
    
    # 3. Map to Phases (from GSE215865 definitions)
    phases = pd.read_csv(phase_file)
    # Handle Ensembl/Symbol mismatch?
    # GSE157859 is TPM, indices are Ensembl?
    # Preprocess saved `GSE157859_immune_subset.csv`.
    # Let's check IDs in phase file vs GSE157859.
    # Phase file uses IDs from `activation_ranks.csv` (GSE215865).
    # If both are Ensembl, we are good.
    # Preprocess.py used mapping for both.
    
    # Let's clean IDs to be safe (strip suffix)
    df_onset['CleanID'] = df_onset['Gene'].astype(str).str.split('.').str[0]
    phases['CleanID'] = phases['Gene'].astype(str).str.split('.').str[0]
    
    merged = df_onset.merge(phases[['CleanID', 'Phase']], on='CleanID', how='inner')
    
    print(f"Mapped {len(merged)} genes to phases.")
    
    # 4. Plot
    # Ridgeline or Boxplot
    # X: Time, Y: Phase
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    order = ['Phase I', 'Phase II', 'Phase III', 'Phase IV']
    sns.boxplot(data=merged, x='Onset', y='Phase', order=order, palette='viridis', ax=ax, orient='h')
    
    ax.set_xlabel("Onset Time (Staged Infection)")
    ax.set_ylabel("Immune Response Phase (Defined in Blood Cohort)")
    ax.set_title("Temporal Progression of Immune Phases (GSE157859)")
    
    plt.tight_layout()
    output_path = os.path.join(output_dir, "GSE157859_phase_timeline.png")
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Saved timeline to {output_path}")

if __name__ == "__main__":
    plot_temporal_phase_progression(
        "data/processed/GSE157859_immune_subset.csv",
        "visuals/symposium_final/phase_definitions.csv",
        "visuals/symposium_final"
    )
