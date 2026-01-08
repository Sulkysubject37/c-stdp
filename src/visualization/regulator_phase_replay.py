import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import numpy as np
import os
from scipy.stats import zscore

def generate_spotlight(regulator_base, dataset_name, adj_file, input_file, phase_file, output_path):
    print(f"Generating Spotlight for {regulator_base} in {dataset_name}...")
    
    # 1. Load Data
    adj = pd.read_csv(adj_file, index_col=0)
    phases = pd.read_csv(phase_file)
    df_expr = pd.read_csv(input_file, index_col=0)
    
    # Find full regulator ID
    regulator_id = None
    for g in adj.index:
        if str(g).startswith(regulator_base):
            regulator_id = g
            break
    
    if not regulator_id:
        print(f"  Regulator {regulator_base} not found.")
        return

    # Find Targets and their Phases
    targets = adj.loc[regulator_id]
    targets = targets[targets > 0.01].index.tolist()
    
    gene_to_phase = pd.Series(phases.Phase.values, index=phases.Gene).to_dict()
    
    # 2. Determine Activation Times
    # For simplicity, we use the Rank/Percentile as the time axis [0, 100]
    gene_to_pct = pd.Series(phases.Percentile.values, index=phases.Gene).to_dict()
    
    reg_time = gene_to_pct.get(regulator_id.split('.')[0], 0)
    
    # Activation times for phases (earliest target in that phase)
    phase_activation = {
        'Phase I': 999, 'Phase II': 999, 'Phase III': 999, 'Phase IV': 999
    }
    
    # The regulator itself activates its own phase
    reg_phase = gene_to_phase.get(regulator_id.split('.')[0], 'Phase I')
    phase_activation[reg_phase] = reg_time
    
    for t in targets:
        clean_t = str(t).split('.')[0]
        p = gene_to_phase.get(clean_t)
        pct = gene_to_pct.get(clean_t)
        if p in phase_activation:
            if pct < phase_activation[p]:
                phase_activation[p] = pct
                
    print(f"  Phase Activation Times: {phase_activation}")

    # 3. Animation
    fig, ax = plt.subplots(figsize=(8, 6))
    
    phases_list = ['Phase I', 'Phase II', 'Phase III', 'Phase IV']
    pos = {p: (i, 0) for i, p in enumerate(phases_list)}
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    def update(frame):
        ax.clear()
        current_time = frame
        
        # Draw 4 large circles for phases
        for i, p in enumerate(phases_list):
            is_active = phase_activation[p] <= current_time
            color = colors[i] if is_active else 'lightgray'
            alpha = 1.0 if is_active else 0.3
            
            circle = plt.Circle((i, 0), 0.4, color=color, alpha=alpha, ec='black')
            ax.add_patch(circle)
            ax.text(i, 0, p, ha='center', va='center', fontweight='bold', color='white' if is_active else 'black')
            
            # Special label for regulator
            if p == reg_phase and is_active:
                ax.text(i, 0.5, f"[{regulator_base} ACTIVES]", ha='center', color='gold', fontweight='bold')

        ax.set_xlim(-1, 4)
        ax.set_ylim(-1, 1)
        ax.set_title(f"Regulator Spotlight: {regulator_base} Influence ({dataset_name})\nT = {current_time:.1f}% Activation Order")
        ax.axis('off')

    ani = animation.FuncAnimation(fig, update, frames=np.linspace(0, 105, 50), interval=100)
    
    ani.save(output_path, writer='ffmpeg', fps=10)
    plt.close()
    print(f"  Saved spotlight to {output_path}")

if __name__ == "__main__":
    PHASE_FILE = "visuals/symposium_final/phase_definitions.csv"
    
    # HEXIM1
    generate_spotlight(
        "ENSG00000186834", "GSE215865",
        "visuals/real_data/GSE215865_Immune/inferred_grn_adj.csv",
        "data/processed/GSE215865_immune_subset.csv",
        PHASE_FILE,
        "visuals/symposium_final/hexim1_spotlight.mp4"
    )
    
    # HMOX1
    generate_spotlight(
        "ENSG00000100292", "GSE157859",
        "visuals/real_data/GSE157859_Immune/inferred_grn_adj.csv",
        "data/processed/GSE157859_immune_subset.csv",
        PHASE_FILE,
        "visuals/symposium_final/hmox1_spotlight.mp4"
    )
