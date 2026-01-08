import pandas as pd
import numpy as np
import os
import sys

def compute_influence(adj_file, phase_file, output_dir):
    print("Computing Phase-to-Phase Influence...")
    
    if not os.path.exists(adj_file) or not os.path.exists(phase_file):
        print("Files not found.")
        return
        
    # Load Adjacency
    adj = pd.read_csv(adj_file, index_col=0)
    
    # Load Phases
    phases = pd.read_csv(phase_file)
    # Map Gene -> Phase
    # Note: Adj genes might have version suffixes. Phase file might not?
    # Let's check.
    # rank_file was used for phase definition.
    # In 'gse215865_rank_transform.py', valid_genes came from df.index.
    # df.index came from `data/processed/GSE215865_immune_subset.csv`.
    # `preprocess.py` saved Ensembl IDs (with suffixes stripped in logic, but let's verify).
    # Preprocess stripped suffixes: "df.index = [x.split('.')[0] for x in df.index]"
    # So `phase_definitions.csv` has clean IDs.
    # Adjacency matrix from `run_real_data_cstdp.py` (vectorized) used `genes_sub` from `df.index` of the same file.
    # So they should match.
    
    gene_to_phase = pd.Series(phases.Phase.values, index=phases.Gene).to_dict()
    
    # Aggregate
    # Create Matrix 4x4
    phase_order = ['Phase I', 'Phase II', 'Phase III', 'Phase IV']
    inf_matrix = pd.DataFrame(0.0, index=phase_order, columns=phase_order)
    
    # Iterate Adjacency (only non-zero)
    # This might be slow if dense. Vectorize?
    # Adj is 200x200 (subset). Fast.
    # Wait, 'inferred_grn_adj.csv' is 200x200?
    # Yes, from `run_real_data_cstdp.py`.
    # But `phase_definitions.csv` has 3488 genes.
    # We only care about the 200 genes in the Adjacency matrix.
    
    for u in adj.index:
        for v in adj.columns:
            w = adj.loc[u, v]
            if w > 0:
                p_u = gene_to_phase.get(u, 'Unknown')
                p_v = gene_to_phase.get(v, 'Unknown')
                
                if p_u in phase_order and p_v in phase_order:
                    inf_matrix.loc[p_u, p_v] += w
                    
    # Normalize?
    # "Sum weights per phase pair. Normalize for visualization only."
    # I'll save raw weights.
    
    print("\nPhase Influence Matrix (Raw Sum of Weights):")
    print(inf_matrix)
    
    # Dominant Flows
    print("\nDominant Flows (> 1.0 total weight):")
    for u in phase_order:
        for v in phase_order:
            if inf_matrix.loc[u, v] > 1.0:
                print(f"  {u} -> {v}: {inf_matrix.loc[u, v]:.2f}")
                
    output_path = os.path.join(output_dir, "phase_influence.csv")
    inf_matrix.to_csv(output_path)
    print(f"Saved influence matrix to {output_path}")

if __name__ == "__main__":
    compute_influence(
        "visuals/real_data/GSE215865_Immune/inferred_grn_adj.csv",
        "visuals/symposium_final/phase_definitions.csv",
        "visuals/symposium_final"
    )
