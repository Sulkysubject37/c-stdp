import pandas as pd
import numpy as np
import os
import sys

sys.path.append(os.getcwd())

def select_subgraph_for_viz():
    print("--- Selecting Minimal Subgraph for Visualization ---")
    
    # Paths
    data_path = "data/processed/GSE215865_subset.csv"
    adj_path = "results/visuals/GSE215865/inferred_grn_adj.csv" # Need to check if this exists or where it was saved
    
    # Check if adj exists in results/visuals/GSE215865/ (based on run_real_data_cstdp.py output)
    # The output dir was f"visuals/real_data/{dataset_name}" in the script, which got moved to results/visuals?
    # No, I moved `analysis/visuals` to `results/visuals`.
    # `run_real_data_cstdp.py` wrote to `visuals/real_data/...`.
    # Wait, the refactor moved `pipelines` to `scripts/inference`.
    # Did I move the *output* of the pipelines?
    # The logs said `visuals/real_data/GSE215865` was created.
    # I moved `analysis/visuals` to `results/visuals`.
    # I did NOT move `visuals/real_data`.
    # Let me check where the adj file is.
    
    # I will rely on finding the file.
    
    real_data_path = "visuals/real_data/GSE215865/inferred_grn_adj.csv"
    
    if not os.path.exists(real_data_path):
        print(f"Adjacency file not found at {real_data_path}. Checking alternatives...")
        # Check results/visuals just in case
        if os.path.exists("results/visuals/GSE215865/inferred_grn_adj.csv"):
            real_data_path = "results/visuals/GSE215865/inferred_grn_adj.csv"
        else:
            print("Cannot find adjacency matrix. Run inference first.")
            return

    print(f"Loading adjacency from {real_data_path}...")
    adj_df = pd.read_csv(real_data_path, index_col=0)
    
    # Targets
    # PHC2: ENSG00000134686
    # NSA2: ENSG00000092853
    
    regulators = ["ENSG00000134686", "ENSG00000092853"]
    # Verify they exist in the subset
    valid_regs = [r for r in regulators if r in adj_df.index]
    
    if not valid_regs:
        print("Top regulators PHC2/NSA2 not found in subset. Using top out-degree nodes.")
        out_degree = adj_df.sum(axis=1).sort_values(ascending=False)
        valid_regs = out_degree.head(2).index.tolist()
    
    print(f"Selected Regulators: {valid_regs}")
    
    # Find top targets for these regulators
    subgraph_nodes = set(valid_regs)
    
    for reg in valid_regs:
        # Get row
        targets = adj_df.loc[reg].sort_values(ascending=False)
        # Take top 3 non-zero targets
        top_targets = targets[targets > 0].head(3).index.tolist()
        subgraph_nodes.update(top_targets)
        
    final_nodes = list(subgraph_nodes)
    print(f"Selected {len(final_nodes)} nodes for visualization: {final_nodes}")
    
    # Save selection
    os.makedirs("results/brian2", exist_ok=True)
    with open("results/brian2/selected_genes.txt", "w") as f:
        for gene in final_nodes:
            f.write(gene + "\n")
            
    # Also save the mini-adjacency
    mini_adj = adj_df.loc[final_nodes, final_nodes]
    mini_adj.to_csv("results/brian2/mini_adj.csv")
    print("Saved subgraph to results/brian2/")

if __name__ == "__main__":
    select_subgraph_for_viz()
