import pandas as pd
import os
import sys

def preprocess_immune_subset(input_file, mapping_file, immune_file, output_file):
    print("Loading data...")
    # Read matrix
    df = pd.read_csv(input_file, compression='gzip', index_col=0)
    print(f"Original shape: {df.shape}")
    
    # Strip versions from index
    df.index = [x.split('.')[0] for x in df.index]
    
    # Load mapping
    print("Loading mapping...")
    map_df = pd.read_csv(mapping_file)
    # Create dict: Ensembl -> Symbol
    ens_to_sym = pd.Series(map_df.symbol.values, index=map_df.ensembl_id).to_dict()
    
    # Load immune universe
    print("Loading immune universe...")
    with open(immune_file, 'r') as f:
        immune_genes = set([line.strip() for line in f])
        
    # Map index to symbols
    # We want to keep rows where the mapped symbol is in immune_genes
    
    keep_indices = []
    mapped_symbols = []
    
    for ens_id in df.index:
        sym = ens_to_sym.get(ens_id)
        if sym and sym in immune_genes:
            keep_indices.append(True)
            mapped_symbols.append(sym)
        else:
            keep_indices.append(False)
            
    df_immune = df[keep_indices]
    
    print(f"Filtered shape: {df_immune.shape}")
    
    retained_count = df_immune.shape[0]
    percent = (retained_count / df.shape[0]) * 100
    print(f"Percentage retained: {percent:.2f}%")
    
    if retained_count < 50:
        print(f"STOP: Retained genes ({retained_count}) < 50.")
        sys.exit(1)
    if retained_count > 5000:
        print(f"STOP: Retained genes ({retained_count}) > 5000.")
        # The prompt says > 5000 is a stop condition.
        # My universe is 4352, so it shouldn't exceed 4352 unless I have dups or logic error.
        # But wait, if my universe was larger, this would be a blocker.
        # However, 4352 is < 5000, so I should be fine.
        # Just in case:
        if retained_count > 5000:
             sys.exit(1)

    # Save
    df_immune.to_csv(output_file)
    print(f"Saved to {output_file}")

if __name__ == "__main__":
    INPUT_FILE = "data/raw/GSE215865/GSE215865_rnaseq_logCPM_matrix.csv.gz"
    MAPPING_FILE = "data/processed/gene_id_mapping.csv"
    IMMUNE_FILE = "data/processed/immune_universe.txt"
    OUTPUT_FILE = "data/processed/GSE215865_immune_subset.csv"
    
    preprocess_immune_subset(INPUT_FILE, MAPPING_FILE, IMMUNE_FILE, OUTPUT_FILE)
