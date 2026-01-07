import pandas as pd
import mygene
import os
import sys

def harmonize_gene_ids(input_file, output_map_file):
    """
    Reads Ensembl IDs from input matrix, strips versions, maps to Symbols using mygene.
    """
    print(f"Loading IDs from {input_file}...")
    try:
        # Just read the index
        df = pd.read_csv(input_file, compression='gzip', index_col=0, usecols=[0])
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    ensembl_ids_with_version = df.index.tolist()
    # Strip versions
    ensembl_ids_clean = [x.split('.')[0] for x in ensembl_ids_with_version]
    
    unique_ids = list(set(ensembl_ids_clean))
    print(f"Found {len(unique_ids)} unique Ensembl IDs.")

    # Check cache
    if os.path.exists(output_map_file):
        print(f"Loading mapping from cache: {output_map_file}")
        mapping_df = pd.read_csv(output_map_file)
    else:
        print("Querying mygene...")
        mg = mygene.MyGeneInfo()
        # Query in batches if necessary, but mygene handles list queries well
        # We want 'symbol'
        results = mg.querymany(unique_ids, scopes='ensembl.gene', fields='symbol', species='human', returnall=False)
        
        # Convert to dataframe
        res_data = []
        for res in results:
            if 'symbol' in res:
                res_data.append({'ensembl_id': res['query'], 'symbol': res['symbol']})
            else:
                res_data.append({'ensembl_id': res['query'], 'symbol': None}) # Keep track of unmapped
        
        mapping_df = pd.DataFrame(res_data)
        # Save cache
        mapping_df.to_csv(output_map_file, index=False)
        print(f"Mapping saved to {output_map_file}")

    # Stats
    total = len(mapping_df)
    mapped = mapping_df['symbol'].notna().sum()
    unmapped = total - mapped
    
    print(f"Total Genes: {total}")
    print(f"Mapped Genes: {mapped}")
    print(f"Unmapped Genes: {unmapped}")

if __name__ == "__main__":
    INPUT_FILE = "data/raw/GSE215865/GSE215865_rnaseq_logCPM_matrix.csv.gz"
    OUTPUT_MAP = "data/processed/gene_id_mapping.csv"
    
    # Ensure processed dir exists
    os.makedirs("data/processed", exist_ok=True)
    
    harmonize_gene_ids(INPUT_FILE, OUTPUT_MAP)
