import pandas as pd
import mygene
import os
import sys

def annotate_druggability(adj_file, mapping_file, dataset_name):
    print(f"\n--- Checking Drug Interactions for {dataset_name} ---")
    
    if not os.path.exists(adj_file):
        print(f"File not found: {adj_file}")
        return

    adj_df = pd.read_csv(adj_file, index_col=0)
    # Top 5 regulators
    out_degree = adj_df.sum(axis=1).sort_values(ascending=False).head(5)
    top_ids = out_degree.index.tolist()
    
    # Map to Symbols
    map_df = pd.read_csv(mapping_file)
    ens_to_sym = pd.Series(map_df.symbol.values, index=map_df.ensembl_id).to_dict()
    
    top_symbols = []
    for ident in top_ids:
        clean_id = ident.split('.')[0]
        sym = ens_to_sym.get(clean_id, ident)
        top_symbols.append(sym)
        
    print(f"Top 5 Regulators: {top_symbols}")
    
    mg = mygene.MyGeneInfo()
    # Query 'pharmgkb', 'drugbank', 'chembl' scopes? Actually 'chembl' is not a scope in mygene directly always?
    # We query by symbol, asking for 'drugbank', 'chembl' fields.
    try:
        results = mg.querymany(top_symbols, scopes='symbol', fields='drugbank,chembl', species='human')
        
        found_drugs = False
        for res in results:
            gene = res.get('query')
            db = res.get('drugbank')
            chem = res.get('chembl')
            
            drugs = []
            if db:
                if isinstance(db, list):
                    drugs.extend([d.get('drug_name', 'Unknown') for d in db])
                else:
                    drugs.append(db.get('drug_name', 'Unknown'))
                    
            if chem:
                # ChEMBL usually returns ID. MyGene might not give names.
                # Just noting presence
                drugs.append("ChEMBL Hits")
                
            if drugs:
                found_drugs = True
                # Limit to 3 drugs for display
                print(f"  {gene}: {', '.join(drugs[:3])} ... ({len(drugs)} total)")
            else:
                print(f"  {gene}: No direct drug annotations found.")
                
        if not found_drugs:
            print("WARNING: Sparse drug coverage. No major drug targets found in top 5.")
            # We don't stop the pipeline, just stop this step.
            return

    except Exception as e:
        print(f"Drug annotation error: {e}")

if __name__ == "__main__":
    MAPPING = "data/processed/gene_id_mapping.csv"
    
    print("NOTE: This is an annotation layer only. Presence of a drug interaction does NOT imply therapeutic efficacy for COVID-19.")
    
    annotate_druggability("visuals/real_data/GSE215865_Immune/inferred_grn_adj.csv", MAPPING, "GSE215865")
    annotate_druggability("visuals/real_data/GSE157859_Immune/inferred_grn_adj.csv", MAPPING, "GSE157859")
