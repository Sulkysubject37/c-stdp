import pandas as pd
import gseapy as gp
import os
import sys

def annotate_regulators(adj_file, mapping_file, output_name):
    print(f"\n--- Annotating Regulators for {output_name} ---")
    if not os.path.exists(adj_file):
        print(f"File not found: {adj_file}")
        return

    adj_df = pd.read_csv(adj_file, index_col=0)
    # Get top 10 regulators by out-degree
    out_degree = adj_df.sum(axis=1).sort_values(ascending=False).head(10)
    top_ids = out_degree.index.tolist()
    
    # Map to Symbols
    map_df = pd.read_csv(mapping_file)
    ens_to_sym = pd.Series(map_df.symbol.values, index=map_df.ensembl_id).to_dict()
    
    top_symbols = []
    for ident in top_ids:
        clean_id = ident.split('.')[0]
        sym = ens_to_sym.get(clean_id, ident) # Fallback to original if not found
        top_symbols.append(sym)
        
    print(f"Top 10 Regulators (Symbols): {top_symbols}")
    
    # Pathway Enrichment for these 10 genes
    # Note: 10 genes is small for standard enrichment, but we can check their participation in pathways
    try:
        enr = gp.enrichr(gene_list=top_symbols,
                         gene_sets=['KEGG_2021_Human', 'WikiPathway_2021_Human', 'Reactome_2022'],
                         organism='human',
                         outdir=None)
        
        results = enr.results
        # Filter for relevant immune keywords
        immune_keywords = ['interferon', 'cytokine', 'nf-kb', 'viral', 'inflammation', 'toll', 't cell', 'b cell']
        
        matches = results[results['Term'].str.lower().str.contains('|'.join(immune_keywords))]
        
        if not matches.empty:
            print("\nKey Pathway Matches for Top Regulators:")
            # Show Term and Genes
            subset = matches[['Term', 'Overlap', 'P-value', 'Genes']].head(10)
            print(subset.to_string(index=False))
        else:
            print("No strong matches with immune keywords in top pathways.")
            # Show top 5 regardless
            print("\nTop 5 General Pathways:")
            print(results[['Term', 'Overlap', 'P-value']].head(5).to_string(index=False))

    except Exception as e:
        print(f"Enrichment error: {e}")

if __name__ == "__main__":
    MAPPING = "data/processed/gene_id_mapping.csv"
    
    # Primary
    annotate_regulators("visuals/real_data/GSE215865_Immune/inferred_grn_adj.csv", MAPPING, "GSE215865")
    
    # Supplementary
    annotate_regulators("visuals/real_data/GSE157859_Immune/inferred_grn_adj.csv", MAPPING, "GSE157859")
