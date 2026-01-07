import gseapy as gp
import os
import sys

def define_immune_universe(output_file):
    print("Fetching gene sets from gseapy...")
    
    libraries = ['KEGG_2021_Human', 'Reactome_2022', 'GO_Biological_Process_2023']
    immune_keywords = [
        'immune', 'inflammation', 'inflammatory', 'cytokine', 
        'interferon', 'interleukin', 'nf-kappab', 'leukocyte', 
        'lymphocyte', 't cell', 'b cell', 'macrophage', 
        'innate', 'adaptive', 'viral', 'bacterial', 'defense'
    ]
    
    immune_genes = set()
    
    for lib in libraries:
        try:
            print(f"Fetching {lib}...")
            gs = gp.get_library(name=lib)
            
            # Filter pathways
            count = 0
            for term, genes in gs.items():
                term_lower = term.lower()
                if any(k in term_lower for k in immune_keywords):
                    immune_genes.update(genes)
                    count += 1
            print(f"  Found {count} immune-related pathways in {lib}")
            
        except Exception as e:
            print(f"  Failed to fetch {lib}: {e}")

    print(f"Total unique immune universe size: {len(immune_genes)}")
    
    if len(immune_genes) < 100:
        print("WARNING: Immune universe seems too small. Check connection or keywords.")
        # Fallback/Error handling would go here, but strict mode says STOP if source unavailable.
        # However, getting <100 might just mean filters are too strict or fetch failed silently.
        
    # Save
    with open(output_file, 'w') as f:
        for gene in sorted(list(immune_genes)):
            f.write(f"{gene}\n")
    
    print(f"Immune universe saved to {output_file}")

if __name__ == "__main__":
    OUTPUT_FILE = "data/processed/immune_universe.txt"
    os.makedirs("data/processed", exist_ok=True)
    define_immune_universe(OUTPUT_FILE)
