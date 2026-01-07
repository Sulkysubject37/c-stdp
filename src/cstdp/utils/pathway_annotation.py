# Pathway Annotation Layer for Symposium Visualization
# Maps gene lists to static pathway definitions (KEGG/Reactome)
# Strictly observational; no enrichment statistics generated.

from typing import List, Dict, Set

# Static Pathway Definitions (Curated for COVID-19 Context)
PATHWAY_DB = {
    "KEGG_hsa05171_Coronavirus_disease_COVID19": {
        "ACE2", "TMPRSS2", "FURIN", "CTSL", "NFKB1", "RELA", "IL6", "TNF", "CXCL10", "STAT1", "IFNB1"
    },
    "KEGG_hsa04668_TNF_signaling_pathway": {
        "TNF", "TNFRSF1A", "TRADD", "TRAF2", "MAP3K7", "NFKB1", "RELA", "CCL2", "IL1B", "IL6"
    },
    "KEGG_hsa04630_JAK_STAT_signaling_pathway": {
        "JAK1", "JAK2", "TYK2", "STAT1", "STAT2", "IFNA1", "IFNB1", "IL6"
    },
    "REACTOME_R-HSA-168256_Immune_System": {
        "NFKB1", "RELA", "IL6", "TNF", "IL1B", "CXCL8", "OAS1", "MX1", "ISG15", "IFIH1", "DDX58"
    }
}

class PathwayAnnotator:
    @staticmethod
    def annotate_genes(gene_list: List[str]) -> Dict[str, List[str]]:
        """
        Intersects input gene list with curated pathways.
        Returns: Dictionary {Pathway_Name: [Found_Genes]}
        """
        results = {}
        gene_set = set(gene_list)
        
        for pathway, members in PATHWAY_DB.items():
            intersection = gene_set.intersection(members)
            if intersection:
                results[pathway] = list(intersection)
                
        return results

    @staticmethod
    def print_summary_table(annotations: Dict[str, List[str]]):
        """Prints a markdown-compatible table of annotations."""
        print("| Pathway | Overlap Count | Genes |")
        print("|:---|:---:|:---|")
        for pathway, genes in annotations.items():
            print(f"| {pathway} | {len(genes)} | {', '.join(sorted(genes))} |")
