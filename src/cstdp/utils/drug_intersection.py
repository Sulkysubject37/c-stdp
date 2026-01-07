# Drug Intersection Layer
# Maps regulators to known DrugBank/ChEMBL targets.
# WARNING: For observational/symposium context only. NOT for clinical decision making.

from typing import List, Dict

# Static Drug Knowledge Base (COVID-19 relevant)
DRUG_DB = {
    "IL6": [("Tocilizumab", "Approved"), ("Sarilumab", "Approved")],
    "IL6R": [("Tocilizumab", "Approved")],
    "JAK1": [("Baricitinib", "Approved"), ("Ruxolitinib", "Investigational")],
    "JAK2": [("Baricitinib", "Approved")],
    "TNF": [("Adalimumab", "Approved"), ("Infliximab", "Approved")],
    "TMPRSS2": [("Camostat", "Investigational"), ("Nafamostat", "Investigational")],
    "ACE2": [("Recombinant ACE2", "Investigational")],
    "CTSL": [("Amantadine", "Investigational")], # Weak evidence, strictly observational
    "FURIN": [("Decanoyl-RVKR-CMK", "Experimental")]
}

class DrugIntersector:
    @staticmethod
    def query_drug_targets(gene_list: List[str]) -> Dict[str, List[tuple]]:
        """
        Checks if genes have known drug interactions in the local DB.
        Returns: {Gene: [(Drug, Status), ...]}}
        """
        hits = {}
        for gene in gene_list:
            if gene in DRUG_DB:
                hits[gene] = DRUG_DB[gene]
        return hits

    @staticmethod
    def report_druggability(hits: Dict[str, List[tuple]]):
        """Prints formatted report."""
        if not hits:
            print("No direct druggable targets found in the provided list.")
            return

        print("\n### Druggability Intersection Report")
        print("Disclaimer: Intersects with known drug targets. Does NOT suggest treatment.\n")
        print("| Target Gene | Drug Candidate | Status |")
        print("|:---|:---|:---|")
        for gene, drugs in hits.items():
            for (drug, status) in drugs:
                print(f"| {gene} | {drug} | {status} |")
