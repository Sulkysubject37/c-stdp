# Immune Gene Sets for COVID-19 Host Response Analysis
# Curated for International Symposium on Infectious Disease
# Sources: KEGG (hsa05171), Reactome (R-HSA-168256), Literature (Zhou et al., 2020)

INTERFERON_RESPONSE = {
    "IFNB1", "IFNA1", "IFNA2", "IFNE", "IFNK",  # Type I Interferons
    "IFNL1", "IFNL2", "IFNL3",                  # Type III Interferons
    "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", # ISGs
    "IFIH1", "DDX58", "DHX58",                  # Viral Sensors (MDA5, RIG-I)
    "STAT1", "STAT2", "IRF3", "IRF7", "IRF9"    # Signaling
}

INFLAMMATORY_RESPONSE = {
    "NFKB1", "RELA", "RELB", "NFKB2",           # NF-kB Family
    "TNF", "IL6", "IL1B", "IL1A", "IL18",       # Pro-inflammatory Cytokines
    "CCL2", "CCL3", "CCL5", "CXCL10", "CXCL8",  # Chemokines
    "JAK1", "JAK2", "TYK2",                     # JAK-STAT
    "NLRP3", "PYCARD"                           # Inflammasome
}

VIRAL_ENTRY_FACTORS = {
    "ACE2", "TMPRSS2", "CTSL", "FURIN", "BSG"   # SARS-CoV-2 Entry
}

# Combined Universe
IMMUNE_GENE_UNIVERSE = INTERFERON_RESPONSE | INFLAMMATORY_RESPONSE | VIRAL_ENTRY_FACTORS

def get_immune_universe() -> set:
    """Returns the set of all curated immune genes."""
    return IMMUNE_GENE_UNIVERSE
