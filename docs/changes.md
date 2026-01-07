# C-STDP Symposium Adaptation Log

## Phase 1: Biologically Informed Input Selection
**Date:** 2026-01-08

### Changes
- **Immune Gene Universe:** Created `src/cstdp/utils/immune_gene_sets.py` defining core sets for Interferon, NF-kB, and Viral Entry.
- **Preprocessing Logic:** Modified `scripts/data_prep/preprocess_primary.py` to prioritize `immune_focused` selection over raw variance.

### Rationale (Clinical Relevance)
- Moving from statistical variance to biological relevance ensures the inferred network represents the host immune response to SARS-CoV-2, not just metabolic noise. This is critical for the Infectious Disease Symposium audience.

### Non-Claims

- We do not claim these are the *only* relevant genes, but a curated subset of known drivers.



## Phase 2: Terminology Reframing

**Date:** 2026-01-08



### Changes

- **Mapping Document:** Created `docs/terminology_mapping.md` to map neuroscience jargon (Spike, STDP) to clinical terms (Sudden Onset, Temporal Event Mining).



### Rationale (Clinical Relevance)

- Ensures the methodology is accessible to virologists and clinicians. Avoids confusion with "Spike protein" (viral) vs "Spike" (algorithm).



### Non-Claims



- We are not changing the underlying mathematics, only the presentation layer.







## Phase 3: Downstream Biological Contextualization



**Date:** 2026-01-08







### Changes



- **Pathway Annotation:** Created `src/cstdp/utils/pathway_annotation.py` to map genes to static KEGG/Reactome subsets.



- **Druggability:** Created `src/cstdp/utils/drug_intersection.py` to highlight overlaps with known COVID-19 therapeutics (e.g., Tocilizumab/IL6).







### Rationale (Clinical Relevance)



- Moves the output from "Graph Theory" to "Translational Insight". Answers "What can we do about it?"







### Non-Claims



- **STRICT:** We report *intersection* with drug targets. We do NOT predict efficacy or suggest off-label use. This is hypothesis-generating context only.




