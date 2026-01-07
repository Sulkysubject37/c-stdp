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
