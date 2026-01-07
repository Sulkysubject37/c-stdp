# C-STDP: Temporal Event Mining of Host Immune Response

## Overview
Understanding the host immune response to infection requires more than a snapshot of differentially expressed genes; it demands a temporal map of activation events. This project pivots the Causal Spike-Timing Dependent Plasticity (C-STDP) algorithm—originally biologically inspired—into a purely statistical tool for mining temporal precedence in gene expression. The goal is to reconstruct the "infection cascade": the precise order in which immune factors activate, allowing us to distinguish early drivers from downstream effectors.

## Input Data
The pipeline accepts bulk RNA-seq expression matrices (LogCPM or TPM). Crucially, the algorithm treats the sample axis as a pseudo-time dimension, assuming that the samples represent a progression of the infection state or a time-series experiment.

## Immune-Aware Preprocessing
To focus the analysis on relevant biological signals and reduce computational noise, we implement a strict immune-aware filtering step:
1.  **Harmonization:** Gene identifiers (Ensembl) are stripped of version suffixes and harmonized to standard Gene Symbols using `mygene`.
2.  **Universe Construction:** We define an "Immune Universe" by aggregating pathways from KEGG, Reactome, and GO Biological Process that contain keywords such as "cytokine", "interferon", "inflammation", and "viral defense".
3.  **Filtering:** Only genes present in this universe are retained for analysis. This reduces the feature space (from ~25k to ~4k genes) to a biologically coherent subset, improving the signal-to-noise ratio for detecting regulatory interactions.

## Sudden Onset Detection
Traditional amplitude-based thresholding often fails in log-transformed RNA-seq data due to dynamic range compression. To address this, we implemented a **Derivative-based Z-score Encoding**:
1.  We compute the temporal derivative (rate of change) of gene expression.
2.  We apply gene-wise Z-score normalization to these derivatives.
3.  A "Sudden Onset Event" is defined as a normalized derivative exceeding a strict threshold (Sigma > 1.5).

This approach identifies genes that undergo a rapid acceleration in expression—interpreted biologically as the moment of transcriptional activation or "burst"—rather than just high absolute abundance.

## C-STDP Learning Rule
The core learning algorithm leverages the principle of **temporal precedence**. For any pair of genes (A and B), we analyze the relative timing of their onset events:
- If Gene A consistently activates *before* Gene B within a short time window, the connection A $\rightarrow$ B is strengthened (Potentiation).
- If Gene A activates *after* Gene B, the connection is weakened (Depression).

This creates a directed weight matrix where the magnitude of a link represents the strength and consistency of the temporal delay between the two genes.

## Network Extraction
From the learned weight matrix, we extract a Gene Regulatory Network (GRN):
- **Regulators:** Defined as genes with high "Out-Degree" (sum of outgoing weights). These represent nodes that consistently precede a large number of other immune factors.
- **Sparsity:** We apply a strict sparsity constraint (e.g., retaining only the top 10-25% of strongest links or subsetting to the most active genes) to filter out spurious correlations and focus on the backbone of the immune response.

## Validation Strategy
The algorithm's reliability is bolstered by:
- **Synthetic Validation:** Recovery of ground-truth networks in simulated data.
- **Permutation Tests:** Verifying that the inferred structure vanishes when temporal order is shuffled.
- **Negative Controls:** Ensuring unrelated housekeeping genes do not appear as top regulators.

## Interpretation Boundaries
While powerful, this method infers **temporal precedence**, which is a proxy for, but not identical to, biological causality.
- It can identify: Likely upstream regulators, activation sequences, and temporal modules.
- It cannot claim: Direct physical binding (e.g., transcription factor binding) or therapeutic efficacy without experimental validation.
- Results are dependent on the resolution of the input time-series; very fast molecular events occurring between sampling points cannot be resolved.

## Summary
By rigorously filtering for immune-relevant genes and focusing on the rate of change rather than absolute levels, C-STDP provides a novel lens for infectious disease data. It moves beyond "what changed" to "when it changed," enabling the reconstruction of the infection cascade from transcriptomic data.
