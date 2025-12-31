# C-STDP Final Research Report

## 1. What C-STDP provably does
- **Event-Based Causal Inference:** C-STDP successfully identifies directed relationships between genes by focusing on the temporal precedence of discrete activation events (spikes), rather than continuous intensity correlations.
- **Robustness to bursty signals:** The algorithm demonstrates superior performance compared to linear baselines like Granger Causality on sparse, bursty synthetic expression data.
- **Traceable Interpretability:** Every edge in the inferred Gene Regulatory Network (GRN) can be traced back to specific spike-pair coincidences, providing a granular "causal evidence" log for each interaction.
- **Directional Sensitivity:** C-STDP consistently identifies the correct direction of influence ($i \to j$) within its temporal integration window.

## 2. What C-STDP does not claim
- **Biological Truth:** The identified regulators (e.g., in GSE215865) are candidates for hypothesis generation and require experimental validation (e.g., ChIP-seq or knock-out).
- **Absolute Edge Topological Truth:** topological overlap between cohorts remains low, indicating that individual edges are sensitive to noise and data heterogeneity.
- **Completeness:** The current model focuses on pair-based interactions and does not explicitly account for higher-order combinatorial regulation.

## 3. Where causality is valid
- **Temporal Precedence:** Causality is valid under the assumption that a regulator must activate *before* its target within a biologically plausible time window ($\tau$).
- **Synthetic Validation:** Validated on simulated Delay-Differential Equation models where ground truth is known.
- **Permutation Control:** Provably fails when temporal structure is destroyed, confirming the causal inference is driven by time-series order.

## 4. Where assumptions dominate
- **Pseudo-Time:** In real-world datasets like GSE215865, samples are treated as a sequence. This assumes a consistent biological progression across the sample index.
- **Parameter Fixedness:** The model assumes fixed time constants ($\tau$) and amplitudes ($A$) for all gene pairs, whereas real biological delays may vary by mechanism.
- **Spike Thresholding:** The transformation of continuous data to spikes is dependent on the `sigma` parameter, which influences sensitivity vs specificity.

## 5. When results should not be trusted
- **Small Sample Sizes:** In datasets with very few time points, random coincidences may lead to spurious causal links.
- **Highly Correlated Global Trends:** If all genes respond simultaneously to a global stimulus, STDP may struggle to distinguish specific regulatory links from global coincidences.
- **Out-of-Window Delays:** Causal interactions with delays significantly longer than the STDP integration window (e.g., $> 3\tau$) will not be captured.

---
*End of Report*
