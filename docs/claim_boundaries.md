# C-STDP Claim Boundaries

## 1. Validated Claims
- **Temporal Sensitivity:** C-STDP inferred networks are provably dependent on the specific temporal ordering of samples. Randomized data yields no structure.
- **Behavioral Robustness:** The algorithm produces consistent network properties (sparsity, regulator ranking) across parameter perturbations, independent cohorts, and different datasets.
- **Noise Rejection:** C-STDP effectively filters out uncorrelated noise (Poisson) and reduces the influence of non-temporal correlations (shuffled profiles) compared to real data.
- **Interpretability:** Every inferred edge is supported by a traceable log of temporal coincidence events (\Delta t \approx \tau).

## 2. Limitations (What we cannot claim)
- **Biological Causality:** We cannot claim that inferred edges represent physical molecular interactions (e.g., binding) without experimental validation (ChIP-seq, perturbation). The "causality" is strictly *predictive temporal precedence*.
- **Topological Precision:** The low edge overlap between cohorts (Jaccard < 0.05) indicates that individual edge predictions are noisy and sensitive to sample heterogeneity.
- **Pseudo-Time Validity:** All results assume the sample index correlates with biological progression. If the "pseudo-time" ordering is incorrect, the "causal" inference is invalid.

## 3. Requirements for Upgrading Claims
- **True Time-Series:** Longitudinal data from the *same* subject (not cross-sectional pseudo-time) is needed to validate tracking of dynamics.
- **Perturbation Data:** Knock-out or overexpression experiments are required to prove functional causality beyond observation.
- **Multi-Omic anchoring:** Linking transcriptomic spikes to protein-level events would strengthen the biological plausibility of the delays.
