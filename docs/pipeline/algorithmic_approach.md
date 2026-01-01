# C-STDP: Algorithmic Approach

## Core Concept
C-STDP (Causal Spike-Timing Dependent Plasticity) is a biologically inspired method for inferring directed, causal gene regulatory networks (GRNs) from transcriptomic time-series data. Unlike correlation-based methods (e.g., Pearson, Granger) that rely on continuous intensity matching, C-STDP transforms gene expression profiles into discrete "events" (spikes) and infers causality based on the precise temporal ordering of these events.

The core hypothesis is that **upregulation of a regulator must consistently precede the upregulation of its target**.

## 1. Spike Encoding (Event Generation)
Continuous gene expression data $X_i(t)$ is converted into a discrete spike train $S_i$ for each gene $i$. We focus on the *onset* of activation rather than the steady-state level.

**Algorithm:**
1.  Compute the temporal derivative $\frac{dX_i}{dt}$ (rate of change).
2.  Calculate an adaptive threshold $\theta_i$ for each gene:
    $$ \theta_i = \mu(\frac{dX_i}{dt}) + \sigma \cdot \text{std}(\frac{dX_i}{dt}) $$ 
    where $\sigma$ is a sensitivity parameter (default 1.5).
3.  Generate a spike at time $t_k$ if the derivative exceeds the threshold:
    $$ S_i = \{t_k \mid \frac{dX_i}{dt}(t_k) > \theta_i\} $$ 

This acts as a high-pass temporal filter, capturing sharp "bursts" of transcriptional activity while ignoring slow drifts or constant high expression.

## 2. C-STDP Learning Rule
We construct a fully connected weighted graph where $w_{ij}$ represents the strength of the causal link $i \to j$. Weights are updated based on the relative timing of spikes between gene pairs.

**Pairwise Update Rule:**
For every pair of spikes $(t_i, t_j)$ occurring within a temporal window:
-   Let $\Delta t = t_j - t_i$ (Time of Post - Time of Pre).

$$ \Delta w_{ij} = \begin{cases} A_+ \exp(-\frac{|\Delta t|}{\tau_+}) & \text{if } \Delta t > 0 \text{ (Causal: } i \text{ before } j) \\ -A_- \exp(-\frac{|\Delta t|}{\tau_-}) & \text{if } \Delta t < 0 \text{ (Anti-Causal: } j \text{ before } i) \\ 0 & \text{if } \Delta t = 0 \end{cases} $$ 

**Parameters:**
-   $A_+$: Amplitude of potentiation (Reward for correct order).
-   $A_-$: Amplitude of depression (Penalty for wrong order). We typically set $A_- > A_+$ to enforce sparsity.
-   $\tau_+, \tau_-$: Time constants defining the "causal window". Interaction implies delays must be biologically plausible (e.g., within $\tau$).

## 3. Network Inference
1.  Initialize weights $W$ to zero.
2.  Iterate through all spike pairs in the dataset.
3.  Accumulate $\Delta w$ for each edge.
4.  Normalize weights:
    $$ w_{ij} \leftarrow \text{clip}(w_{ij}, 0, w_{max}) $$ 
5.  Thresholding: Edges with $w_{ij} > w_{thresh}$ are retained as the final GRN.

## Why this works
-   **Noise Robustness:** Random noise rarely produces consistent ordered pairs across many samples. C-STDP effectively integrates evidence over time.
-   **Directionality:** The asymmetry of the rule ($i \to j$ vs $j \to i$) explicitly resolves direction, which correlation matrices cannot do.
-   **Sparsity:** The strong depression term ($A_-$) actively pruning "reverse" or "coincident" edges leads to sparse, interpretable networks.
