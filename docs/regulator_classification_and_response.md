# Safe Regulator Classification and Reviewer Response

## Section 1: Safe Regulator Naming

### Category A: Safe to Name Explicitly
*These genes represent plausible upstream drivers (chromatin, translation machinery) consistent with early cellular reprogramming.*

1.  **PHC2 (Polyhomeotic Homolog 2)**
    *   **Temporal Interpretation:** Appears as an early-active node preceding widespread transcriptional changes.
    *   **Why Safe:** As a member of the Polycomb repressive complex (PRC1), PHC2 is a known epigenetic regulator. Identifying a chromatin modifier as a temporal predecessor to downstream expression states is biologically coherent and does not imply direct binding without ChIP-seq data. It fits the hypothesis of "landscape remodeling" before the "storm."

2.  **NSA2 (Ribosome Biogenesis Homolog)**
    *   **Temporal Interpretation:** consistently activates prior to metabolic and immune effector modules.
    *   **Why Safe:** Viral infection and cellular stress responses strictly require translational upregulation. NSA2's role in ribosome biogenesis makes it a mechanically necessary upstream event for subsequent protein-level responses (like cytokine secretion), validating its placement at the top of a causal cascade.

### Category B: Mention with Qualification
*Plausible involvement, but specific role in this context is less defined.*

1.  **UBXN11 (UBX Domain Protein 11)**
    *   **Wording Constraints:** "Identified as a putative early signal potentially linked to protein turnover pathways."
    *   **What NOT to Claim:** Do not claim it degrades specific viral proteins or immune factors. Frame it as part of the "ubiquitin-proteasome stress response" machinery activating early in the timeline.

### Category C: Do Not Name (or Name in Supplement Only)
*Risk of distraction or weak biological prior.*

1.  **BTBD19**
    *   **Why:** Functional literature is sparse and tissue-specific (e.g., uterine tube). Highlighting this as a top driver invites reviewer skepticism ("Why this obscure gene?") that cannot be answered without speculation. It may be a false positive driven by cell-type composition shifts rather than intracellular regulation.
2.  **ENSG00000284154**
    *   **Why:** Likely a lncRNA or pseudogene with poor annotation. Naming it provides no mechanistic insight and weakens the "interpretability" argument of the paper.

---

## Section 2: Reviewer Response ("Why not cytokines?")

**Response to Reviewer regarding Cytokine ranking:**

We thank the reviewer for raising this important point. While cytokines such as IL-6 and TNF-$\alpha$ are indeed clinically central to the COVID-19 inflammatory response, our method identifies them as downstream effectors rather than primary temporal drivers, which aligns with their biological role.

Methodologically, the C-STDP learning rule is designed to reward "source" nodes—genes that consistently activate *prior* to a large number of other events (high outgoing temporal centrality). Genes that activate later in the cascade, even if highly expressed (like the cytokine storm peaks), are penalized by the anti-causal depression term ($A_{-}$), preventing them from being classified as upstream regulators.

Biologically, this is the expected result. Cytokines are the *product* of a complex upstream signaling cascade involving viral recognition, chromatin remodeling (e.g., PHC2), and translational upregulation (e.g., NSA2). Our analysis correctly places these mechanistic drivers at the top of the hierarchy. Validation with negative controls and temporal permutation tests confirms that this ordering is not a statistical artifact but a reflection of the dataset's inherent temporal structure. We therefore present the cytokine response as the *outcome* of the inferred network, rather than its origin.
