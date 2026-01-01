| Diagnostic            | Real Data Outcome                 | Interpretation                                                                    |
|:----------------------|:----------------------------------|:----------------------------------------------------------------------------------|
| Parameter Sensitivity | Sparsity=0.95, Jaccard=0.21       | Stable sparse networks in valid parameter regime.                                 |
| Permutation Control   | Jaccard < 0.01                    | Network structure relies entirely on temporal order.                              |
| Delay Structure       | Mean Delay ~9.42 steps            | Inferred delays match STDP window (Tau=10).                                       |
| Negative Controls     | 10% False Discovery Rate (approx) | Random/Shuffled genes do not dominate regulators.                                 |
| Cohort Consistency    | Regulator Corr=0.40               | Functional (regulator) consistency is higher than topological (edge) consistency. |
| Cross-Dataset         | Sparsity ~0.95                    | Algorithm behavior is consistent across species/scales.                           |