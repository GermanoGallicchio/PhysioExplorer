
**PhysioExplorer** is a set of functions to extract patterns from multivariate data. Multivariate data refer to any dataset with multiple variables (e.g., many EEG time-frequency points). When the number of variables is much larger than the number of observations and when the variables are highly correlated, the data are sometimes called "megavariate" (see Eriksson et al., 2013). PhysioExplorer can perform any combination of _analysis_ and _objective_ described below in both multivariate and megavariate contexts. 

## pe_cfg.analysis (Analysis)
### 'cluster'
**Cluster-level analysis** (Groppe et al., 2011; Maris & Oostenveld, 2007) is a two-step procedure: (1) compute univariate test statistics at each location and threshold them, (2) form spatial/temporal/spectral clusters of suprathreshold points. Clusters can be defined in a 3-dimensional space (e.g, time-frequency-channel, frequency-frequency-channel) or a lower-dimensional subset (e.g., time-channel, time-frequency, frequency-channel, time). At the heart of the code is a_cluster forming algorithm that combines adjacency criteria (e.g., spatial-temporal-spectral) with the results of univariate statistical testing (e.g., p-values). The code forms clusters on the observed data and, depending on the _objective_ many sets of surrogate data artificially created under the null hypothesis of exchangeability of group/condition labels (permutataion) or many replicates, each with sampling variability, of the original data (bootstrap). The surrogate data are sampled through the Monte-Carlo approach. 

### 'svd_scm'
**SVD-based Symmetric Covariance Mapping** handles multi/megavariate data structures natively (in one step) to find associations between two sets of variables. SVD stands for singular value decomposition. (Note: I think the term SVD-SCM is more representative of the underlying maths and therefore field-agnostic, but (a) nobody apart from me uses it, (b) it does not give justice to its history, and (c) the same mathematical solution can be reached without SVD but more painfully through iterative least squares approximations. More historically-aware and better-known, although field-specific, terms are PLS-C (i.e., Partial Least Squares - Correlation; _REFERENCE_)






## pe_cfg.objective (Objective)
### 'permutation'

**Permutation** is for null-hypothesis testing. In each Monte-Carlo iteration, group/condition labels are shuffled with participants sampled with replacement (bootstrap) or not (permutation), and the statistics are recomputed. The code compares the observed cluster metrics (e.g., cluster mass, singular value) with the null-hypothesis distribution to evaluate their statistical significance. (Note: for cluster analysis, inference is done at the cluster level and not at the point level.)
### 'bootstrap'.

**Bootstrap** is for stability estimation.

```mermaid
  graph LR;
    A(pe_cfg.analysis)
    B(cluster)
    C(svd_scm)
    D(pe_cfg.objective)
    E(permutationH0testing)
    F(bootstrapStability)

    A-->B;
    A-->C;
    D-->E;
    D-->F;

```









### Wish list (future updates)
- mediation
- 4d space (time x 3d MNI space)
- alternative cluster forming not based on statistical significance
- write tutorials on how to use PhysioExplorer



## Cite as
...DOI to arrive...
https://github.com/GermanoGallicchio/PhysioExplorer


# References
Eriksson, L., Byrne, T., Johansson, E., Trygg, J., & Vikström, C. (2013). Multi-and megavariate data analysis basic principles and applications. Umetrics Academy.

Groppe, D. M., Urbach, T. P., & Kutas, M. (2011). Mass univariate analysis of event‐related brain potentials/fields I: A critical tutorial review. Psychophysiology, 48(12), 1711-1725.

Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG-and MEG-data. Journal of neuroscience methods, 164(1), 177-190.
