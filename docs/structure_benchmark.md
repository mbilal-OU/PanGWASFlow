# Population-structure benchmark

PanGWASFlow includes a deterministic synthetic stress test designed to demonstrate why microbial GWAS must account for clonal or lineage structure.

The benchmark contains three feature classes:

- `feature_000`: a true phenotype-associated feature generated independently of lineage.
- `feature_001` to `feature_035`: lineage markers with no direct phenotype effect. They become associated in the naive scan because both feature frequency and phenotype prevalence differ between latent lineages.
- remaining features: null features.

Population structure is estimated from the feature matrix by standardized principal-component analysis. The adjusted scan uses logistic regression with the selected PCs as covariates. Benjamini-Hochberg correction is applied to both naive and adjusted scans.

The latent lineage label is retained only because this is a synthetic validation benchmark. It is used to verify that PCA recovers the intended structure and is not passed to the adjusted association model.

The benchmark is a software and statistical validation dataset, not a biological result.
