# Scientific guardrails

PanGWASFlow is designed to make common microbial-GWAS failure modes visible rather than hide them.

## Association is not causation

A significant association does not by itself establish mechanism, direct causality, or clinical relevance. Significant features are candidates for follow-up and should be interpreted with annotation, lineage context, linkage, and independent evidence.

## Population structure is a first-class variable

Clonal population structure can create strong associations when both genotype and phenotype are lineage-correlated. The workflow therefore reports naive and structure-adjusted results separately. A large change after adjustment is itself informative and should not be silently discarded.

## Multiple testing is explicit

Raw P values and Benjamini-Hochberg adjusted q values are both retained. Future modules may add permutation or family-wise error procedures for representations with strong dependency structures.

## Phenotype quality limits inference

Binary labels, MIC values, environmental traits, and host-associated phenotypes have different measurement properties. Missingness, class imbalance, uncertain labels, and inconsistent measurement protocols must be reported before association testing.

## Feature frequency matters

Very rare variants can produce unstable effect estimates. Future biological workflows will expose minimum allele or feature-frequency filters in configuration and report the number of tested features after filtering.

## Genomic representations are not interchangeable

SNPs, accessory genes, k-mers, and unitigs capture different forms of variation and different linkage structures. PanGWASFlow keeps their results separate and compares them at the interpretation stage.

## Reproducibility

Every result should be reproducible from versioned inputs, configuration, software environments, and code. README figures will be generated from committed or CI-produced result tables rather than drawn as illustrative substitutes.
