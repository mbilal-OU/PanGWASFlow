# Pneumococcal penicillin-resistance case study

This case study is PanGWASFlow's first published-data validation beyond synthetic benchmarks.

## Dataset

The source dataset is the public pyseer tutorial collection for penicillin resistance in *Streptococcus pneumoniae*:

- Figshare DOI: `10.6084/m9.figshare.7588832`
- Figshare article ID: `7588832`
- Figshare version resolved by CI: `10.6084/m9.figshare.7588832.v1`
- License: CC BY 4.0
- Published archive: `pyseer_tutorial.tar.bz2`
- Archive size: 395,485,754 bytes
- Published MD5: `d7d5db5f931d7d8233199d416355f2d1`
- Genomes represented in the accessory-gene matrix: 616
- Phenotyped isolates: 603
- Accessory-gene families in the Rtab input: 10,944

PanGWASFlow does not vendor the source dataset. GitHub Actions resolves the Figshare record, verifies the published archive checksum and size, and extracts only the phenotype, accessory-gene matrix, and Mash distance matrix required for this validation path.

## Scientific scope

The current case-study layer evaluates accessory-gene association with binary penicillin resistance. It is a reproducibility and software-validation example, not a claim of a new resistance mechanism.

The analysis uses the published `mash.tsv` genomic distance matrix to estimate population structure by classical multidimensional scaling. Eight MDS components are included as fixed-effect covariates, matching the structure dimensionality used in the reference pyseer tutorial. The association implementation and feature filtering are PanGWASFlow's own, so numerical values are not expected to be identical to pyseer.

The published pyseer tutorial also reports that its SNP analysis recovers the established penicillin-binding-protein loci `pbp2x`, `pbp1a`, and `pbp2b`. Those provide positive controls for a later SNP-coordinate validation layer.

## PanGWASFlow validation stages

1. Resolve and verify the public Figshare record.
2. Verify archive size and MD5 before analysis.
3. Extract `resistances.pheno`, `gene_presence_absence.Rtab`, and `mash.tsv` only.
4. Align the 603 phenotyped isolates to the 616-genome feature and distance inputs.
5. Preserve Rtab missing-value semantics and apply explicit feature QC.
6. Retain features between 1% and 99% prevalence.
7. Compute eight classical-MDS population-structure covariates from Mash distances.
8. Run unadjusted Fisher and structure-adjusted logistic association scans.
9. Report genomic inflation, multiple-testing adjusted results, and complete machine-readable association tables.
10. Generate figures only from the audited CI outputs.

## Interpretation guardrails

Accessory-gene hits are association signals, not proof of causality. Residual genomic inflation, linked gene families, lineage effects, phenotype quality, and annotation context must be evaluated before biological interpretation. The README will only display real-data figures after the CI outputs pass that audit.

## References

- Lees JA et al. Sequence element enrichment analysis to determine the genetic basis of bacterial phenotypes. *Nature Communications* (2016). DOI: `10.1038/ncomms12797`.
- Lees JA et al. pyseer: a comprehensive tool for microbial pangenome-wide association studies. *Bioinformatics* (2018). DOI: `10.1093/bioinformatics/bty539`.
- Pyseer GWAS tutorial: `https://pyseer.readthedocs.io/en/master/tutorial.html`.
