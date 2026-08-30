# Pneumococcal penicillin-resistance case study

This case study is being built as PanGWASFlow's first public-data validation beyond synthetic benchmarks.

## Dataset

The source dataset is the public pyseer tutorial collection for penicillin resistance in *Streptococcus pneumoniae*:

- Figshare DOI: `10.6084/m9.figshare.7588832`
- Figshare article ID: `7588832`
- Population described by the pyseer tutorial: 616 pneumococcal genomes collected in Massachusetts
- Phenotyped isolates reported by the tutorial: 603
- Relevant published inputs include `resistances.pheno`, `gene_presence_absence.Rtab`, `snps.vcf.gz`, `core_genome_aln.tree`, and reference annotation files.

PanGWASFlow does not vendor the source dataset. The GitHub Actions validation resolves the Figshare record at runtime, records the published file manifest, and will download only the inputs required for the analysis being validated.

## Scientific scope

The first case-study layer will focus on a reproducible genotype-to-phenotype association analysis and on diagnostics for population structure, feature frequency, and multiple testing. It is intended as a reproducibility and software-validation example, not as a claim of a new resistance mechanism.

The published pyseer tutorial reports that its SNP analysis recovers the established penicillin-binding-protein loci `pbp2x`, `pbp1a`, and `pbp2b`. These published results provide positive controls for evaluating a future PanGWASFlow SNP-coordinate workflow. Accessory-gene results will be treated as association signals requiring biological interpretation rather than causal proof.

## Planned PanGWASFlow stages

1. Resolve and record the public Figshare manifest.
2. Fetch the smallest published files needed for each validation path.
3. Normalize binary phenotype and genomic feature inputs without changing missing-value meaning.
4. Report sample overlap and feature QC before association testing.
5. Compare baseline and population-structure-adjusted results.
6. Generate data-driven PCA, QQ, association, and top-hit figures.
7. Preserve a compact result preview and exact source provenance in the repository.

## References

- Lees JA et al. Sequence element enrichment analysis to determine the genetic basis of bacterial phenotypes. *Nature Communications* (2016). DOI: `10.1038/ncomms12797`.
- Lees JA et al. pyseer: a comprehensive tool for microbial pangenome-wide association studies. *Bioinformatics* (2018). DOI: `10.1093/bioinformatics/bty539`.
- Pyseer GWAS tutorial: `https://pyseer.readthedocs.io/en/master/tutorial.html`.
