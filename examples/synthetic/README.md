# Synthetic structure benchmark

This example is a deterministic statistical validation dataset used to test population-structure correction in PanGWASFlow.

Run it with:

```bash
pangwasflow demo --outdir results/demo --seed 42 --samples 400 --features 120 --pcs 2
```

The benchmark deliberately contains one true phenotype-associated feature, lineage-correlated neutral features, and null features. The expected behavior is that lineage-driven associations dominate the naive scan, then collapse after PCA-based structure adjustment while the true feature remains significant.

Generated result previews are published automatically from the validated workflow after changes reach `main`.
