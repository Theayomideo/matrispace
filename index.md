# matrispace

R package companion to
[MatriSpace](https://github.com/Theayomideo/matrispace-app) for
identifying, quantifying, and interpreting spatially-resolved
extracellular matrix (ECM) gene expression patterns in spatial
transcriptomics data.

If you use `matrispace` in your publications, please cite our paper:
**doi** [10.1242/jcs.264998](https://doi.org/10.1242/jcs.264998)

[![Badge](https://img.shields.io/badge/Paper-Journal_of_Cell_Science-red)](https://doi.org/10.1242/jcs.264998)

The [online MatriSpace
app](https://izzilab-matrispace.share.connect.posit.cloud/) provides an
interactive interface over 240 curated datasets (180 cancer and 60
healthy tissue sections) and user uploads. The `matrispace` package
exposes the same workflow in R for scripted analyses, custom modeling,
and figure generation.

## Motivation

The extracellular matrix is a spatially organized signaling environment
that shapes tissue architecture, niche identity, inflammation, fibrosis,
and tumor progression. Spatial transcriptomics enables the study of ECM
gene expression in its native tissue context, but most general-purpose
workflows are not built around matrisome organization or ECM niches.
`matrispace` brings the MatriSpace analytical framework into R so that
matrisome-centered analyses can be reproduced and extended in local
pipelines.

## Workflow

The package follows the same three stages as the MatriSpace web
application: **data input**, **matrisome profiling**, and **feature
analysis**.

### Data input

Accepts `Seurat` objects directly and converts `SpatialExperiment`
objects via
[`convert_spe_to_seurat()`](https://Theayomideo.github.io/matrispace/reference/convert_spe_to_seurat.md).
[`prepare_object()`](https://Theayomideo.github.io/matrispace/reference/prepare_object.md)
standardizes gene symbols, computes matrisome signatures, and classifies
each spot into an ECM niche.

### Matrisome profiling

- [`score_matrisome()`](https://Theayomideo.github.io/matrispace/reference/score_matrisome.md)
  — UCell scores for the six matrisome categories and curated
  subfamilies (Basement Membrane, Laminins, Matricellular, Mucins,
  Peri-vascular ECM).
- [`annotate_ecm_niches()`](https://Theayomideo.github.io/matrispace/reference/annotate_ecm_niches.md)
  — ECM niche classification (e.g. Interstitial, Basement membrane),
  with continuous per-spot niche scores.
- [`score_matrisome_pair_coexpression()`](https://Theayomideo.github.io/matrispace/reference/score_matrisome_pair_coexpression.md)
  /
  [`find_matrisome_pair_enrichment()`](https://Theayomideo.github.io/matrispace/reference/find_matrisome_pair_enrichment.md)
  — spatial co-expression of unique heterotypic MatriComDB pairs in
  which at least one member is a matrisome gene.

### Feature analysis

- Spatial visualization:
  [`plot_spatial_feature()`](https://Theayomideo.github.io/matrispace/reference/plot_spatial_feature.md),
  [`plot_matrisome()`](https://Theayomideo.github.io/matrispace/reference/plot_matrisome.md),
  [`plot_spatial_blend()`](https://Theayomideo.github.io/matrispace/reference/plot_spatial_blend.md),
  [`plot_lisa()`](https://Theayomideo.github.io/matrispace/reference/plot_lisa.md).
- Local spatial association:
  [`compute_lisa()`](https://Theayomideo.github.io/matrispace/reference/compute_lisa.md).
- Global spatial autocorrelation:
  [`compute_morans_i()`](https://Theayomideo.github.io/matrispace/reference/compute_morans_i.md).
- Differential expression:
  - [`find_spatially_corrected_ecm_markers()`](https://Theayomideo.github.io/matrispace/reference/find_spatially_corrected_ecm_markers.md)
    — within-sample, two-stage Wilcoxon screen followed by a `spaMM`
    Matérn spatial mixed model.
  - `find_ecm_de_samples(pseudobulk = TRUE)` — cross-sample pseudobulk
    testing with `DESeq2`, with optional cluster and blocking
    adjustment.

## Installation

Requires R (\>= 4.1.0).

``` r

install.packages("pak")
pak::pak("Theayomideo/matrispace")
```

## Quick start

``` r

library(matrispace)

seurat_obj <- readRDS("my_visium_data.rds")
seurat_obj <- prepare_object(seurat_obj)
seurat_obj <- score_matrisome(seurat_obj)
seurat_obj <- annotate_ecm_niches(seurat_obj)
```

See the package vignettes for full differential-expression,
matrisome-pair co-expression, and spatial-statistics examples.

## See also

- [MatriSpace web
  app](https://izzilab-matrispace.share.connect.posit.cloud/)
- [MatriSpace app
  repository](https://github.com/Theayomideo/matrispace-app)
- [The Matrisome
  Project](https://sites.google.com/uic.edu/matrisome/home)
- [MatriCom](https://github.com/Izzilab/MatriCom)
