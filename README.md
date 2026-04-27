# matrispace

<!-- badges: start -->
[![R-CMD-check](https://github.com/Theayomideo/matrispace/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Theayomideo/matrispace/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL (>= 3)](https://img.shields.io/badge/license-GPL%20(%3E%3D%203)-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

R package companion to [MatriSpace](https://github.com/Theayomideo/matrispace-app) for identifying, quantifying, and interpreting spatially-resolved extracellular matrix (ECM) gene expression patterns in spatial transcriptomics data.

The [online MatriSpace app](http://matrinet.shinyapps.io/matrispace) provides an interactive interface over curated datasets and user uploads. The `matrispace` package exposes the same workflow in R for scripted analyses, custom modeling, and figure generation.

## Motivation

The extracellular matrix is a spatially organized signaling environment that shapes tissue architecture, niche identity, inflammation, fibrosis, and tumor progression. Spatial transcriptomics enables the study of ECM gene expression in its native tissue context, but most general-purpose workflows are not built around matrisome organization or ECM niches. `matrispace` brings the MatriSpace analytical framework into R so that matrisome-centered analyses can be reproduced and extended in local pipelines.

## Workflow

The package follows the same three stages as the MatriSpace web application: **data input**, **matrisome profiling**, and **feature analysis**.

### Data input

Accepts `Seurat` objects directly and converts `SpatialExperiment` objects via `convert_spe_to_seurat()`. `prepare_object()` standardizes gene symbols, computes matrisome signatures, and classifies each spot into an ECM niche.

### Matrisome profiling

- `score_matrisome()` — UCell scores for the six matrisome categories and curated subfamilies (Basement Membrane, Laminins, Matricellular, Mucins, Peri-vascular ECM).
- `annotate_ecm_niches()` — ECM niche classification (e.g. Interstitial, Basement membrane), with continuous per-spot niche scores.
- `score_lr_activity()` / `find_lr_enrichment()` — ECM-focused ligand-receptor co-expression within and across ECM niches.

### Feature analysis

- Spatial visualization: `plot_spatial_feature()`, `plot_matrisome()`, `plot_spatial_blend()`, `plot_lisa()`.
- Local spatial association: `compute_lisa()`.
- Global spatial autocorrelation: `compute_morans_i()`.
- Differential expression:
  - `find_spatially_corrected_ecm_markers()` — within-sample, two-stage Wilcoxon screen followed by a `spaMM` Matérn spatial mixed model.
  - `find_ecm_de_samples(pseudobulk = TRUE)` — cross-sample pseudobulk testing with `DESeq2`, with optional cluster and blocking adjustment.

## Installation

Requires R (>= 4.1.0).

```r
install.packages("pak")
pak::pak("Theayomideo/matrispace")
```

## Quick start

```r
library(matrispace)

seurat_obj <- readRDS("my_visium_data.rds")
seurat_obj <- prepare_object(seurat_obj)
seurat_obj <- score_matrisome(seurat_obj)
seurat_obj <- annotate_ecm_niches(seurat_obj)
```

See the package vignettes for full differential-expression, ligand-receptor, and spatial-statistics examples.

## See also

- [MatriSpace web app](http://matrinet.shinyapps.io/matrispace)
- [MatriSpace app repository](https://github.com/Theayomideo/matrispace-app)
- [The Matrisome Project](https://sites.google.com/uic.edu/matrisome/home)
- [MatriCom](https://github.com/Izzilab/MatriCom)
