# Getting started with matrispace

## Overview

`matrispace` is the R package companion to MatriSpace. It mirrors the
same workflow structure used in the web app: data input, matrisome
profiling, and feature analysis.

The package is designed for local, reproducible analysis of spatial
transcriptomics data with a focus on extracellular matrix biology.

``` r

library(matrispace)
library(Seurat)
library(ggplot2)
```

## Load data

You can start from a Seurat object or convert a SpatialExperiment object
with
[`convert_spe_to_seurat()`](https://Theayomideo.github.io/matrispace/reference/convert_spe_to_seurat.md).

``` r

seurat_obj <- readRDS("path/to/your/spatial_object.rds")
```

For examples and testing, the package also ships with a small Seurat
object:

``` r

test_path <- system.file("testdata", "test_seurat.rds", package = "matrispace")
seurat_obj <- readRDS(test_path)
```

## Prepare the object

[`prepare_object()`](https://Theayomideo.github.io/matrispace/reference/prepare_object.md)
runs the standard MatriSpace preprocessing steps:

- gene/signature translation
- matrisome scoring
- ECM niche classification
- ECM niche scoring

``` r

seurat_obj <- prepare_object(seurat_obj, verbose = TRUE)
```

## Matrisome profiling

The main profiling functions mirror the MatriSpace app workflow.

``` r

seurat_obj <- score_matrisome(seurat_obj)
seurat_obj <- annotate_ecm_niches(seurat_obj)
seurat_obj <- score_matrisome_pair_coexpression(seurat_obj)
```

Built-in reference objects are available directly:

``` r

data(matrisome_db)
data(matrisome_signatures)
data(ecm_signatures)
data(matrisome_pairs)

nrow(matrisome_db)
#> [1] 1029
length(matrisome_signatures)
#> [1] 25
names(ecm_signatures)
#> [1] "Interstitial" "Basement"     "Vascular"
```

## Visualization

You can generate the same style of ECM-focused spatial views used
throughout MatriSpace.

``` r

plot_spatial_feature(seurat_obj, "Collagens")
plot_spatial_blend(seurat_obj, c("COL1A1", "FN1"))
```

## Differential expression

The package now exposes two complementary differential-expression
workflows.

### Within-sample spatial correction

Use
[`find_spatially_corrected_ecm_markers()`](https://Theayomideo.github.io/matrispace/reference/find_spatially_corrected_ecm_markers.md)
when comparing niches or compartments within a tissue section and
spatial autocorrelation should be modeled explicitly.

``` r

spatial_de <- find_spatially_corrected_ecm_markers(
  seurat_obj,
  group.by = "CompositionCluster_CC",
  ident.1 = "CC1",
  ident.2 = "CC3"
)
```

### Cross-sample pseudobulk testing

Use
[`find_ecm_de_samples()`](https://Theayomideo.github.io/matrispace/reference/find_ecm_de_samples.md)
with `pseudobulk = TRUE` when comparing conditions across biological
samples.

``` r

pseudobulk_de <- find_ecm_de_samples(
  seurat_obj,
  condition.col = "condition",
  condition.1 = "Inflamed",
  condition.2 = "Non-inflamed",
  pseudobulk = TRUE,
  sample.col = "sample_id",
  cluster.col = "CompositionCluster_CC",
  blocking.col = "patient_id"
)
```

## Spatial statistics

`matrispace` also provides direct access to local and global spatial
analyses.

``` r

lisa_df <- compute_lisa(seurat_obj, "COL1A1", "gene", "FN1", "gene")
morans <- compute_morans_i(seurat_obj)
```

## Related resources

- [MatriSpace web
  app](https://izzilab-matrispace.share.connect.posit.cloud/)
- [MatriSpace app
  repository](https://github.com/Theayomideo/matrispace-app)
- [The Matrisome
  Project](https://sites.google.com/uic.edu/matrisome/home)
