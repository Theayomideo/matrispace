# Compute Moran's I spatial autocorrelation

Calculates global Moran's I statistic for matrisome scores to assess
spatial clustering patterns.

## Usage

``` r
compute_morans_i(seurat_obj, features = NULL)
```

## Arguments

- seurat_obj:

  A Seurat object with matrisome scores

- features:

  Character vector of feature names to test. If NULL, automatically
  detects columns ending in "\_robust_score".

## Value

Data frame with Moran's I and p-values for each feature

## Details

Requires the MERINGUE package. If not installed, returns NULL with a
warning.

## Examples

``` r
if (FALSE) { # \dontrun{
autocorr <- compute_morans_i(seurat_obj)
} # }
```
