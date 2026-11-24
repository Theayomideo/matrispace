# Annotate ECM niches using ScType

Classifies each spot into an ECM niche (e.g. Interstitial, Basement
membrane) using the ScType algorithm with ECM-specific marker genes.

## Usage

``` r
annotate_ecm_niches(seurat_obj, markers_path = NULL)
```

## Arguments

- seurat_obj:

  A Seurat object with SCT assay

- markers_path:

  Path to Excel file with ECM niche markers. If NULL, uses built-in
  markers from inst/extdata.

## Value

Seurat object with `ecm_niche` column in metadata (and
`ecm_domain_annotation` retained as an alias for backward
compatibility).

## Details

The function requires scaled SCT data. If the object lacks an SCT assay,
run
[`Seurat::SCTransform()`](https://satijalab.org/seurat/reference/SCTransform.html)
first.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- annotate_ecm_niches(seurat_obj)
table(seurat_obj$ecm_niche)
} # }
```
