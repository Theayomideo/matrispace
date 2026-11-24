# Prepare Seurat object for matrisome analysis

Master preprocessing function that checks and adds required components:

1.  SCTransform (if needed for ScType)

2.  Matrisome feature scores (UCell)

3.  ECM niche classification (ScType)

4.  ECM niche signatures

## Usage

``` r
prepare_object(
  seurat_obj,
  mat_feat_sigs = NULL,
  ecm_ucell_sigs = NULL,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object

- mat_feat_sigs:

  List of matrisome gene signatures (defaults to package data)

- ecm_ucell_sigs:

  List of ECM niche signatures (defaults to package data)

- verbose:

  Print progress messages

## Value

Preprocessed Seurat object

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- prepare_object(seurat_obj)
} # }
```
