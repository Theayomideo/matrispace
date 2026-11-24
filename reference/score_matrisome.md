# Calculate matrisome feature scores

Adds UCell-based module scores for the current MatriSpace matrisome
signatures to a Seurat object. Signatures include main categories
(Collagens, Glycoproteins, etc.) and curated subcategories such as
Basement Membrane, Laminins, and Peri-vascular ECM.

## Usage

``` r
score_matrisome(seurat_obj, signature_list = NULL)
```

## Arguments

- seurat_obj:

  A Seurat object

- signature_list:

  Named list of gene vectors. Defaults to package data.

## Value

Seurat object with matrisome scores in metadata

## Details

The function:

1.  Translates gene symbols using HGNChelper

2.  Calculates UCell scores for each signature

3.  Creates display-named columns for UI compatibility

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- score_matrisome(seurat_obj)
} # }
```
