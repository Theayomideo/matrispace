# Translate gene signatures using HGNChelper

Checks gene symbols and corrects outdated/alias symbols to match what's
present in the Seurat object.

## Usage

``` r
translate_signatures(seurat_obj, signature_list)
```

## Arguments

- seurat_obj:

  A Seurat object

- signature_list:

  Named list of gene vectors

## Value

Named list of corrected gene vectors

## Examples

``` r
if (FALSE) { # \dontrun{
sigs <- translate_signatures(seurat_obj, my_signatures)
} # }
```
