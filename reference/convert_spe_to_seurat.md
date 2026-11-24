# Convert SpatialExperiment to Seurat object

Converts a SpatialExperiment object to a Seurat object with spatial
image data preserved. Handles multi-sample objects by selecting the
first sample.

## Usage

``` r
convert_spe_to_seurat(
  spe,
  spatial_cols = c(tissue = "in_tissue", row = "array_row", col = "array_col", imagerow =
    "pxl_row_in_fullres", imagecol = "pxl_col_in_fullres"),
  symbol_col = NULL,
  verbose = TRUE
)
```

## Arguments

- spe:

  A SpatialExperiment object

- spatial_cols:

  Named character vector mapping spatial columns:

  - tissue: Column for tissue/in_tissue flag

  - row: Array row coordinate

  - col: Array column coordinate

  - imagerow: Pixel row coordinate (required)

  - imagecol: Pixel column coordinate (required)

- symbol_col:

  Column in rowData containing gene symbols. Auto-detected if NULL.

- verbose:

  Print progress messages

## Value

A Seurat object with VisiumV1 image

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- convert_spe_to_seurat(spe)
} # }
```
