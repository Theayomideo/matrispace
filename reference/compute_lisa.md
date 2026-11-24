# Compute LISA statistics for feature pairs

Calculates Local Indicators of Spatial Association (LISA) to identify
spatial clustering patterns between two features.

## Usage

``` r
compute_lisa(
  seurat_obj,
  feature1,
  feature1_type,
  feature2 = NULL,
  feature2_type = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- feature1:

  First feature name (gene or metadata column)

- feature1_type:

  Type of feature1: "gene", "matrisome gene", "signature", or "any gene"

- feature2:

  Second feature name (optional, for bivariate LISA)

- feature2_type:

  Type of feature2

## Value

Data frame with coordinates, feature values, and LISA classifications

## Details

LISA classifications are based on:

- High.High: High feature1, high neighbor feature2

- High.Low: High feature1, low neighbor feature2

- Low.High: Low feature1, high neighbor feature2

- Low.Low: Low feature1, low neighbor feature2

## Examples

``` r
if (FALSE) { # \dontrun{
lisa_results <- compute_lisa(seurat_obj, "COL1A1", "gene", "COL1A2", "gene")
} # }
```
