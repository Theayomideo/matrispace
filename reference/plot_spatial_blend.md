# Create blended spatial plot for two features

Displays two features on the same spatial plot using a 2D color
gradient.

## Usage

``` r
plot_spatial_blend(
  seurat_obj,
  features,
  feature1_name = NULL,
  feature2_name = NULL,
  bottom_left = "#000000",
  bottom_right = "#FF0000",
  top_left = "#00FF00",
  top_right = "#FFFF00"
)
```

## Arguments

- seurat_obj:

  A Seurat object

- features:

  Character vector of exactly 2 feature names

- feature1_name:

  Display name for feature 1

- feature2_name:

  Display name for feature 2

- bottom_left:

  Color for low/low (default: black)

- bottom_right:

  Color for high feature1/low feature2 (default: red)

- top_left:

  Color for low feature1/high feature2 (default: green)

- top_right:

  Color for high/high (default: yellow)

## Value

A patchwork plot with individual plots and blended result

## Examples

``` r
if (FALSE) { # \dontrun{
p <- plot_spatial_blend(seurat_obj, c("COL1A1", "FN1"))
} # }
```
