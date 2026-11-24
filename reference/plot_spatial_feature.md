# Create spatial feature plot

Generates a spatial feature plot using ggplot2 with RdYlBu color
palette.

## Usage

``` r
plot_spatial_feature(
  seurat_obj,
  feature,
  annotation_col = NULL,
  pt.size = 1.5,
  title = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object

- feature:

  Feature name (gene or metadata column)

- annotation_col:

  Annotation column for hover text (optional)

- pt.size:

  Point size (default: 1.5)

- title:

  Plot title (defaults to feature name)

## Value

A ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
p <- plot_spatial_feature(seurat_obj, "COL1A1")
} # }
```
