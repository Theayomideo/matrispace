# Add scale bar to spatial plot

Adds a scale bar to a ggplot spatial plot based on Visium spot spacing.

## Usage

``` r
add_scalebar(p, seurat_obj, bar_length_um = 1000)
```

## Arguments

- p:

  A ggplot object

- seurat_obj:

  The Seurat object used to create the plot

- bar_length_um:

  Scale bar length in micrometers (default: 1000)

## Value

The plot with scale bar added

## Examples

``` r
if (FALSE) { # \dontrun{
p <- plot_spatial_feature(seurat_obj, "COL1A1")
p <- add_scalebar(p, seurat_obj)
} # }
```
