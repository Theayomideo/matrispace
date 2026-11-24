# Create matrisome category spatial plots

Generates spatial plots for a matrisome category with robust and
log-scaled scores.

## Usage

``` r
plot_matrisome(seurat_obj, category_name, internal_name, pt.size = 1.5)
```

## Arguments

- seurat_obj:

  A Seurat object with matrisome scores

- category_name:

  Display name (e.g., "Collagens")

- internal_name:

  Internal column prefix

- pt.size:

  Point size

## Value

A patchwork object combining both plots, or NULL if data unavailable

## Examples

``` r
if (FALSE) { # \dontrun{
p <- plot_matrisome(seurat_obj, "Collagens", "collagens")
} # }
```
