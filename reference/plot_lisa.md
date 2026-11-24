# Create LISA clustering plot

Creates a split-view visualization showing spatial LISA clusters and
quantitative summary bars.

## Usage

``` r
plot_lisa(
  lisa_data,
  selected_types = NULL,
  ann_col = "id",
  f1_name = "F1",
  f2_name = "F2"
)
```

## Arguments

- lisa_data:

  Data frame from `compute_lisa`

- selected_types:

  Annotation types to include (NULL for all)

- ann_col:

  Column name for annotations

- f1_name:

  Feature 1 display name

- f2_name:

  Feature 2 display name

## Value

A patchwork plot

## Examples

``` r
if (FALSE) { # \dontrun{
lisa_data <- compute_lisa(seurat_obj, "COL1A1", "gene", "COL1A2", "gene")
p <- plot_lisa(lisa_data, f1_name = "COL1A1", f2_name = "COL1A2")
} # }
```
