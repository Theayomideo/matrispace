# Compute expression metrics for matrisome categories

Calculates raw counts, log-scaled, and robust-normalized expression
values for a set of genes within each spot.

## Usage

``` r
compute_expression_metrics(
  gene_list,
  seurat_obj,
  category_name,
  counts_matrix = NULL,
  matrisome_db = NULL
)
```

## Arguments

- gene_list:

  Character vector of gene names

- seurat_obj:

  A Seurat object

- category_name:

  Name for the category (used in column names)

- counts_matrix:

  Pre-fetched counts matrix (optional)

- matrisome_db:

  Matrisome gene database with synonym info

## Value

Data frame with expression metrics and coordinates

## Examples

``` r
if (FALSE) { # \dontrun{
metrics <- compute_expression_metrics(
  gene_list = c("COL1A1", "COL1A2"),
  seurat_obj = my_obj,
  category_name = "Collagens"
)
} # }
```
