# Aggregate reciprocal rows into unique matrisome pairs

Identifies reciprocal rows (e.g., A-B and B-A) and combines them into
one unordered heterotypic pair.

## Usage

``` r
aggregate_matrisome_pairs(stats_df, means_matrix)
```

## Arguments

- stats_df:

  Long-format stats from `find_matrisome_pair_enrichment`

- means_matrix:

  Wide matrix with features as rows, clusters as columns

## Value

List with coalesced `stats` data frame and `means` matrix

## Examples

``` r
if (FALSE) { # \dontrun{
means_matrix <- AverageExpression(
  seurat_obj,
  assays = "MATRISOMEPAIR"
)$MATRISOMEPAIR
coalesced <- aggregate_matrisome_pairs(pair_stats, means_matrix)
} # }
```
