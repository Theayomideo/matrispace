# Find enriched matrisome-pair co-expression per group

Calculates enrichment statistics (avg_log2FC and percentage difference)
for each matrisome pair across cell or spot groups.

## Usage

``` r
find_matrisome_pair_enrichment(
  seurat_obj,
  assay_name = "MATRISOMEPAIR",
  group_by
)
```

## Arguments

- seurat_obj:

  A Seurat object with a MATRISOMEPAIR assay

- assay_name:

  Assay name (default: "MATRISOMEPAIR")

- group_by:

  Metadata column for grouping (e.g., "ecm_niche")

## Value

Data frame with columns: feature, cluster, avg_log2FC, perc_difference,
pct.1, pct.2

## Examples

``` r
if (FALSE) { # \dontrun{
pair_stats <- find_matrisome_pair_enrichment(
  seurat_obj,
  group_by = "ecm_niche"
)
} # }
```
