# Calculate spatial matrisome-pair co-expression scores

Computes spatially weighted co-expression scores for unique heterotypic
MatriComDB pairs at each spot. The method measures whether one member of
a pair is expressed in a spot while its partner is expressed in
neighboring spots; it does not infer signaling direction or
communication activity.

## Usage

``` r
score_matrisome_pair_coexpression(
  seurat_obj,
  pair_db = NULL,
  adj_matrix = NULL,
  assay = NULL,
  layer = "data",
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object

- pair_db:

  MatriComDB pair data with `Gene1` and `Gene2` columns identifying the
  two pair members. Defaults to package data.

- adj_matrix:

  Sparse adjacency matrix for spatial neighbors. If NULL, computed from
  coordinates.

- assay:

  Assay to use (default: "SCT")

- layer:

  Layer within assay (default: "data")

- verbose:

  Print progress messages

## Value

Seurat object with a "MATRISOMEPAIR" assay added.

## Details

### Scoring Formula

\$\$Score\_{spot} = \sqrt{Gene2\_{spot}} \times \sum\_{neighbors}
\sqrt{Gene1\_{neighbor}}\$\$

Pair members are alphabetically ordered before scoring so reciprocal
rows collapse to one stable feature name. Homomeric rows are excluded.

**Product structure:** A score is non-zero only when both pair members
are detected in the relevant spot and neighborhood.

**Square Root Transformation:**

- Variance stabilization for skewed expression data

- Reduced influence of highly expressed genes

- Improved robustness to outliers

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- score_matrisome_pair_coexpression(seurat_obj)
} # }
```
