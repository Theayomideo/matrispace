# Find spatially corrected ECM markers between two groups

Runs a two-stage workflow designed for within-sample comparisons where
spatial autocorrelation matters. ECM genes are first screened with a
Wilcoxon rank-sum test, then significant hits are re-tested with a spaMM
spatial mixed model using Matern covariance.

## Usage

``` r
find_spatially_corrected_ecm_markers(
  seurat_obj,
  group.by,
  ident.1,
  ident.2,
  assay = NULL,
  layer = "data",
  min.pct = 0.05,
  logfc.threshold = 0.25,
  wilcox.padj.threshold = 0.05,
  categories = NULL,
  sp_topfrac = 1,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object

- group.by:

  Metadata column defining the comparison groups

- ident.1:

  Group of interest

- ident.2:

  Reference group

- assay:

  Assay to use. If NULL, auto-detects.

- layer:

  Expression layer to test. Defaults to `"data"`.

- min.pct:

  Minimum fraction of spots expressing the gene in either group.
  Defaults to `0.05`.

- logfc.threshold:

  Minimum absolute log2 fold-change required to send a gene to the
  spatial model. Defaults to `0.25`.

- wilcox.padj.threshold:

  Adjusted p-value threshold for Wilcoxon pre-screening. Defaults to
  `0.05`.

- categories:

  Optional matrisome categories to restrict testing to.

- sp_topfrac:

  Optional fraction of Wilcoxon-significant genes to pass to spaMM,
  ranked by adjusted p-value. Defaults to `1`.

- verbose:

  Print progress messages

## Value

Data frame containing Wilcoxon statistics, spaMM statistics, and
matrisome annotations.

## Details

This function mirrors the differential expression strategy used in the
MatriSpace case study: Wilcoxon pre-screening on ECM genes followed by
spatial correction with
`spaMM::fitme(gene_expr ~ group + Matern(1 | x + y))`.

## Examples

``` r
if (FALSE) { # \dontrun{
spatial_de <- find_spatially_corrected_ecm_markers(
  seurat_obj,
  group.by = "CompositionCluster_CC",
  ident.1 = "CC1",
  ident.2 = "CC3"
)
} # }
```
