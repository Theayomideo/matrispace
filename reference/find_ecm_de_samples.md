# Differential expression of ECM genes between samples or conditions

Compares ECM gene expression between different samples or experimental
conditions. Supports both spot-level and pseudo-bulk analysis
approaches.

## Usage

``` r
find_ecm_de_samples(
  seurat_obj,
  condition.col,
  condition.1,
  condition.2,
  cluster.col = NULL,
  assay = NULL,
  test.use = "wilcox",
  min.pct = 0.05,
  logfc.threshold = 0.25,
  categories = NULL,
  pseudobulk = FALSE,
  sample.col = NULL,
  blocking.col = NULL,
  min.cells.per.sample = 10,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object

- condition.col:

  Metadata column containing condition/sample labels

- condition.1:

  First condition to compare

- condition.2:

  Second condition to compare

- cluster.col:

  Optional metadata column for cluster-aware analysis. In spot-level
  mode, DE is run separately within each cluster. In pseudobulk mode,
  counts are aggregated at the sample-by-cluster level and cluster is
  added to the DESeq2 design.

- assay:

  Assay to use. If NULL, auto-detects.

- test.use:

  Statistical test: "wilcox" (default), "MAST", "DESeq2", "LR"

- min.pct:

  Minimum fraction expressing (default: 0.05)

- logfc.threshold:

  Minimum log2 fold change (default: 0.25)

- categories:

  Matrisome categories to include. If NULL, uses all ECM genes.

- pseudobulk:

  If TRUE, aggregates to pseudo-bulk per sample before DE. Recommended
  for multi-sample experiments.

- sample.col:

  Required if pseudobulk=TRUE. Column identifying individual samples.

- blocking.col:

  Optional metadata column added to the DESeq2 design in pseudobulk mode
  to control for paired subjects, batches, or other blocking effects.

- min.cells.per.sample:

  Minimum cells per sample for pseudo-bulk (default: 10)

- verbose:

  Print progress messages

## Value

Data frame with DE results. If cluster.col is provided, includes a
'cluster' column for cluster-stratified results.

## Details

### Analysis Modes

**Spot-level (pseudobulk=FALSE):** Treats each spot as independent
observation. Fast but may inflate significance due to pseudo-replication
within samples.

**Pseudo-bulk (pseudobulk=TRUE):** Aggregates counts by sample, then
compares between conditions with DESeq2. When `cluster.col` is supplied,
profiles are aggregated at the sample-by-cluster level and cluster is
included in the design. Use `blocking.col` to adjust for repeated
subjects or batch effects.

### Cluster-Stratified Analysis

When cluster.col is provided, DE is run separately within each cluster.
This reveals cluster-specific expression changes between conditions.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic comparison between conditions
de_results <- find_ecm_de_samples(seurat_obj,
                                   condition.col = "condition",
                                   condition.1 = "tumor",
                                   condition.2 = "normal")

# Cluster-stratified analysis
stratified_de <- find_ecm_de_samples(seurat_obj,
                                      condition.col = "condition",
                                      condition.1 = "tumor",
                                      condition.2 = "normal",
                                      cluster.col = "cell_type")

# Pseudo-bulk analysis (recommended for multiple samples)
pseudobulk_de <- find_ecm_de_samples(seurat_obj,
                                      condition.col = "condition",
                                      condition.1 = "treated",
                                      condition.2 = "control",
                                      pseudobulk = TRUE,
                                      sample.col = "sample_id",
                                      blocking.col = "patient_id")
} # }
```
