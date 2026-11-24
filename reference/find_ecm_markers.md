# Find differentially expressed ECM genes between groups

Identifies differentially expressed matrisome genes between two groups
(clusters, conditions, etc.) using Seurat's FindMarkers framework.

## Usage

``` r
find_ecm_markers(
  seurat_obj,
  ident.1,
  ident.2 = NULL,
  group.by = NULL,
  assay = NULL,
  test.use = "wilcox",
  min.pct = 0.05,
  logfc.threshold = 0.25,
  categories = NULL,
  only.pos = FALSE,
  p.adjust.method = "BH",
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object

- ident.1:

  Identity of the first group (cluster name, cell names, or condition)

- ident.2:

  Identity of the second group. If NULL, compares ident.1 vs all other
  cells.

- group.by:

  Metadata column to use for grouping. If NULL, uses active identities.

- assay:

  Assay to use. If NULL, auto-detects (SCT \> Spatial \> first
  available).

- test.use:

  Statistical test: "wilcox" (default), "MAST", "DESeq2", "LR", "t"

- min.pct:

  Minimum fraction of cells expressing the gene in either group
  (default: 0.05)

- logfc.threshold:

  Minimum log2 fold change threshold (default: 0.25)

- categories:

  Matrisome categories to include. Can be divisions ("Core matrisome")
  or subcategories ("Collagens"). If NULL, uses all ECM genes.

- only.pos:

  Only return genes with positive log2FC (upregulated in ident.1)

- p.adjust.method:

  Method for p-value adjustment (default: "BH")

- verbose:

  Print progress messages

## Value

Data frame with columns:

- gene:

  Gene symbol

- p_val:

  Raw p-value

- avg_log2FC:

  Average log2 fold change

- pct.1:

  Fraction of cells in ident.1 expressing the gene

- pct.2:

  Fraction of cells in ident.2 expressing the gene

- p_val_adj:

  Adjusted p-value

- category:

  Matrisome division (Core/Associated)

- subcategory:

  Matrisome subcategory

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare fibroblasts vs all other cells
markers <- find_ecm_markers(seurat_obj, ident.1 = "Fibroblasts",
                            group.by = "cell_type")

# Compare two specific clusters
markers <- find_ecm_markers(seurat_obj, ident.1 = "1", ident.2 = "2",
                            group.by = "seurat_clusters")

# Filter to collagens only
collagen_de <- find_ecm_markers(seurat_obj, ident.1 = "Stroma",
                                 categories = "Collagens")
} # }
```
