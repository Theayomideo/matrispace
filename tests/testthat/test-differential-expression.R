# Tests for ECM differential expression functions

# Helper function to create a mock Seurat object with ECM genes
create_mock_seurat <- function(n_cells = 100, n_clusters = 3, n_conditions = 2) {
  # Get some real ECM genes
  data(matrisome_db, package = "matrispace")
  ecm_genes <- head(matrisome_db$gene, 50)

  # Add some non-ECM genes
  other_genes <- paste0("OtherGene", 1:50)
  all_genes <- c(ecm_genes, other_genes)

  # Create count matrix with differential expression
  set.seed(42)
  counts <- matrix(
    rpois(length(all_genes) * n_cells, lambda = 5),
    nrow = length(all_genes),
    ncol = n_cells,
    dimnames = list(all_genes, paste0("Cell", 1:n_cells))
  )

  # Create cluster assignments
  clusters <- factor(rep(paste0("Cluster", 1:n_clusters), length.out = n_cells))

  # Create condition assignments
  conditions <- factor(rep(c("Control", "Treatment")[1:n_conditions], length.out = n_cells))

  # Create sample assignments (2 samples per condition)
  samples <- factor(paste0("Sample", rep(1:(n_conditions * 2), each = n_cells / (n_conditions * 2))))

  # Add differential expression for some ECM genes in Cluster1
  cluster1_cells <- which(clusters == "Cluster1")
  counts[1:10, cluster1_cells] <- counts[1:10, cluster1_cells] + rpois(10 * length(cluster1_cells), lambda = 20)

  # Add differential expression between conditions for some genes
  treatment_cells <- which(conditions == "Treatment")
  counts[11:20, treatment_cells] <- counts[11:20, treatment_cells] + rpois(10 * length(treatment_cells), lambda = 15)

  # Create Seurat object
  obj <- Seurat::CreateSeuratObject(counts = counts)

  # Normalize to populate the data layer (important for Seurat v5)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)

  obj$seurat_clusters <- clusters
  obj$condition <- conditions
  obj$sample_id <- samples

  # Add mock spatial coordinates
  obj$imagerow <- runif(n_cells, 0, 100)
  obj$imagecol <- runif(n_cells, 0, 100)

  obj
}

# Test helper function
test_that(".get_ecm_genes returns valid gene list", {
  genes <- matrispace:::.get_ecm_genes()
  expect_type(genes, "character")
  expect_gt(length(genes), 1000)

  # Test category filtering
  collagen_genes <- matrispace:::.get_ecm_genes(categories = "Collagens")
  expect_type(collagen_genes, "character")
  expect_gt(length(collagen_genes), 0)
  expect_lt(length(collagen_genes), length(genes))

  # Test invalid category
  expect_warning(
    matrispace:::.get_ecm_genes(categories = "NonexistentCategory"),
    "No genes found"
  )
})

test_that(".annotate_ecm_results adds category columns", {
  df <- data.frame(
    gene = c("COL1A1", "FN1", "MMP2"),
    p_val = c(0.01, 0.02, 0.03)
  )

  result <- matrispace:::.annotate_ecm_results(df)

  expect_true("category" %in% colnames(result))
  expect_true("subcategory" %in% colnames(result))
  expect_equal(nrow(result), 3)
})

test_that(".annotate_ecm_results handles empty input", {
  df <- data.frame(gene = character(), p_val = numeric())
  result <- matrispace:::.annotate_ecm_results(df)
  expect_equal(nrow(result), 0)
})

# Test find_ecm_markers
test_that("find_ecm_markers works with basic input", {
  skip_if_not_installed("Seurat")

  obj <- create_mock_seurat()

  result <- find_ecm_markers(
    obj,
    ident.1 = "Cluster1",
    group.by = "seurat_clusters",
    verbose = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true("gene" %in% colnames(result))
  expect_true("avg_log2FC" %in% colnames(result))
  expect_true("p_val_adj" %in% colnames(result))
  expect_true("category" %in% colnames(result))
})

test_that("find_ecm_markers filters by category", {
  skip_if_not_installed("Seurat")

  obj <- create_mock_seurat()

  # Get all ECM markers
  all_result <- find_ecm_markers(
    obj,
    ident.1 = "Cluster1",
    group.by = "seurat_clusters",
    logfc.threshold = 0,
    min.pct = 0,
    verbose = FALSE
  )

  # Get only Core matrisome (broader category that should match our mock data)
  core_result <- find_ecm_markers(
    obj,
    ident.1 = "Cluster1",
    group.by = "seurat_clusters",
    categories = "Core matrisome",
    logfc.threshold = 0,
    min.pct = 0,
    verbose = FALSE
  )

  # Core result should have fewer or equal genes than all ECM
  expect_lte(nrow(core_result), nrow(all_result))
})

test_that("find_ecm_markers handles pairwise comparison", {
  skip_if_not_installed("Seurat")

  obj <- create_mock_seurat()

  result <- find_ecm_markers(
    obj,
    ident.1 = "Cluster1",
    ident.2 = "Cluster2",
    group.by = "seurat_clusters",
    verbose = FALSE
  )

  expect_s3_class(result, "data.frame")
})

test_that("find_ecm_markers validates input", {
  skip_if_not_installed("Seurat")

  obj <- create_mock_seurat()

  # Invalid Seurat object
  expect_error(
    find_ecm_markers("not_a_seurat", ident.1 = "Cluster1"),
    "must be a Seurat object"
  )

  # Invalid cluster name (error comes from Seurat's FindMarkers)
  expect_error(
    find_ecm_markers(obj, ident.1 = "NonexistentCluster", group.by = "seurat_clusters"),
    "Cannot find|not found"
  )
})

# Test find_all_ecm_markers
test_that("find_all_ecm_markers works with all clusters", {
  skip_if_not_installed("Seurat")

  obj <- create_mock_seurat()

  result <- find_all_ecm_markers(
    obj,
    group.by = "seurat_clusters",
    verbose = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true("cluster" %in% colnames(result))
  expect_true("gene" %in% colnames(result))

  # Should have results from multiple clusters
  expect_gt(length(unique(result$cluster)), 1)
})

test_that("find_all_ecm_markers requires multiple clusters", {
  skip_if_not_installed("Seurat")

  obj <- create_mock_seurat(n_clusters = 1)

  expect_error(
    find_all_ecm_markers(obj, group.by = "seurat_clusters", verbose = FALSE),
    "at least 2 clusters"
  )
})

# Test find_ecm_de_samples
test_that("find_ecm_de_samples works for condition comparison", {
  skip_if_not_installed("Seurat")

  obj <- create_mock_seurat()

  result <- find_ecm_de_samples(
    obj,
    condition.col = "condition",
    condition.1 = "Treatment",
    condition.2 = "Control",
    verbose = FALSE
  )

  expect_s3_class(result, "data.frame")
  expect_true("gene" %in% colnames(result))
  expect_true("avg_log2FC" %in% colnames(result))
})

test_that("find_ecm_de_samples works with cluster stratification", {
  skip_if_not_installed("Seurat")

  obj <- create_mock_seurat(n_cells = 200)

  result <- find_ecm_de_samples(
    obj,
    condition.col = "condition",
    condition.1 = "Treatment",
    condition.2 = "Control",
    cluster.col = "seurat_clusters",
    verbose = FALSE
  )

  expect_s3_class(result, "data.frame")

  # Should have cluster column if results exist
  if (nrow(result) > 0) {
    expect_true("cluster" %in% colnames(result))
  }
})

test_that("find_ecm_de_samples validates condition column", {
  skip_if_not_installed("Seurat")

  obj <- create_mock_seurat()

  expect_error(
    find_ecm_de_samples(
      obj,
      condition.col = "nonexistent",
      condition.1 = "A",
      condition.2 = "B"
    ),
    "not found in metadata"
  )
})

test_that("find_ecm_de_samples validates condition values", {
  skip_if_not_installed("Seurat")

  obj <- create_mock_seurat()

  expect_error(
    find_ecm_de_samples(
      obj,
      condition.col = "condition",
      condition.1 = "NonexistentCondition",
      condition.2 = "Control"
    ),
    "not found"
  )
})

test_that("find_ecm_de_samples pseudobulk requires sample.col", {
  skip_if_not_installed("Seurat")

  obj <- create_mock_seurat()

  expect_error(
    find_ecm_de_samples(
      obj,
      condition.col = "condition",
      condition.1 = "Treatment",
      condition.2 = "Control",
      pseudobulk = TRUE
    ),
    "sample.col is required"
  )
})

test_that("find_ecm_de_samples pseudobulk mode runs", {
  skip_if_not_installed("Seurat")

  # Need more cells for pseudobulk
  obj <- create_mock_seurat(n_cells = 200)

  result <- find_ecm_de_samples(
    obj,
    condition.col = "condition",
    condition.1 = "Treatment",
    condition.2 = "Control",
    pseudobulk = TRUE,
    sample.col = "sample_id",
    verbose = FALSE
  )

  expect_s3_class(result, "data.frame")
})

# Integration test
test_that("DE workflow integrates with other matrispace functions", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("UCell")

  obj <- create_mock_seurat()

  # This should work without errors
  markers <- find_ecm_markers(
    obj,
    ident.1 = "Cluster1",
    group.by = "seurat_clusters",
    verbose = FALSE
  )

  all_markers <- find_all_ecm_markers(
    obj,
    group.by = "seurat_clusters",
    verbose = FALSE
  )

  expect_s3_class(markers, "data.frame")
  expect_s3_class(all_markers, "data.frame")
})
