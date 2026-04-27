# Tests for spatial analysis functions

test_that("compute_lisa runs on bundled Visium object", {
  skip_if_not_installed("Seurat")

  test_path <- system.file("testdata", "test_seurat.rds", package = "matrispace")
  skip_if(test_path == "", "bundled test_seurat.rds not available")

  obj <- readRDS(test_path)
  result <- compute_lisa(obj, "COL1A1", "gene", "FN1", "gene")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), ncol(obj))
  expect_true(all(c("feature1", "feature2", "LISA") %in% colnames(result)))
  expect_true(all(result$LISA %in%
                    c("High-High", "Low-Low", "High-Low", "Low-High",
                      "not.significant", "not.applicable")))
})

test_that("compute_morans_i returns NULL without MERINGUE", {
  skip_if_not_installed("Seurat")
  skip_if(requireNamespace("MERINGUE", quietly = TRUE),
          "MERINGUE installed; this test covers the missing-dependency path")

  counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- paste0("Cell", 1:10)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  expect_warning(result <- compute_morans_i(obj), "MERINGUE")
  expect_null(result)
})
