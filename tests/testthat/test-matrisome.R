# Tests for matrisome scoring functions

test_that("matrisome_db loads correctly", {
  data(matrisome_db)
  expect_s3_class(matrisome_db, "data.frame")
  expect_true("gene" %in% colnames(matrisome_db))
  expect_gt(nrow(matrisome_db), 1000)
})

test_that("matrisome_signatures loads correctly", {
  data(matrisome_signatures)
  expect_type(matrisome_signatures, "list")
  expect_gt(length(matrisome_signatures), 20)
  expect_true(all(sapply(matrisome_signatures, is.character)))
})

test_that("ecm_signatures loads correctly", {
  data(ecm_signatures)
  expect_type(ecm_signatures, "list")
  expect_true("Vascular" %in% names(ecm_signatures))
  expect_true("Basement" %in% names(ecm_signatures))
  expect_true("Interstitial" %in% names(ecm_signatures))
})

test_that("translate_signatures handles empty input", {
  skip_if_not_installed("Seurat")

  # Create minimal Seurat object
  counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- paste0("Cell", 1:10)
  obj <- Seurat::CreateSeuratObject(counts = counts)

  result <- translate_signatures(obj, list(test = c("NonexistentGene")))
  expect_length(result, 0)
})

test_that("palette_annotation returns valid colors", {
  colors <- palette_annotation(c("A", "B", "C", "not.assigned"))
  expect_length(colors, 4)
  expect_true(all(grepl("^#", colors)))
  expect_equal(unname(colors["not.assigned"]), "#7f7f7f")
})

test_that("palette_ecm returns correct domain colors", {
  colors <- palette_ecm(c("Vascular ECM", "Basement ECM", "Interstitial ECM", "not.assigned"))
  expect_length(colors, 4)
  expect_true(grepl("#d42626", colors["Vascular ECM"]))
  expect_true(grepl("#6a3d9a", colors["Basement ECM"]))
  expect_true(grepl("#2b9e2b", colors["Interstitial ECM"]))
})
