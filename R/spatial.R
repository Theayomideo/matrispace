#' @title Spatial Statistics Functions
#' @name spatial
#' @description Functions for computing spatial statistics including LISA
#'   (Local Indicators of Spatial Association) and Moran's I autocorrelation.
NULL

#' Compute LISA statistics for feature pairs
#'
#' Calculates Local Indicators of Spatial Association (LISA) to identify
#' spatial clustering patterns between two features.
#'
#' @param seurat_obj A Seurat object
#' @param feature1 First feature name (gene or metadata column)
#' @param feature1_type Type of feature1: "gene", "matrisome gene", "signature", or "any gene"
#' @param feature2 Second feature name (optional, for bivariate LISA)
#' @param feature2_type Type of feature2
#'
#' @return Data frame with coordinates, feature values, and LISA classifications
#'
#' @details
#' LISA classifications are based on:
#' \itemize{
#'   \item High.High: High feature1, high neighbor feature2
#'   \item High.Low: High feature1, low neighbor feature2
#'   \item Low.High: Low feature1, high neighbor feature2
#'   \item Low.Low: Low feature1, low neighbor feature2
#' }
#'
#' @examples
#' \dontrun{
#' lisa_results <- compute_lisa(seurat_obj, "COL1A1", "gene", "COL1A2", "gene")
#' }
#'
#' @export
compute_lisa <- function(seurat_obj, feature1, feature1_type,
                         feature2 = NULL, feature2_type = NULL) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  m <- seurat_obj@meta.data

  # Handle primary feature
  if (feature1 == "") {
    m$feature1 <- 0
  } else if (feature1_type %in% c("matrisome gene", "any gene", "gene")) {
    feature1_data <- .safe_get_assay_data(seurat_obj)
    if (!is.null(feature1_data) && feature1 %in% rownames(feature1_data)) {
      m$feature1 <- feature1_data[feature1, ]
    } else {
      warning("Could not retrieve data for feature1")
      m$feature1 <- 0
    }
  } else if (feature1_type %in% c("matrisome signature", "signature")) {
    if (feature1 %in% colnames(seurat_obj@meta.data)) {
      m$feature1 <- seurat_obj@meta.data[, feature1]
    } else {
      m$feature1 <- 0
    }
  }

  # Handle secondary feature
  if (is.null(feature2_type) || is.null(feature2) || feature2 == "" || feature2_type == "none") {
    m$feature2 <- 0
  } else if (feature2_type %in% c("matrisome gene", "any gene", "gene")) {
    feature2_data <- .safe_get_assay_data(seurat_obj)
    if (!is.null(feature2_data) && feature2 %in% rownames(feature2_data)) {
      m$feature2 <- feature2_data[feature2, ]
    } else {
      warning("Could not retrieve data for feature2")
      m$feature2 <- 0
    }
  } else if (feature2_type %in% c("matrisome signature", "signature")) {
    if (feature2 %in% colnames(seurat_obj@meta.data)) {
      m$feature2 <- seurat_obj@meta.data[, feature2]
    } else {
      m$feature2 <- 0
    }
  }

  # Get spatial coordinates
  coords <- .safe_get_tissue_coordinates(seurat_obj)
  m <- cbind(m, coords)

  # Clean up negative values
  m$feature1[m$feature1 < 0] <- 0
  m$feature2[m$feature2 < 0] <- 0

  # Calculate LISA if secondary feature has variance
  if (!is.null(feature2_type) && feature2_type != "none" && feature2 != "" &&
      stats::var(m$feature2, na.rm = TRUE) > 0) {

    # Calculate spatial weight matrix
    Wij <- as.matrix(stats::dist(as.matrix(cbind(m$row, m$col))))
    Wij[Wij == 0] <- 1e-9
    Wij <- 1 / Wij
    diag(Wij) <- 0
    if (sum(Wij) > 0) Wij <- Wij / sum(Wij)
    diag(Wij) <- 0

    # Standardize features (Z-score)
    x <- m$feature1
    y <- m$feature2

    x_scaled <- if (stats::sd(x, na.rm = TRUE) > 0) as.numeric(scale(x)) else rep(0, length(x))
    y_scaled <- if (stats::sd(y, na.rm = TRUE) > 0) as.numeric(scale(y)) else rep(0, length(y))
    x_scaled[is.na(x_scaled)] <- 0
    y_scaled[is.na(y_scaled)] <- 0

    # Calculate LISA clustering
    lisa.clust <- as.character(interaction(x_scaled > 0, Wij %*% y_scaled > 0))
    lisa.clust <- gsub("TRUE", "High", lisa.clust)
    lisa.clust <- gsub("FALSE", "Low", lisa.clust)
    m$LISA <- lisa.clust
  } else {
    m$LISA <- "not.applicable"
  }

  return(m)
}

#' Compute Moran's I spatial autocorrelation
#'
#' Calculates global Moran's I statistic for matrisome scores to assess
#' spatial clustering patterns.
#'
#' @param seurat_obj A Seurat object with matrisome scores
#' @param features Character vector of feature names to test. If NULL,
#'   automatically detects columns ending in "_robust_score".
#'
#' @return Data frame with Moran's I and p-values for each feature
#'
#' @details
#' Requires the MERINGUE package. If not installed, returns NULL with a warning.
#'
#' @examples
#' \dontrun{
#' autocorr <- compute_morans_i(seurat_obj)
#' }
#'
#' @export
compute_morans_i <- function(seurat_obj, features = NULL) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  if (!requireNamespace("MERINGUE", quietly = TRUE)) {
    warning("Package 'MERINGUE' is required for Moran's I calculation")
    return(NULL)
  }

  tryCatch({
    coords <- .safe_get_tissue_coordinates(seurat_obj)
    if (nrow(coords) < 3) return(NULL)

    weight_matrix <- MERINGUE::getSpatialNeighbors(coords)

    if (is.null(features)) {
      features <- grep("_robust_score", colnames(seurat_obj@meta.data), value = TRUE)
    }

    if (length(features) == 0) {
      warning("No features found for autocorrelation analysis")
      return(NULL)
    }

    results_list <- lapply(features, function(col_name) {
      scores <- seurat_obj@meta.data[[col_name]]
      scores <- scores[is.finite(scores)]

      if (length(scores) < 3 || stats::var(scores) == 0) {
        return(NULL)
      }

      names(scores) <- rownames(seurat_obj@meta.data)[
        is.finite(seurat_obj@meta.data[[col_name]])
      ]

      test_result <- MERINGUE::moranTest(scores, weight_matrix)

      data.frame(
        category = col_name,
        morans_I = test_result["observed"],
        p_value = test_result["p.value"],
        stringsAsFactors = FALSE
      )
    })

    results_df <- do.call(rbind, results_list[!sapply(results_list, is.null)])
    return(results_df)

  }, error = function(e) {
    warning(paste("Error in autocorrelation calculation:", e$message))
    return(NULL)
  })
}
