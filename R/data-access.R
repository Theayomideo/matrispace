#' @title Internal Seurat Access Helpers
#' @name data-access
#' @description Internal helpers for validating Seurat inputs and retrieving
#'   assays, feature matrices, and spatial coordinates across Seurat versions.
#' @keywords internal
NULL

#' Validate Seurat input
#'
#' @param object Object to validate
#' @param arg Argument name used in error messages
#' @keywords internal
#' @noRd
.assert_seurat_object <- function(object, arg = "object") {
  if (!inherits(object, "Seurat")) {
    stop(arg, " must be a Seurat object")
  }
}

#' Emit a message when verbose output is enabled
#'
#' @param verbose Logical flag
#' @param ... Message parts
#' @keywords internal
#' @noRd
.message_if <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    message(...)
  }
}

#' Resolve which assay to use
#'
#' @param object A Seurat object
#' @param assay Assay name. If NULL, uses SCT > Spatial > default assay > first available.
#' @return Character scalar naming the assay
#' @keywords internal
#' @noRd
.resolve_assay <- function(object, assay = NULL) {
  .assert_seurat_object(object)

  assays <- SeuratObject::Assays(object)
  if (length(assays) == 0) {
    stop("No assays found in Seurat object")
  }

  if (!is.null(assay)) {
    if (!assay %in% assays) {
      stop(
        "Assay '", assay, "' not found. Available assays: ",
        paste(assays, collapse = ", ")
      )
    }
    return(assay)
  }

  preferred <- c("SCT", "Spatial", Seurat::DefaultAssay(object), assays)
  preferred <- preferred[!is.na(preferred) & nzchar(preferred)]
  preferred[match(TRUE, preferred %in% assays)]
}

#' Safely access assay data across Seurat versions
#'
#' @param object A Seurat object
#' @param assay Assay name. If NULL, uses the preferred assay.
#' @param slot Data layer: "counts", "data", or "scale.data". If NULL, tries in order.
#' @return A matrix-like object or NULL if unavailable
#' @keywords internal
#' @noRd
.safe_get_assay_data <- function(object, assay = NULL, slot = NULL) {
  assay <- .resolve_assay(object, assay)
  layers_to_try <- if (is.null(slot)) c("data", "counts", "scale.data") else slot

  fetch_assay_data <- function(layer_name) {
    tryCatch(
      Seurat::GetAssayData(object = object, assay = assay, layer = layer_name),
      error = function(layer_error) {
        tryCatch(
          Seurat::GetAssayData(object = object, assay = assay, slot = layer_name),
          error = function(slot_error) {
            NULL
          }
        )
      }
    )
  }

  for (layer_name in layers_to_try) {
    data <- fetch_assay_data(layer_name)
    if (!is.null(data) && length(data) > 0) {
      return(data)
    }
  }

  NULL
}

#' Check whether assay data are present
#'
#' @param object A Seurat object
#' @param assay Assay name (optional)
#' @param slot Data layer to inspect
#' @return Logical scalar
#' @keywords internal
#' @noRd
.safe_has_data <- function(object, assay = NULL, slot = "data") {
  data <- .safe_get_assay_data(object, assay = assay, slot = slot)
  !is.null(data) && all(dim(data) > 0)
}

#' Get gene names from Seurat object
#'
#' @param object A Seurat object
#' @param assay Assay name (optional)
#' @param slot Data layer (optional)
#' @return Character vector of gene names
#' @keywords internal
#' @noRd
.safe_get_rownames <- function(object, assay = NULL, slot = NULL) {
  data <- .safe_get_assay_data(object, assay = assay, slot = slot)
  if (is.null(data)) {
    return(character(0))
  }
  rownames(data)
}

#' Get number of features in Seurat object
#'
#' @param object A Seurat object
#' @param assay Assay name (optional)
#' @return Integer count of features
#' @keywords internal
#' @noRd
.safe_get_nfeatures <- function(object, assay = NULL) {
  data <- .safe_get_assay_data(object, assay = assay)
  if (is.null(data)) {
    return(0L)
  }
  nrow(data)
}

#' Safely retrieve tissue coordinates
#'
#' @param object A Seurat object
#' @return Data frame with row, col, imagerow, and imagecol columns when possible
#' @keywords internal
#' @noRd
.safe_get_tissue_coordinates <- function(object) {
  .assert_seurat_object(object)

  coords <- tryCatch(
    Seurat::GetTissueCoordinates(object),
    error = function(e) NULL
  )

  if (is.null(coords) || nrow(coords) == 0) {
    stop("Could not retrieve spatial coordinates from the Seurat object")
  }

  if (!"row" %in% colnames(coords)) {
    coords$row <- if ("y" %in% colnames(coords)) coords$y else seq_len(nrow(coords))
  }
  if (!"col" %in% colnames(coords)) {
    coords$col <- if ("x" %in% colnames(coords)) coords$x else seq_len(nrow(coords))
  }
  if (!"imagerow" %in% colnames(coords)) {
    coords$imagerow <- if ("y" %in% colnames(coords)) coords$y else coords$row
  }
  if (!"imagecol" %in% colnames(coords)) {
    coords$imagecol <- if ("x" %in% colnames(coords)) coords$x else coords$col
  }

  coords
}
