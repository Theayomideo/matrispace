#' @title ECM Niche Annotation Functions
#' @name ecm-annotation
#' @description Functions for classifying spots into ECM niches using the
#'   ScType algorithm.
#'
#' ScType functions adapted from:
#' Ianevski et al., 2022 (https://github.com/IanevskiAleksandr/sc-type)
#' Licensed under GPL-3.0
NULL

#' Annotate ECM niches using ScType
#'
#' Classifies each spot into an ECM niche (e.g. Interstitial, Basement
#' membrane) using the ScType algorithm with ECM-specific marker genes.
#'
#' @param seurat_obj A Seurat object with SCT assay
#' @param markers_path Path to Excel file with ECM niche markers.
#'   If NULL, uses built-in markers from inst/extdata.
#'
#' @return Seurat object with `ecm_niche` column in metadata
#'   (and `ecm_domain_annotation` retained as an alias for backward
#'   compatibility).
#'
#' @details
#' The function requires scaled SCT data. If the object lacks an SCT assay,
#' run \code{Seurat::SCTransform()} first.
#'
#' @examples
#' \dontrun{
#' seurat_obj <- annotate_ecm_niches(seurat_obj)
#' table(seurat_obj$ecm_niche)
#' }
#'
#' @export
annotate_ecm_niches <- function(seurat_obj, markers_path = NULL) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  # Find markers file
  if (is.null(markers_path)) {
    markers_path <- system.file(
      "extdata", "ECM_domains_transformed4ScType.xlsx",
      package = "matrispace"
    )
    if (markers_path == "") {
      warning("ECM niche marker database not found. Skipping annotation.")
      return(seurat_obj)
    }
  }

  if (!file.exists(markers_path)) {
    warning("ECM niche marker database not found at: ", markers_path)
    return(seurat_obj)
  }

  # Check for required dependencies
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required for ECM annotation")
  }

  # Check for SCT scale.data
  if (!"SCT" %in% SeuratObject::Assays(seurat_obj)) {
    stop("SCT assay not found. Run Seurat::SCTransform() first.")
  }

  # Prepare gene sets
  gs_list <- .gene_sets_prepare(markers_path, "ECM")

  # Get scale data
  scale_data <- .safe_get_assay_data(seurat_obj, assay = "SCT", slot = "scale.data")
  if (is.null(scale_data) || nrow(scale_data) == 0) {
    stop("No scale.data found in SCT assay. Re-run SCTransform.")
  }

  # Calculate scores
  es.max <- .sctype_score(
    scRNAseqData = scale_data,
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )

  # Get top score for each spot
  spot_results <- do.call("rbind", lapply(rownames(seurat_obj@meta.data), function(coord) {
    es.max.subset <- es.max[, coord, drop = FALSE]
    es.max.coord <- sort(rowSums(as.matrix(es.max.subset)), decreasing = TRUE)
    utils::head(data.frame(
      spot = coord,
      type = names(es.max.coord),
      scores = es.max.coord
    ), 1)
  }))

  spot_results$type[as.numeric(as.character(spot_results$scores)) <= 0] <- "not.assigned"

  # Add to metadata
  niche_calls <- spot_results$type[
    match(rownames(seurat_obj@meta.data), spot_results$spot)
  ]
  seurat_obj$ecm_niche <- niche_calls
  # Alias retained for backward compatibility with downstream code that
  # still reads `ecm_domain_annotation`. New code should use `ecm_niche`.
  seurat_obj$ecm_domain_annotation <- niche_calls

  return(seurat_obj)
}

#' Annotate ECM domains (deprecated)
#'
#' Deprecated alias for [annotate_ecm_niches()]. Use the new name in new code.
#'
#' @inheritParams annotate_ecm_niches
#' @return Same as [annotate_ecm_niches()].
#' @export
annotate_ecm_domains <- function(seurat_obj, markers_path = NULL) {
  .Deprecated("annotate_ecm_niches")
  annotate_ecm_niches(seurat_obj, markers_path = markers_path)
}

#' Prepare gene sets from ScType database file
#'
#' @param path_to_db_file Path to Excel database file
#' @param cell_type Cell/tissue type to filter
#' @return List with gs_positive and gs_negative gene sets
#' @keywords internal
#' @noRd
.gene_sets_prepare <- function(path_to_db_file, cell_type) {
  cell_markers <- openxlsx::read.xlsx(path_to_db_file)
  cell_markers <- cell_markers[cell_markers$tissueType == cell_type, ]
  cell_markers$geneSymbolmore1 <- gsub(" ", "", cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 <- gsub(" ", "", cell_markers$geneSymbolmore2)

  # Correct gene symbols (positive markers)
  cell_markers$geneSymbolmore1 <- sapply(seq_len(nrow(cell_markers)), function(i) {
    markers_all <- gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i], ",")))
    markers_all <- toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all <- sort(markers_all)

    if (length(markers_all) > 0) {
      suppressMessages({
        markers_all <- unique(stats::na.omit(
          HGNChelper::checkGeneSymbols(markers_all)$Suggested.Symbol
        ))
      })
      paste0(markers_all, collapse = ",")
    } else {
      ""
    }
  })

  # Correct gene symbols (negative markers)
  cell_markers$geneSymbolmore2 <- sapply(seq_len(nrow(cell_markers)), function(i) {
    markers_all <- gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i], ",")))
    markers_all <- toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all <- sort(markers_all)

    if (length(markers_all) > 0) {
      suppressMessages({
        markers_all <- unique(stats::na.omit(
          HGNChelper::checkGeneSymbols(markers_all)$Suggested.Symbol
        ))
      })
      paste0(markers_all, collapse = ",")
    } else {
      ""
    }
  })

  cell_markers$geneSymbolmore1 <- gsub("///", ",", cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 <- gsub("///", ",", cell_markers$geneSymbolmore2)

  gs <- lapply(seq_len(nrow(cell_markers)), function(j) {
    gsub(" ", "", unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]), ",")))
  })
  names(gs) <- cell_markers$cellName

  gs2 <- lapply(seq_len(nrow(cell_markers)), function(j) {
    gsub(" ", "", unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]), ",")))
  })
  names(gs2) <- cell_markers$cellName

  list(gs_positive = gs, gs_negative = gs2)
}

#' Calculate ScType scores
#'
#' @param scRNAseqData Expression matrix (genes x cells)
#' @param scaled Whether matrix is already scaled
#' @param gs List of positive marker gene sets
#' @param gs2 List of negative marker gene sets
#' @param gene_names_to_uppercase Convert gene names to uppercase
#' @return Matrix of ScType scores (cell types x cells)
#' @keywords internal
#' @noRd
.sctype_score <- function(scRNAseqData, scaled = TRUE, gs, gs2 = NULL,
                          gene_names_to_uppercase = TRUE) {

  # Marker sensitivity
  marker_stat <- sort(table(unlist(gs)), decreasing = TRUE)
  marker_sensitivity <- data.frame(
    score_marker_sensitivity = scales::rescale(
      as.numeric(marker_stat),
      to = c(0, 1),
      from = c(length(gs), 1)
    ),
    gene_ = names(marker_stat),
    stringsAsFactors = FALSE
  )

  # Convert gene names to uppercase
  if (gene_names_to_uppercase) {
    rownames(scRNAseqData) <- toupper(rownames(scRNAseqData))
  }

  # Subselect genes only found in data
  names_gs_cp <- names(gs)
  names_gs_2_cp <- names(gs2)

  gs <- lapply(seq_along(gs), function(d_) {
    gene_ind <- rownames(scRNAseqData) %in% as.character(gs[[d_]])
    rownames(scRNAseqData)[gene_ind]
  })

  gs2 <- lapply(seq_along(gs2), function(d_) {
    gene_ind <- rownames(scRNAseqData) %in% as.character(gs2[[d_]])
    rownames(scRNAseqData)[gene_ind]
  })

  names(gs) <- names_gs_cp
  names(gs2) <- names_gs_2_cp

  cell_markers_genes_score <- marker_sensitivity[
    marker_sensitivity$gene_ %in% unique(unlist(gs)),
  ]

  # Z-scale if not already
  if (!scaled) {
    Z <- t(scale(t(scRNAseqData)))
  } else {
    Z <- scRNAseqData
  }

  # Weight by marker sensitivity
  for (jj in seq_len(nrow(cell_markers_genes_score))) {
    gene <- cell_markers_genes_score[jj, "gene_"]
    sensitivity <- cell_markers_genes_score[jj, "score_marker_sensitivity"]
    Z[gene, ] <- Z[gene, ] * sensitivity
  }

  # Subselect only marker genes
  Z <- Z[unique(c(unlist(gs), unlist(gs2))), ]

  # Combine scores
  es <- do.call("rbind", lapply(names(gs), function(gss_) {
    sapply(seq_len(ncol(Z)), function(j) {
      gs_z <- Z[gs[[gss_]], j]
      gz_2 <- Z[gs2[[gss_]], j] * -1
      sum_t1 <- sum(gs_z) / sqrt(length(gs_z))
      sum_t2 <- sum(gz_2) / sqrt(length(gz_2))
      if (is.na(sum_t2)) sum_t2 <- 0
      sum_t1 + sum_t2
    })
  }))

  dimnames(es) <- list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all), ]

  es.max
}
