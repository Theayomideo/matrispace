#' @title Preprocessing Functions
#' @name preprocessing
#' @description Functions for preparing Seurat and SpatialExperiment objects
#'   for matrisome analysis.
NULL

#' Convert SpatialExperiment to Seurat object
#'
#' Converts a SpatialExperiment object to a Seurat object with spatial image
#' data preserved. Handles multi-sample objects by selecting the first sample.
#'
#' @param spe A SpatialExperiment object
#' @param spatial_cols Named character vector mapping spatial columns:
#'   \itemize{
#'     \item tissue: Column for tissue/in_tissue flag
#'     \item row: Array row coordinate
#'     \item col: Array column coordinate
#'     \item imagerow: Pixel row coordinate (required)
#'     \item imagecol: Pixel column coordinate (required)
#'   }
#' @param symbol_col Column in rowData containing gene symbols. Auto-detected if NULL.
#' @param verbose Print progress messages
#'
#' @return A Seurat object with VisiumV1 image
#'
#' @examples
#' \dontrun{
#' seurat_obj <- convert_spe_to_seurat(spe)
#' }
#'
#' @export
convert_spe_to_seurat <- function(
    spe,
    spatial_cols = c(
      "tissue" = "in_tissue",
      "row" = "array_row",
      "col" = "array_col",
      "imagerow" = "pxl_row_in_fullres",
      "imagecol" = "pxl_col_in_fullres"
    ),
    symbol_col = NULL,
    verbose = TRUE
) {

  # Check dependencies
  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("Package 'SpatialExperiment' is required for SPE conversion")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required for SPE conversion")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required for SPE conversion")
  }

  # Validate input class
  if (!inherits(spe, "SpatialExperiment")) {
    stop("Input must be a SpatialExperiment object")
  }

  # Multi-sample handling: auto-subset first sample
  sample_ids <- unique(spe$sample_id)
  if (length(sample_ids) > 1) {
    if (verbose) message(sprintf(
      "Multi-sample SPE detected (%d samples). Using first: '%s'",
      length(sample_ids), sample_ids[1]
    ))
    spe <- spe[, spe$sample_id == sample_ids[1]]
  }
  sample_id <- unique(spe$sample_id)

  # Check counts assay exists
  if (!"counts" %in% SummarizedExperiment::assayNames(spe)) {
    stop("No 'counts' assay found. Available: ",
         paste(SummarizedExperiment::assayNames(spe), collapse = ", "))
  }

  # Check required spatial_cols elements
  required_elements <- c("imagerow", "imagecol")
  missing_elements <- setdiff(required_elements, names(spatial_cols))
  if (length(missing_elements) > 0) {
    stop("spatial_cols must contain: ", paste(missing_elements, collapse = ", "))
  }

  # Combine colData and spatialCoords
  col_info <- cbind(
    SummarizedExperiment::colData(spe),
    SpatialExperiment::spatialCoords(spe)
  )

  # Check imagerow/imagecol exist
  for (coord in c("imagerow", "imagecol")) {
    col_name <- spatial_cols[coord]
    if (is.na(col_name) || !col_name %in% colnames(col_info)) {
      stop(sprintf(
        "'%s' column '%s' not found. Available: %s",
        coord, col_name, paste(colnames(col_info), collapse = ", ")
      ))
    }
  }

  # Check imgData exists
  img_data <- SpatialExperiment::imgData(spe)
  if (nrow(img_data) == 0) {
    stop("No image data found in imgData(spe)")
  }
  if (!"lowres" %in% img_data$image_id) {
    stop("No 'lowres' image found. Available: ",
         paste(img_data$image_id, collapse = ", "))
  }

  # Auto-detect symbol column
  rd_cols <- colnames(SummarizedExperiment::rowData(spe))
  common_symbol_cols <- c("symbol", "gene_name", "Symbol", "gene_short_name", "SYMBOL")

  if (is.null(symbol_col)) {
    detected <- intersect(common_symbol_cols, rd_cols)
    if (length(detected) > 0) {
      symbol_col <- detected[1]
      if (verbose) message(sprintf("Auto-detected symbol column: '%s'", symbol_col))
    }
  } else if (!symbol_col %in% rd_cols) {
    warning(sprintf("Specified symbol_col '%s' not found in rowData", symbol_col))
    symbol_col <- NULL
  }

  has_symbols <- !is.null(symbol_col) &&
    !all(is.na(SummarizedExperiment::rowData(spe)[[symbol_col]]))

  if (verbose) {
    message(sprintf("Input: %d genes, %d spots", nrow(spe), ncol(spe)))
  }

  # Convert
  SPOT_DIAMETER <- 55e-6

  if (verbose) message("Converting to Seurat object...")

  # Remove altExps
  spe_clean <- spe
  if (length(SingleCellExperiment::altExpNames(spe_clean)) > 0) {
    for (ae in SingleCellExperiment::altExpNames(spe_clean)) {
      SingleCellExperiment::altExp(spe_clean, ae) <- NULL
    }
  }

  seur <- Seurat::as.Seurat(spe_clean, counts = "counts", data = NULL)
  assay_name <- Seurat::DefaultAssay(seur)

  # Update gene names to symbols if available
  if (has_symbols) {
    if (verbose) message(sprintf("Mapping gene symbols from '%s' column...", symbol_col))
    ensembl_ids <- rownames(seur)
    symbols <- SummarizedExperiment::rowData(spe)[[symbol_col]]
    new_names <- ifelse(is.na(symbols) | symbols == "", ensembl_ids, symbols)
    new_names <- make.unique(new_names)

    rownames(seur@assays[[assay_name]]@counts) <- new_names
    rownames(seur@assays[[assay_name]]@data) <- new_names
  }

  if (verbose) message(sprintf("Adding coordinates and image for sample %s...", sample_id))

  # Helper to get column with fallback
  get_col <- function(key, fallback) {
    col_name <- spatial_cols[key]
    if (!is.na(col_name) && col_name %in% colnames(col_info)) {
      return(col_info[[col_name]])
    }
    return(fallback)
  }

  # Build coordinates
  coords <- data.frame(
    tissue = as.integer(get_col("tissue", rep(1L, ncol(spe)))),
    row = get_col("row", seq_len(ncol(spe))),
    col = get_col("col", seq_len(ncol(spe))),
    imagerow = col_info[[spatial_cols["imagerow"]]],
    imagecol = col_info[[spatial_cols["imagecol"]]],
    row.names = colnames(spe)
  )

  # Convert image to array
  this_img <- array(
    t(grDevices::col2rgb(SpatialExperiment::imgRaster(spe))),
    dim = c(dim(SpatialExperiment::imgRaster(spe)), 3)
  ) / 256

  # Get scale factor
  sf <- img_data$scaleFactor[img_data$image_id == "lowres"]
  safe_key <- paste0(gsub("[^[:alnum:]]", "", sample_id), "_")

  # Create VisiumV1 object
  seur@images[[sample_id]] <- methods::new(
    Class = "VisiumV1",
    image = this_img,
    scale.factors = Seurat::scalefactors(
      spot = NA, fiducial = NA, hires = NA, lowres = sf
    ),
    coordinates = coords,
    spot.radius = SPOT_DIAMETER / sf,
    assay = assay_name,
    key = safe_key
  )

  if (verbose) message("Done!")
  return(seur)
}

#' Prepare Seurat object for matrisome analysis
#'
#' Master preprocessing function that checks and adds required components:
#' \enumerate{
#'   \item SCTransform (if needed for ScType)
#'   \item Matrisome feature scores (UCell)
#'   \item ECM niche classification (ScType)
#'   \item ECM niche signatures
#' }
#'
#' @param seurat_obj A Seurat object
#' @param mat_feat_sigs List of matrisome gene signatures (defaults to package data)
#' @param ecm_ucell_sigs List of ECM niche signatures (defaults to package data)
#' @param verbose Print progress messages
#'
#' @return Preprocessed Seurat object
#'
#' @examples
#' \dontrun{
#' seurat_obj <- prepare_object(seurat_obj)
#' }
#'
#' @export
prepare_object <- function(seurat_obj,
                           mat_feat_sigs = NULL,
                           ecm_ucell_sigs = NULL,
                           verbose = TRUE) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  # Load package data if not provided
  if (is.null(mat_feat_sigs)) {
    mat_feat_sigs <- matrisome_signatures
  }
  if (is.null(ecm_ucell_sigs)) {
    ecm_ucell_sigs <- ecm_signatures
  }

  has_niche <- any(c("ecm_niche", "ecm_domain_annotation") %in%
                     colnames(seurat_obj@meta.data))

  # Step 1: Check & Run SCTransform
  needs_sctype <- !has_niche
  has_scale_data <- .safe_has_data(seurat_obj, assay = "SCT", slot = "scale.data")

  if (needs_sctype && !has_scale_data) {
    .message_if(verbose, "SCT data not found. Running SCTransform...")
    seurat_obj <- Seurat::SCTransform(
      seurat_obj,
      assay = Seurat::DefaultAssay(seurat_obj),
      verbose = FALSE
    )
  }
  if ("SCT" %in% SeuratObject::Assays(seurat_obj)) {
    Seurat::DefaultAssay(seurat_obj) <- "SCT"
  }

  # Step 2: Add Matrisome Feature Signatures
  if (!"collagens" %in% colnames(seurat_obj@meta.data)) {
    .message_if(verbose, "Matrisome feature scores not found. Calculating now...")
    seurat_obj <- score_matrisome(seurat_obj, mat_feat_sigs)
  }

  # Step 3: Add ECM niche classification
  if (!has_niche) {
    .message_if(verbose, "ECM niche annotation not found. Running ScType...")
    seurat_obj <- annotate_ecm_niches(seurat_obj)
  }

  # Step 4: Add ECM Niche Signatures
  if (!"Interstitial_UCell" %in% colnames(seurat_obj@meta.data)) {
    .message_if(verbose, "ECM niche scores not found. Running UCell...")
    translated_sigs <- translate_signatures(seurat_obj, ecm_ucell_sigs)
    if (length(translated_sigs) > 0) {
      seurat_obj <- UCell::AddModuleScore_UCell(seurat_obj, features = translated_sigs)
    }
  }

  .message_if(verbose, "Preprocessing complete!")
  return(seurat_obj)
}

#' Translate gene signatures using HGNChelper
#'
#' Checks gene symbols and corrects outdated/alias symbols to match
#' what's present in the Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param signature_list Named list of gene vectors
#'
#' @return Named list of corrected gene vectors
#'
#' @examples
#' \dontrun{
#' sigs <- translate_signatures(seurat_obj, my_signatures)
#' }
#'
#' @export
translate_signatures <- function(seurat_obj, signature_list) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  if (length(signature_list) == 0) {
    return(list())
  }

  seurat_genes <- rownames(seurat_obj)
  translated_signatures <- list()

  for (sig_name in names(signature_list)) {
    original_genes <- signature_list[[sig_name]]
    if (length(original_genes) == 0) {
      next
    }

    # Check gene symbols and get corrections
    correction_map <- suppressMessages(
      HGNChelper::checkGeneSymbols(original_genes, unmapped.as.na = FALSE)
    )

    corrected_genes <- correction_map$Suggested.Symbol
    missing_suggestions <- is.na(corrected_genes) | corrected_genes == ""
    corrected_genes[missing_suggestions] <- original_genes[missing_suggestions]

    # Filter for genes actually present in the data
    valid_genes <- unique(corrected_genes[corrected_genes %in% seurat_genes])

    if (length(valid_genes) > 0) {
      translated_signatures[[sig_name]] <- valid_genes
    }
  }
  return(translated_signatures)
}

#' Combine spatial coordinates with metadata
#'
#' @param obj A Seurat object
#' @return Data frame with coordinates and metadata
#' @keywords internal
#' @noRd
.prepare_plot_data <- function(obj) {
  spot_coords <- .safe_get_tissue_coordinates(obj)
  metadata <- obj@meta.data
  cbind(spot_coords, metadata)
}

#' Clean Seurat metadata column names
#'
#' Renames legacy column naming conventions to standard display names.
#'
#' @param seurat_obj A Seurat object
#' @return Seurat object with cleaned metadata column names
#' @keywords internal
#' @noRd
.clean_metadata <- function(seurat_obj) {
  # Column renaming mappings
  column_mappings <- list(
    list(old = "glycoproteins.Spatial1", new = "glycoproteins_Spatial", display = "ECM Glycoproteins"),
    list(old = "collagens.Spatial1", new = "collagens_Spatial", display = "Collagens"),
    list(old = "proteoglycans.Spatial1", new = "proteoglycans_Spatial", display = "Proteoglycans"),
    list(old = "affiliated.Spatial1", new = "affiliated_Spatial", display = "ECM-affiliated Proteins"),
    list(old = "regulators.Spatial1", new = "regulators_Spatial", display = "ECM Regulators"),
    list(old = "secreted.Spatial1", new = "secreted_Spatial", display = "Secreted Factors"),
    list(old = "basement_mem.Spatial1", new = "basement_mem_Spatial", display = "Basement Membrane"),
    list(old = "hemostasis.Spatial1", new = "hemostasis_Spatial", display = "Hemostasis"),
    list(old = "elastic_fibers.Spatial1", new = "elastic_fibers_Spatial", display = "Elastic fibers"),
    list(old = "gf_binding.Spatial1", new = "gf_binding_Spatial", display = "Growth Factor-binding")
  )

  for (mapping in column_mappings) {
    if (mapping$old %in% colnames(seurat_obj@meta.data)) {
      colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == mapping$old] <- mapping$display
    } else if (mapping$new %in% colnames(seurat_obj@meta.data)) {
      colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == mapping$new] <- mapping$display
    }
  }

  return(seurat_obj)
}
