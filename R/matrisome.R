#' @title Matrisome Scoring Functions
#' @name matrisome
#' @description Functions for calculating matrisome category and subcategory
#'   expression scores using UCell.
NULL

#' Calculate matrisome feature scores
#'
#' Adds UCell-based module scores for the current MatriSpace matrisome
#' signatures to a Seurat object. Signatures include main categories
#' (Collagens, Glycoproteins, etc.) and curated subcategories such as
#' Basement Membrane, Laminins, and Peri-vascular ECM.
#'
#' @param seurat_obj A Seurat object
#' @param signature_list Named list of gene vectors. Defaults to package data.
#'
#' @return Seurat object with matrisome scores in metadata
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Translates gene symbols using HGNChelper
#'   \item Calculates UCell scores for each signature
#'   \item Creates display-named columns for UI compatibility
#' }
#'
#' @examples
#' \dontrun{
#' seurat_obj <- score_matrisome(seurat_obj)
#' }
#'
#' @export
score_matrisome <- function(seurat_obj, signature_list = NULL) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  if (is.null(signature_list)) {
    signature_list <- matrisome_signatures
  }

  translated_sigs <- translate_signatures(seurat_obj, signature_list)
  if (length(translated_sigs) == 0) {
    warning("No signature genes matched the Seurat object")
    return(seurat_obj)
  }

  # UCell scoring
  seurat_obj <- UCell::AddModuleScore_UCell(
    seurat_obj,
    features = translated_sigs,
    name = "_score"
  )

  # Clean up column names
  new_colnames <- gsub("_score1$", "", colnames(seurat_obj@meta.data))
  colnames(seurat_obj@meta.data) <- new_colnames

  # Map internal names to display names
  display_name_map <- c(
    "ecm_glycoproteins_score" = "ECM Glycoproteins",
    "collagens_score" = "Collagens",
    "proteoglycans_score" = "Proteoglycans",
    "ecm_regulators_score" = "ECM Regulators",
    "secreted_factors_score" = "Secreted Factors",
    "ecm-affiliated_proteins_score" = "ECM-affiliated Proteins",
    "basement_membrane_score" = "Basement Membrane",
    "hemostasis_score" = "Hemostasis",
    "elastic_fibers_score" = "Elastic fibers",
    "growth_factor-binding_score" = "Growth Factor-binding",
    "laminin_-_basement_membrane_score" = "Laminin",
    "matricellular_score" = "Matricellular proteins",
    "syndecan_score" = "Syndecan",
    "glypican_score" = "Glypican",
    "annexin_score" = "Annexins",
    "cathepsin_score" = "Cathepsins",
    "ccn_family_score" = "CCNs",
    "cystatin_score" = "Cystatins",
    "facit_score" = "FACITs",
    "fibulin_score" = "Fibulins",
    "galectin_score" = "Galectins",
    "mucin_score" = "Mucins",
    "plexin_score" = "Plexins",
    "semaphorin_score" = "Semaphorins"
  )

  # Copy UCell columns to display-named columns
  for (internal_name in names(display_name_map)) {
    if (internal_name %in% colnames(seurat_obj@meta.data)) {
      display_name <- display_name_map[[internal_name]]
      seurat_obj@meta.data[[display_name]] <- seurat_obj@meta.data[[internal_name]]
    }
  }

  return(seurat_obj)
}

#' Compute expression metrics for matrisome categories
#'
#' Calculates raw counts, log-scaled, and robust-normalized expression
#' values for a set of genes within each spot.
#'
#' @param gene_list Character vector of gene names
#' @param seurat_obj A Seurat object
#' @param category_name Name for the category (used in column names)
#' @param counts_matrix Pre-fetched counts matrix (optional)
#' @param matrisome_db Matrisome gene database with synonym info
#'
#' @return Data frame with expression metrics and coordinates
#'
#' @examples
#' \dontrun{
#' metrics <- compute_expression_metrics(
#'   gene_list = c("COL1A1", "COL1A2"),
#'   seurat_obj = my_obj,
#'   category_name = "Collagens"
#' )
#' }
#'
#' @export
compute_expression_metrics <- function(gene_list, seurat_obj, category_name,
                                       counts_matrix = NULL, matrisome_db = NULL) {
  .assert_seurat_object(seurat_obj, "seurat_obj")

  if (is.null(counts_matrix)) {
    counts_matrix <- .safe_get_assay_data(seurat_obj, slot = "counts")
  }
  if (is.null(matrisome_db)) {
    matrisome_db <- matrisome_db
  }

  available_genes <- c()

  # Check primary names first, then synonyms
  for (gene in gene_list) {
    if (gene %in% rownames(counts_matrix)) {
      available_genes <- c(available_genes, gene)
    } else {
      gene_row <- matrisome_db[matrisome_db$gene == gene, ]
      if (nrow(gene_row) == 1 && !is.na(gene_row$Synonyms)) {
        synonyms <- strsplit(gene_row$Synonyms, "\\||;|,|\\s+|/|-")[[1]]
        matching_synonym <- synonyms[synonyms %in% rownames(counts_matrix)]
        if (length(matching_synonym) > 0) {
          available_genes <- c(available_genes, matching_synonym[1])
        }
      }
    }
  }

  if (length(available_genes) == 0) {
    warning(paste("No matching genes found for", category_name))
    return(NULL)
  }

  # Get spot coordinates
  coords <- .safe_get_tissue_coordinates(seurat_obj)

  # Calculate metrics
  base_counts <- Matrix::colSums(counts_matrix[available_genes, , drop = FALSE])

  sums_log <- log1p(base_counts)
  log_range <- max(sums_log) - min(sums_log)
  log_scaled <- if (isTRUE(all.equal(log_range, 0))) {
    rep(0, length(base_counts))
  } else {
    (sums_log - min(sums_log)) / log_range
  }

  robust_spread <- stats::IQR(base_counts)
  robust <- if (!is.finite(robust_spread) || robust_spread == 0) {
    rep(0, length(base_counts))
  } else {
    (base_counts - stats::median(base_counts)) / robust_spread
  }

  # Build data frame
  df <- data.frame(
    spot = names(base_counts),
    base_counts = base_counts,
    log_scaled = log_scaled,
    robust = robust,
    imagerow = coords$imagerow,
    imagecol = coords$imagecol
  )

  return(df)
}
