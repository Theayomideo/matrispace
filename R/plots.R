#' @title Visualization Functions
#' @name plots
#' @description Static plotting functions using ggplot2 and Seurat for
#'   matrisome spatial transcriptomics analysis.
#'
#' @import ggplot2
NULL

# Global variable declarations for ggplot2 aes() - avoids R CMD check NOTE
utils::globalVariables(c(
  "imagecol", "imagerow", "value", "robust", "log_scaled",
  "LISA_dynamic", "proportion", "id", "f1", "f2", "blend_color", "color"
))

#' Create spatial feature plot
#'
#' Generates a spatial feature plot using ggplot2 with RdYlBu color palette.
#'
#' @param seurat_obj A Seurat object
#' @param feature Feature name (gene or metadata column)
#' @param annotation_col Annotation column for hover text (optional)
#' @param pt.size Point size (default: 1.5)
#' @param title Plot title (defaults to feature name)
#'
#' @return A ggplot object
#'
#' @examples
#' \dontrun{
#' p <- plot_spatial_feature(seurat_obj, "COL1A1")
#' }
#'
#' @export
plot_spatial_feature <- function(seurat_obj, feature,
                                  annotation_col = NULL,
                                  pt.size = 1.5,
                                  title = NULL) {

  if (is.null(title)) title <- feature

  # Try Seurat's SpatialFeaturePlot first
  tryCatch({
    p <- Seurat::SpatialFeaturePlot(
      seurat_obj,
      features = feature,
      pt.size.factor = pt.size
    ) +
      scale_fill_gradientn(
        colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
        name = title
      ) +
      labs(title = title) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "right"
      )
    return(p)
  }, error = function(e) {
    # Fall back to manual plotting
    .manual_spatial_plot(seurat_obj, feature, pt.size, title)
  })
}

#' Manual spatial plot fallback
#' @keywords internal
#' @noRd
.manual_spatial_plot <- function(seurat_obj, feature, pt.size, title) {
  coords <- .safe_get_tissue_coordinates(seurat_obj)

  # Get feature values
  if (feature %in% colnames(seurat_obj@meta.data)) {
    values <- seurat_obj@meta.data[[feature]]
  } else {
    expr_data <- .safe_get_assay_data(seurat_obj)
    if (feature %in% rownames(expr_data)) {
      values <- expr_data[feature, ]
    } else {
      stop(paste("Feature not found:", feature))
    }
  }

  plot_df <- data.frame(
    imagecol = coords$imagecol,
    imagerow = coords$imagerow,
    value = values
  )

  ggplot(plot_df, aes(x = imagecol, y = imagerow, color = value)) +
    geom_point(size = pt.size) +
    scale_color_gradientn(
      colors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
      name = title
    ) +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    labs(title = title) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
}

#' Create matrisome category spatial plots
#'
#' Generates spatial plots for a matrisome category with robust and
#' log-scaled scores.
#'
#' @param seurat_obj A Seurat object with matrisome scores
#' @param category_name Display name (e.g., "Collagens")
#' @param internal_name Internal column prefix
#' @param pt.size Point size
#'
#' @return A patchwork object combining both plots, or NULL if data unavailable
#'
#' @examples
#' \dontrun{
#' p <- plot_matrisome(seurat_obj, "Collagens", "collagens")
#' }
#'
#' @export
plot_matrisome <- function(seurat_obj, category_name, internal_name, pt.size = 1.5) {

  robust_col <- paste0(internal_name, "_robust_score")
  log_col <- paste0(internal_name, "_log_scaled_score")

  # Validate columns exist
  if (!robust_col %in% colnames(seurat_obj@meta.data) ||
      !log_col %in% colnames(seurat_obj@meta.data)) {
    return(NULL)
  }

  # Check for sufficient data
  feature_values <- seurat_obj@meta.data[[robust_col]]
  non_na_count <- sum(!is.na(feature_values))
  if (non_na_count < 20) {
    return(NULL)
  }

  # Get coordinates
  coords <- .safe_get_tissue_coordinates(seurat_obj)

  # Build data frame
  df <- data.frame(
    imagecol = coords$imagecol,
    imagerow = coords$imagerow,
    robust = seurat_obj@meta.data[[robust_col]],
    log_scaled = seurat_obj@meta.data[[log_col]]
  )

  # Seurat spectral palette
  colors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))

  # Robust score plot
  p_robust <- ggplot(df, aes(x = imagecol, y = imagerow, color = robust)) +
    geom_point(size = pt.size) +
    scale_color_gradientn(colors = colors, name = "Robust\nScore") +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    labs(title = paste(category_name, "- Robust Score")) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))

  # Log-scaled plot
  p_log <- ggplot(df, aes(x = imagecol, y = imagerow, color = log_scaled)) +
    geom_point(size = pt.size) +
    scale_color_gradientn(colors = colors, name = "Log-scaled") +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    labs(title = paste(category_name, "- Log-scaled")) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))

  # Combine with patchwork
  patchwork::wrap_plots(p_robust, p_log, ncol = 2)
}

#' Create LISA clustering plot
#'
#' Creates a split-view visualization showing spatial LISA clusters
#' and quantitative summary bars.
#'
#' @param lisa_data Data frame from \code{compute_lisa}
#' @param selected_types Annotation types to include (NULL for all)
#' @param ann_col Column name for annotations
#' @param f1_name Feature 1 display name
#' @param f2_name Feature 2 display name
#'
#' @return A patchwork plot
#'
#' @examples
#' \dontrun{
#' lisa_data <- compute_lisa(seurat_obj, "COL1A1", "gene", "COL1A2", "gene")
#' p <- plot_lisa(lisa_data, f1_name = "COL1A1", f2_name = "COL1A2")
#' }
#'
#' @export
plot_lisa <- function(lisa_data, selected_types = NULL,
                      ann_col = "id", f1_name = "F1", f2_name = "F2") {

  k <- lisa_data
  k$id <- k[[ann_col]]

  # Handle missing feature2
  if (all(k$LISA == "not.applicable")) {
    return(
      ggplot(k, aes(x = imagecol, y = imagerow)) +
        geom_point(color = "grey80", size = 1.5) +
        theme_void() +
        scale_y_reverse() +
        coord_fixed() +
        labs(title = "LISA Clustering",
             subtitle = "Secondary feature not selected")
    )
  }

  # Create dynamic labels
  k$LISA_dynamic <- k$LISA
  k$LISA_dynamic <- gsub("High\\.High", paste0("High ", f1_name, " | High ", f2_name, " Neighbors"), k$LISA_dynamic)
  k$LISA_dynamic <- gsub("High\\.Low", paste0("High ", f1_name, " | Low ", f2_name, " Neighbors"), k$LISA_dynamic)
  k$LISA_dynamic <- gsub("Low\\.High", paste0("Low ", f1_name, " | High ", f2_name, " Neighbors"), k$LISA_dynamic)
  k$LISA_dynamic <- gsub("Low\\.Low", paste0("Low ", f1_name, " | Low ", f2_name, " Neighbors"), k$LISA_dynamic)

  dynamic_levels <- c(
    paste0("High ", f1_name, " | High ", f2_name, " Neighbors"),
    paste0("High ", f1_name, " | Low ", f2_name, " Neighbors"),
    paste0("Low ", f1_name, " | High ", f2_name, " Neighbors"),
    paste0("Low ", f1_name, " | Low ", f2_name, " Neighbors")
  )
  k$LISA_dynamic <- factor(k$LISA_dynamic, levels = dynamic_levels)

  lisa_colors <- c("#d62728", "#ff7f0e", "#9467bd", "#1f77b4")
  names(lisa_colors) <- dynamic_levels

  # Filter by selected cell types
  if (!is.null(selected_types)) {
    k <- k[k$id %in% selected_types, ]
  }
  if (nrow(k) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data for selected annotations"))
  }

  # Spatial plot
  p_spatial <- ggplot(k, aes(x = imagecol, y = imagerow, color = LISA_dynamic)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = lisa_colors, name = "Spatial Association Pattern", drop = FALSE) +
    theme_void() +
    scale_y_reverse() +
    coord_fixed() +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    theme(legend.position = "bottom")

  # Summary bar chart
  summary_data <- k %>%
    dplyr::count(.data$LISA_dynamic, .data$id, .drop = FALSE) %>%
    dplyr::group_by(.data$LISA_dynamic) %>%
    dplyr::mutate(proportion = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(id = .reorder_within(.data$id, .data$proportion, .data$LISA_dynamic))

  p_summary <- ggplot(summary_data, aes(x = proportion, y = id, fill = LISA_dynamic)) +
    geom_col() +
    scale_fill_manual(values = lisa_colors, drop = FALSE) +
    .scale_y_reordered() +
    scale_x_continuous(labels = scales::percent_format(), expand = c(0, 0.01)) +
    facet_wrap(~LISA_dynamic, ncol = 1, scales = "free_y") +
    theme_light(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(size = 9),
      axis.title.y = element_blank()
    ) +
    labs(x = "Proportion")

  # Combine
  patchwork::wrap_plots(p_spatial, p_summary, widths = c(2, 1))
}

#' Create blended spatial plot for two features
#'
#' Displays two features on the same spatial plot using a 2D color gradient.
#'
#' @param seurat_obj A Seurat object
#' @param features Character vector of exactly 2 feature names
#' @param feature1_name Display name for feature 1
#' @param feature2_name Display name for feature 2
#' @param bottom_left Color for low/low (default: black)
#' @param bottom_right Color for high feature1/low feature2 (default: red)
#' @param top_left Color for low feature1/high feature2 (default: green)
#' @param top_right Color for high/high (default: yellow)
#'
#' @return A patchwork plot with individual plots and blended result
#'
#' @examples
#' \dontrun{
#' p <- plot_spatial_blend(seurat_obj, c("COL1A1", "FN1"))
#' }
#'
#' @export
plot_spatial_blend <- function(seurat_obj, features,
                               feature1_name = NULL,
                               feature2_name = NULL,
                               bottom_left = "#000000",
                               bottom_right = "#FF0000",
                               top_left = "#00FF00",
                               top_right = "#FFFF00") {

  if (length(features) != 2) {
    stop("Requires exactly two features")
  }

  if (is.null(feature1_name)) feature1_name <- features[1]
  if (is.null(feature2_name)) feature2_name <- features[2]

  # Get data
  coords <- Seurat::GetTissueCoordinates(seurat_obj)
  dat <- Seurat::FetchData(seurat_obj, vars = features)

  df <- data.frame(
    imagecol = if ("imagecol" %in% names(coords)) coords$imagecol else coords$x,
    imagerow = if ("imagerow" %in% names(coords)) coords$imagerow else coords$y,
    f1 = dat[[features[1]]],
    f2 = dat[[features[2]]]
  )

  # Generate color grid
  side_length <- 100
  col_grid <- .gen_color_grid(side_length, bottom_left, bottom_right, top_left, top_right)

  # Normalize and map to colors
  df$f1_norm <- round((side_length - 1) * df$f1 / max(df$f1, na.rm = TRUE)) + 1
  df$f2_norm <- round((side_length - 1) * df$f2 / max(df$f2, na.rm = TRUE)) + 1

  df$blend_color <- mapply(function(x, y) {
    col_grid[x, y]
  }, df$f1_norm, df$f2_norm)

  # Individual plots
  p1 <- ggplot(df, aes(x = imagecol, y = imagerow, color = f1)) +
    geom_point(size = 1) +
    scale_color_gradient(low = bottom_left, high = bottom_right, name = feature1_name) +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    labs(title = feature1_name)

  p2 <- ggplot(df, aes(x = imagecol, y = imagerow, color = f2)) +
    geom_point(size = 1) +
    scale_color_gradient(low = bottom_left, high = top_left, name = feature2_name) +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    labs(title = feature2_name)

  # Blended plot
  p3 <- ggplot(df, aes(x = imagecol, y = imagerow, color = blend_color)) +
    geom_point(size = 1) +
    scale_color_identity() +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    labs(title = paste(feature1_name, "&", feature2_name))

  # Legend
  legend_grid <- expand.grid(
    f1 = seq(from = min(df$f1), to = max(df$f1), length.out = side_length),
    f2 = seq(from = min(df$f2), to = max(df$f2), length.out = side_length)
  )
  legend_grid$color <- c(col_grid)

  p_legend <- ggplot(legend_grid, aes(x = f1, y = f2, color = color)) +
    geom_point(shape = 15, size = 1.9) +
    scale_color_identity() +
    coord_cartesian(expand = FALSE) +
    theme(
      legend.position = "none",
      aspect.ratio = 1,
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    xlab(feature1_name) +
    ylab(feature2_name)

  patchwork::wrap_plots(p1, p2, p3, p_legend, nrow = 1, widths = c(1, 1, 1, 0.5))
}

#' Generate color grid for blend plots
#' @keywords internal
#' @noRd
.gen_color_grid <- function(side_length, bottom_left, bottom_right, top_left, top_right) {
  grad_gen <- function(start, end, n = side_length) {
    grDevices::colorRampPalette(c(start, end))(n)
  }

  bl_tl <- grad_gen(bottom_left, bottom_right)
  br_tr <- grad_gen(top_left, top_right)

  l <- lapply(seq_along(bl_tl), function(i) {
    grad_gen(bl_tl[i], br_tr[i])
  })

  t(matrix(unlist(l), ncol = side_length, nrow = side_length))
}

#' Add scale bar to spatial plot
#'
#' Adds a scale bar to a ggplot spatial plot based on Visium spot spacing.
#'
#' @param p A ggplot object
#' @param seurat_obj The Seurat object used to create the plot
#' @param bar_length_um Scale bar length in micrometers (default: 1000)
#'
#' @return The plot with scale bar added
#'
#' @examples
#' \dontrun{
#' p <- plot_spatial_feature(seurat_obj, "COL1A1")
#' p <- add_scalebar(p, seurat_obj)
#' }
#'
#' @export
add_scalebar <- function(p, seurat_obj, bar_length_um = 1000) {
  coords <- Seurat::GetTissueCoordinates(seurat_obj)
  if (is.null(coords) || nrow(coords) < 10) {
    return(p)
  }

  # Calculate pixel-to-micron ratio (100um between adjacent spots)
  sample_size <- min(100, nrow(coords))
  coord_cols <- if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
    c("imagecol", "imagerow")
  } else {
    c("x", "y")
  }
  dist_matrix <- as.matrix(stats::dist(coords[seq_len(sample_size), coord_cols]))
  diag(dist_matrix) <- NA
  pixels_per_100um <- min(dist_matrix, na.rm = TRUE)

  bar_length_px <- (bar_length_um / 100) * pixels_per_100um

  # Position in bottom-right, outside plot
  margin_x <- diff(range(coords[[coord_cols[1]]])) * 0.05
  x_end <- max(coords[[coord_cols[1]]]) - margin_x
  x_start <- x_end - bar_length_px

  y_offset <- diff(range(coords[[coord_cols[2]]])) * 0.05
  y_pos <- min(coords[[coord_cols[2]]]) - y_offset

  tick_height <- y_offset * 0.4
  bar_label <- paste(bar_length_um / 1000, "mm")

  # Allow drawing outside plot area
  p$coordinates$clip <- "off"

  p +
    annotate("segment", x = x_start, xend = x_end, y = y_pos, yend = y_pos,
             color = "black", linewidth = 1) +
    annotate("segment", x = x_start, xend = x_start, y = y_pos, yend = y_pos + tick_height,
             color = "black", linewidth = 1) +
    annotate("segment", x = x_end, xend = x_end, y = y_pos, yend = y_pos + tick_height,
             color = "black", linewidth = 1) +
    annotate("text", x = x_start + bar_length_px / 2, y = y_pos,
             label = bar_label, vjust = 2.0, color = "black", size = 4) +
    theme(plot.margin = margin(t = 5, r = 5, b = 30, l = 5, unit = "pt"))
}

#' Reorder within facets
#' @keywords internal
#' @noRd
.reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

#' Scale y reordered
#' @keywords internal
#' @noRd
.scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}
