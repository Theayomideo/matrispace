#' @title Color Palette Functions
#' @name color-palettes
#' @description Functions for generating consistent color palettes for
#'   matrisome visualizations.
NULL

#' Generate annotation color palette
#'
#' Creates a consistent color palette for categorical annotations, with
#' special handling for "not.assigned" values (always grey).
#'
#' @param levels_vector Character vector of annotation levels
#' @return Named character vector of hex colors
#'
#' @details
#' Uses the Glasbey palette for up to 32 distinct colors. For more categories,
#' falls back to viridis interpolation. "not.assigned" is always mapped to grey.
#'
#' @examples
#' \dontrun{
#' colors <- palette_annotation(c("Tumor", "Stroma", "not.assigned"))
#' }
#'
#' @export
palette_annotation <- function(levels_vector) {
  # Ensure input is unique characters and sorted for consistency
  all_levels <- sort(unique(as.character(levels_vector)))

  # Define the neutral color for unassigned spots

not_assigned_color <- "#7f7f7f"

  # Check if "not.assigned" exists
  if ("not.assigned" %in% all_levels) {
    real_levels <- all_levels[all_levels != "not.assigned"]
    num_real_levels <- length(real_levels)

    if (num_real_levels > 0) {
      if (num_real_levels > 32) {
        # Use viridis for many categories
        if (requireNamespace("viridis", quietly = TRUE)) {
          palette_func <- grDevices::colorRampPalette(viridis::viridis(num_real_levels))
          real_palette <- palette_func(num_real_levels)
        } else {
          real_palette <- grDevices::rainbow(num_real_levels)
        }
      } else {
        # Use glasbey palette if scCustomize available
        if (requireNamespace("scCustomize", quietly = TRUE)) {
          real_palette <- scCustomize::DiscretePalette_scCustomize(
            num_colors = num_real_levels,
            palette = "glasbey"
          )
        } else {
          real_palette <- scales::hue_pal()(num_real_levels)
        }
      }
      real_color_map <- stats::setNames(real_palette, real_levels)
    } else {
      real_color_map <- character(0)
    }

    not_assigned_map <- stats::setNames(not_assigned_color, "not.assigned")
    final_color_map <- c(real_color_map, not_assigned_map)

  } else {
    num_all_levels <- length(all_levels)
    if (num_all_levels > 0) {
      if (num_all_levels > 32) {
        if (requireNamespace("viridis", quietly = TRUE)) {
          palette_func <- grDevices::colorRampPalette(viridis::viridis(num_all_levels))
          palette <- palette_func(num_all_levels)
        } else {
          palette <- grDevices::rainbow(num_all_levels)
        }
      } else {
        if (requireNamespace("scCustomize", quietly = TRUE)) {
          palette <- scCustomize::DiscretePalette_scCustomize(
            num_colors = num_all_levels,
            palette = "glasbey"
          )
        } else {
          palette <- scales::hue_pal()(num_all_levels)
        }
      }
      final_color_map <- stats::setNames(palette, all_levels)
    } else {
      final_color_map <- character(0)
    }
  }

  return(final_color_map)
}

#' Generate ECM niche color palette
#'
#' Creates a color palette for ECM niche annotations with biologically
#' meaningful colors: Vascular (red), Basement (purple), Interstitial (green).
#'
#' @param levels_vector Character vector of ECM niche levels
#' @return Named character vector of hex colors
#'
#' @details
#' ECM Niche Color Palette:
#' \itemize{
#'   \item Vascular (#d42626): Red - blood vessels and vascular ECM
#'   \item Basement (#6a3d9a): Purple - structural foundation
#'   \item Interstitial (#2b9e2b): Green - matrix between cells
#'   \item not.assigned (#7f7f7f): Grey - unassigned spots
#' }
#'
#' @examples
#' \dontrun{
#' colors <- palette_ecm(c("Vascular ECM", "Basement ECM", "Interstitial ECM"))
#' }
#'
#' @export
palette_ecm <- function(levels_vector) {
  color_spec <- c(
    Vascular = "#d42626ff",
    Basement = "#6a3d9aff",
    Interstitial = "#2b9e2bff",
    not.assigned = "#7f7f7f"
  )

  color_map <- character(length(levels_vector))
  names(color_map) <- levels_vector

  for (level in levels_vector) {
    if (grepl("Vascular", level, ignore.case = TRUE)) {
      color_map[level] <- color_spec["Vascular"]
    } else if (grepl("Basement", level, ignore.case = TRUE)) {
      color_map[level] <- color_spec["Basement"]
    } else if (grepl("Interstitial", level, ignore.case = TRUE)) {
      color_map[level] <- color_spec["Interstitial"]
    } else if (level == "not.assigned") {
      color_map[level] <- color_spec["not.assigned"]
    } else {
      # Fallback for unrecognized annotations
      color_map[level] <- "#CCCCCC"
    }
  }

  return(color_map)
}
