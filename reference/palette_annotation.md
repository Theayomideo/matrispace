# Generate annotation color palette

Creates a consistent color palette for categorical annotations, with
special handling for "not.assigned" values (always grey).

## Usage

``` r
palette_annotation(levels_vector)
```

## Arguments

- levels_vector:

  Character vector of annotation levels

## Value

Named character vector of hex colors

## Details

Uses the Glasbey palette for up to 32 distinct colors. For more
categories, falls back to viridis interpolation. "not.assigned" is
always mapped to grey.

## Examples

``` r
if (FALSE) { # \dontrun{
colors <- palette_annotation(c("Tumor", "Stroma", "not.assigned"))
} # }
```
