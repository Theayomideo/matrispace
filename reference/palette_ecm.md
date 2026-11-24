# Generate ECM niche color palette

Creates a color palette for ECM niche annotations with biologically
meaningful colors: Vascular (red), Basement (purple), Interstitial
(green).

## Usage

``` r
palette_ecm(levels_vector)
```

## Arguments

- levels_vector:

  Character vector of ECM niche levels

## Value

Named character vector of hex colors

## Details

ECM Niche Color Palette:

- Vascular (#d42626): Red - blood vessels and vascular ECM

- Basement (#6a3d9a): Purple - structural foundation

- Interstitial (#2b9e2b): Green - matrix between cells

- not.assigned (#7f7f7f): Grey - unassigned spots

## Examples

``` r
if (FALSE) { # \dontrun{
colors <- palette_ecm(c("Vascular ECM", "Basement ECM", "Interstitial ECM"))
} # }
```
