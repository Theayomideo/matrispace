# ECM Niche Markers

Marker genes for ECM niche classification using ScType.

## Usage

``` r
ecm_markers
```

## Format

A list with two elements:

- gs_positive:

  Positive marker genes per niche

- gs_negative:

  Negative marker genes per niche

## Examples

``` r
data(ecm_markers)
names(ecm_markers$gs_positive)
#> [1] "Interstitial_ECM" "Basement_ECM"     "Vascular_ECM"    
```
