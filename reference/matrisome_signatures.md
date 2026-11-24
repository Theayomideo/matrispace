# Matrisome UCell Signatures

Named list of gene vectors for MatriSpace matrisome categories and
subcategories, formatted for use with UCell scoring.

## Usage

``` r
matrisome_signatures
```

## Format

A named list where:

- names:

  Signature names (e.g., "collagens", "ecm_glycoproteins")

- values:

  Character vectors of gene symbols

## Details

The six main categories exclude F8 and USH2A; functional signatures
retain them. Includes signatures for:

- 6 main categories: Collagens, ECM Glycoproteins, Proteoglycans, ECM
  Regulators, Secreted Factors, ECM-affiliated Proteins

- curated subcategories and gene families such as Basement Membrane,
  Laminins, Matricellular proteins, Mucins, and Peri-vascular ECM

## Examples

``` r
data(matrisome_signatures)
names(matrisome_signatures)
#>  [1] "ecm_glycoproteins"           "collagens"                  
#>  [3] "proteoglycans"               "ecm_regulators"             
#>  [5] "secreted_factors"            "ecm-affiliated_proteins"    
#>  [7] "basement_membrane"           "hemostasis"                 
#>  [9] "elastic_fibers"              "growth_factor-binding"      
#> [11] "laminin_-_basement_membrane" "matricellular"              
#> [13] "syndecan"                    "glypican"                   
#> [15] "annexin"                     "cathepsin"                  
#> [17] "ccn_family"                  "cystatin"                   
#> [19] "facit"                       "fibulin"                    
#> [21] "galectin"                    "mucin"                      
#> [23] "plexin"                      "semaphorin"                 
#> [25] "perivascular"               
```
