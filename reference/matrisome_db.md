# Matrisome Gene Database

A data frame mapping genes to matrisome categories and subcategories.
Contains 1,027 canonical genes plus F8 and USH2A for functional
signatures.

## Usage

``` r
matrisome_db
```

## Format

A data frame with columns:

- division:

  Main matrisome category (Core matrisome or Matrisome-associated)

- notes:

  Additional notes

- gene:

  HGNC gene symbol

- Gene Name:

  Full gene name

- Synonyms:

  Alternative gene names (pipe-separated)

- HGNC_IDs:

  HGNC identifier

- HGNC_IDs Links:

  Links to HGNC database

- UniProt_IDs:

  UniProt identifiers

- Refseq_IDs:

  RefSeq identifiers

- Notes_ignore:

  Internal notes (ignored)

- ecm_subcategory:

  Subcategory (e.g., Collagens, ECM Glycoproteins)

## Source

The Matrisome Project <https://matrisomeproject.org>

## Examples

``` r
data(matrisome_db)
head(matrisome_db)
#> # A tibble: 6 × 11
#>   division       notes      gene  `Gene Name` Synonyms HGNC_IDs `HGNC_IDs Links`
#>   <chr>          <chr>      <chr> <chr>       <chr>    <chr>    <chr>           
#> 1 Core matrisome ECM Glyco… ABI3… ABI family… FLJ4174… 17265    17265           
#> 2 Core matrisome ECM Glyco… ADIP… adiponecti… ACDC|AC… 13633    13633           
#> 3 Core matrisome ECM Glyco… AEBP1 AE binding… ACLP; F… 303      303             
#> 4 Core matrisome ECM Glyco… AGRN  agrin       FLJ45064 329      329             
#> 5 Core matrisome ECM Glyco… AMBN  ameloblast… -        452      452             
#> 6 Core matrisome ECM Glyco… AMELX amelogenin… AIH1|AL… 461      461             
#> # ℹ 4 more variables: UniProt_IDs <chr>, Refseq_IDs <chr>, Notes_ignore <chr>,
#> #   ecm_subcategory <chr>
```
