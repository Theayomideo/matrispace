# MatriComDB Matrisome Pair Database

Database of unique heterotypic MatriComDB pairs for spatial
co-expression analysis. Every pair contains at least one matrisome gene.
The source data span multiple MatriComDB interaction classes, including
annotated ECM communication records and structural ECM interactions. A
row should not be interpreted as signaling unless its annotation
supports that interpretation.

## Usage

``` r
matrisome_pairs
```

## Format

A data frame with columns:

- Gene1:

  First pair-member gene symbol

- Gene2:

  Second pair-member gene symbol

- Source:

  Database source for the interaction

- Interaction_Type:

  Type of interaction

- Database_Score:

  Confidence score from source database

## Details

Reciprocal records are collapsed to one alphabetically ordered pair and
homomeric records are excluded.

## Examples

``` r
data(matrisome_pairs)
head(matrisome_pairs)
#> # A tibble: 6 × 5
#>   Gene1         Gene2 Interaction_Type Source            Database_Score
#>   <chr>         <chr> <chr>            <chr>                      <dbl>
#> 1 1810065E05RIK VTN   Cell-Matrix      STRING_physical                1
#> 2 2210010C04RIK AMBP  Cell-Matrix      STRING_physical                1
#> 3 4933400A11RIK S100B Cell-Matrix      STRING_physical                1
#> 4 A             PLG   Cell-Matrix      STRING_physical                1
#> 5 A1BG          FGB   Cell-Matrix      STRING_physical                1
#> 6 A2M           ANXA7 Matrix-Matrix    MatrixDB IMEx ext              1
nrow(matrisome_pairs)
#> [1] 16876
```
