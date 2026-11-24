# matrispace: Identification and Visualization of Spatially-Resolved ECM Gene Expression Patterns

matrispace is the R package companion to MatriSpace. It brings the same
extracellular matrix-focused workflow used in the web application to
local, scriptable analysis of spatial transcriptomics data. The package
includes functions for:

- Matrisome scoring using UCell signatures

- ECM niche classification using the ScType algorithm

- Spatial statistics (LISA, Moran's I)

- Spatial co-expression analysis of unique heterotypic matrisome pairs

- Static visualizations for spatial data

## Main Analysis Functions

- [`prepare_object`](https://Theayomideo.github.io/matrispace/reference/prepare_object.md):

  Master preprocessing pipeline

- [`score_matrisome`](https://Theayomideo.github.io/matrispace/reference/score_matrisome.md):

  Calculate matrisome category scores

- [`annotate_ecm_niches`](https://Theayomideo.github.io/matrispace/reference/annotate_ecm_niches.md):

  Classify ECM niches with ScType

- [`compute_lisa`](https://Theayomideo.github.io/matrispace/reference/compute_lisa.md):

  Local spatial association statistics

- [`compute_morans_i`](https://Theayomideo.github.io/matrispace/reference/compute_morans_i.md):

  Global spatial autocorrelation

- [`score_matrisome_pair_coexpression`](https://Theayomideo.github.io/matrispace/reference/score_matrisome_pair_coexpression.md):

  Matrisome-pair co-expression scoring

## Data Objects

- [`matrisome_db`](https://Theayomideo.github.io/matrispace/reference/matrisome_db.md):

  Gene-to-category mapping database

- [`matrisome_signatures`](https://Theayomideo.github.io/matrispace/reference/matrisome_signatures.md):

  UCell gene signatures

- [`ecm_signatures`](https://Theayomideo.github.io/matrispace/reference/ecm_signatures.md):

  ECM niche signatures

- [`matrisome_pairs`](https://Theayomideo.github.io/matrispace/reference/matrisome_pairs.md):

  Unique heterotypic MatriComDB pairs

## Visualization Functions

- [`plot_spatial_feature`](https://Theayomideo.github.io/matrispace/reference/plot_spatial_feature.md):

  Single feature spatial plot

- [`plot_matrisome`](https://Theayomideo.github.io/matrispace/reference/plot_matrisome.md):

  Matrisome category plots

- [`plot_lisa`](https://Theayomideo.github.io/matrispace/reference/plot_lisa.md):

  LISA clustering visualization

- [`plot_spatial_blend`](https://Theayomideo.github.io/matrispace/reference/plot_spatial_blend.md):

  Two-feature blend plot

## References

<https://github.com/Theayomideo/matrispace>

## See also

The companion MatriSpace web application:
<https://izzilab-matrispace.share.connect.posit.cloud/>

## Author

Ayomide Oshinjo, Daiqing Chen, Petar Petrov, Valerio Izzi, Alexandra
Naba
