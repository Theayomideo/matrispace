# Package index

## Package overview

- [`matrispace-package`](https://Theayomideo.github.io/matrispace/reference/matrispace-package.md)
  [`matrispace`](https://Theayomideo.github.io/matrispace/reference/matrispace-package.md)
  : matrispace: Identification and Visualization of Spatially-Resolved
  ECM Gene Expression Patterns

## Preprocessing

- [`prepare_object()`](https://Theayomideo.github.io/matrispace/reference/prepare_object.md)
  : Prepare Seurat object for matrisome analysis
- [`convert_spe_to_seurat()`](https://Theayomideo.github.io/matrispace/reference/convert_spe_to_seurat.md)
  : Convert SpatialExperiment to Seurat object
- [`translate_signatures()`](https://Theayomideo.github.io/matrispace/reference/translate_signatures.md)
  : Translate gene signatures using HGNChelper

## Matrisome and ECM niches

- [`score_matrisome()`](https://Theayomideo.github.io/matrispace/reference/score_matrisome.md)
  : Calculate matrisome feature scores
- [`annotate_ecm_niches()`](https://Theayomideo.github.io/matrispace/reference/annotate_ecm_niches.md)
  : Annotate ECM niches using ScType
- [`compute_expression_metrics()`](https://Theayomideo.github.io/matrispace/reference/compute_expression_metrics.md)
  : Compute expression metrics for matrisome categories

## Differential expression

- [`find_ecm_markers()`](https://Theayomideo.github.io/matrispace/reference/find_ecm_markers.md)
  : Find differentially expressed ECM genes between groups
- [`find_all_ecm_markers()`](https://Theayomideo.github.io/matrispace/reference/find_all_ecm_markers.md)
  : Find ECM markers for all clusters
- [`find_spatially_corrected_ecm_markers()`](https://Theayomideo.github.io/matrispace/reference/find_spatially_corrected_ecm_markers.md)
  : Find spatially corrected ECM markers between two groups
- [`find_ecm_de_samples()`](https://Theayomideo.github.io/matrispace/reference/find_ecm_de_samples.md)
  : Differential expression of ECM genes between samples or conditions

## Spatial statistics

- [`compute_lisa()`](https://Theayomideo.github.io/matrispace/reference/compute_lisa.md)
  : Compute LISA statistics for feature pairs
- [`compute_morans_i()`](https://Theayomideo.github.io/matrispace/reference/compute_morans_i.md)
  : Compute Moran's I spatial autocorrelation

## Matrisome-pair co-expression

- [`score_matrisome_pair_coexpression()`](https://Theayomideo.github.io/matrispace/reference/score_matrisome_pair_coexpression.md)
  : Calculate spatial matrisome-pair co-expression scores
- [`find_matrisome_pair_enrichment()`](https://Theayomideo.github.io/matrispace/reference/find_matrisome_pair_enrichment.md)
  : Find enriched matrisome-pair co-expression per group
- [`aggregate_matrisome_pairs()`](https://Theayomideo.github.io/matrispace/reference/aggregate_matrisome_pairs.md)
  : Aggregate reciprocal rows into unique matrisome pairs

## Visualization

- [`plot_spatial_feature()`](https://Theayomideo.github.io/matrispace/reference/plot_spatial_feature.md)
  : Create spatial feature plot
- [`plot_matrisome()`](https://Theayomideo.github.io/matrispace/reference/plot_matrisome.md)
  : Create matrisome category spatial plots
- [`plot_spatial_blend()`](https://Theayomideo.github.io/matrispace/reference/plot_spatial_blend.md)
  : Create blended spatial plot for two features
- [`plot_lisa()`](https://Theayomideo.github.io/matrispace/reference/plot_lisa.md)
  : Create LISA clustering plot
- [`palette_ecm()`](https://Theayomideo.github.io/matrispace/reference/palette_ecm.md)
  : Generate ECM niche color palette
- [`palette_annotation()`](https://Theayomideo.github.io/matrispace/reference/palette_annotation.md)
  : Generate annotation color palette
- [`add_scalebar()`](https://Theayomideo.github.io/matrispace/reference/add_scalebar.md)
  : Add scale bar to spatial plot

## Reference data

- [`matrisome_db`](https://Theayomideo.github.io/matrispace/reference/matrisome_db.md)
  : Matrisome Gene Database
- [`matrisome_signatures`](https://Theayomideo.github.io/matrispace/reference/matrisome_signatures.md)
  : Matrisome UCell Signatures
- [`ecm_signatures`](https://Theayomideo.github.io/matrispace/reference/ecm_signatures.md)
  : ECM Niche Signatures
- [`matrisome_pairs`](https://Theayomideo.github.io/matrispace/reference/matrisome_pairs.md)
  : MatriComDB Matrisome Pair Database
- [`ecm_markers`](https://Theayomideo.github.io/matrispace/reference/ecm_markers.md)
  : ECM Niche Markers

## Deprecated

- [`annotate_ecm_domains()`](https://Theayomideo.github.io/matrispace/reference/annotate_ecm_domains.md)
  : Annotate ECM domains (deprecated)
