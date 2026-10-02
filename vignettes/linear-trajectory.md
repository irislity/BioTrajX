Evaluating Pseudotime Methods with BioTrajX (Linear Datasets)
================

This article walks through applying BioTrajX’s DOE metrics to a CD8⁺
T-cell exhaustion dataset: head and neck squamous cell carcinoma
tumor-infiltrating and peripheral-blood CD8⁺ T cells (GSE164690),
annotated with ProjecTILs into functional states, of which this vignette
uses four — `CD8.NaiveLike`, `CD8.EM`, `CD8.TPEX`, and `CD8.TEX` — that
define the expected naive-to-exhausted progression; the dataset’s fifth
state, `CD8.CM`, is ignored. The dataset ships with pseudotime from
eight trajectory inference (TI) methods already computed and stored as
`meta.data` columns, so this vignette starts directly from that
pre-computed Seurat object rather than re-running any TI method.

## 1. Load the data

``` r
library(BioTrajX)
library(Seurat)

cd8t <- readRDS("data/cd8t.rds")

# Order the states naive-to-exhausted rather than alphabetically, so plots
# below read low-to-high along the trajectory. CD8.CM is left out of the
# levels (and so ignored throughout this vignette) rather than ordered in.
cd8t$functional.cluster <- factor(
  cd8t$functional.cluster,
  levels = c("CD8.NaiveLike", "CD8.EM", "CD8.TPEX", "CD8.TEX")
)
cd8t
#> An object of class Seurat 
#> 33545 features across 1357 samples within 1 assay 
#> Active assay: RNA (33545 features, 2000 variable features)
#>  3 layers present: data, counts, scale.data
#>  2 dimensional reductions calculated: pca, umap
```

`cd8t$functional.cluster` holds the four ProjecTILs states this vignette
uses, and `cd8t@meta.data` already has one numeric column per TI method:
`Slingshot`, `CytoTRACE`, `Monocle3`, `DPT`, `SCORPIUS`, `TSCAN`,
`PAGA-DPT`, and `Palantir`.

## 2. Visualize labels and pseudotime

One UMAP panel for the ground-truth cell states, plus one panel per TI
method colored by that method’s pseudotime:

``` r
pt_cols <- c("Slingshot", "CytoTRACE", "Monocle3", "DPT",
             "SCORPIUS", "TSCAN", "PAGA-DPT", "Palantir")

umap_df <- as.data.frame(Embeddings(cd8t, "umap"))
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$CellType <- cd8t$functional.cluster[rownames(umap_df)]

p_celltype <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP1, UMAP2, colour = CellType)) +
  ggplot2::geom_point(size = 0.5, alpha = 0.7) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::labs(title = "Cell type", colour = NULL)

method_plots <- lapply(pt_cols, function(method) {
  df <- cbind(umap_df[, c("UMAP1", "UMAP2")], Pseudotime = cd8t@meta.data[[method]])
  ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2, colour = Pseudotime)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.7) +
    # Fixed limits (pseudotime is always on a [0, 1] scale here) so every
    # panel's colour scale is identical — that's what lets patchwork collapse
    # them into one shared legend instead of repeating it on every panel.
    ggplot2::scale_colour_distiller(palette = "Spectral", direction = -1,
                                     na.value = "grey80", limits = c(0, 1)) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::labs(title = method)
})

if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(c(list(p_celltype), method_plots), ncol = 3) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")
} else {
  p_celltype
  method_plots
}
```

![](linear-trajectory_files/figure-gfm/visualize-1.png)<!-- -->

## 3. Marker genes

BioTrajX scores D and E using an “early” marker set (expected high in
naive cells) and a “terminal” marker set (expected high in exhausted
cells). This vignette combines two sources into one marker set: the
`CD8+ Tn` (naive) panel from the pan-cancer T-cell atlas shipped by
[SlimR](https://github.com/zhaoqing-wang/SlimR) for `early`, and the
MSigDB signature `GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN` (pulled
via `get_markers_msigdb()`) for `terminal`. Neither package is a
BioTrajX dependency, so install both separately first:

``` r
install.packages("SlimR")
install.packages("msigdbr")
```

``` r
pctit_s9  <- SlimR::Markers_list_PCTIT
msigdb_s9 <- get_markers_msigdb(
  early      = "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN",
  terminal   = "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN",
  collection = NULL,
  species    = "Homo sapiens"
)

# Combine the SlimR "early" panel with the MSigDB "terminal" signature into one
# marker_set. filter_markers() only dispatches on the class, so building it by
# hand with structure() — rather than the internal, unexported
# BioTrajX:::.marker_set() — works identically to a get_markers_*() result.
ms_s9 <- structure(
  list(early    = pctit_s9[["CD8+ Tn"]]$Markers,
       terminal = msigdb_s9$terminal,
       source   = "SlimR_PCTIT + MSigDB",
       metadata = list(early_set    = "CD8+ Tn",
                        terminal_set = "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN")),
  class = "marker_set"
)

ms_s9 <- filter_markers(ms_s9, cd8t, top_n = NULL, min_detection = 0.10)
early_markers    <- ms_s9$early
terminal_markers <- ms_s9$terminal
c(early_markers = length(early_markers), terminal_markers = length(terminal_markers))
#>    early_markers terminal_markers 
#>               20               23
```

`get_markers_cellmarker()` is available as another database-backed
source of marker genes if you’d rather not depend on SlimR or msigdbr.

## 4. Evaluate a single pseudotime trajectory

With markers in hand, `compute_single_doe_linear()` scores one
pseudotime vector against the three DOE components — Directionality,
Order, and Endpoints:

``` r
res_single <- compute_single_doe_linear(
  expr_or_seurat   = cd8t,
  pseudotime       = cd8t$Monocle3,
  early_markers    = early_markers,
  terminal_markers = terminal_markers,
  cluster_labels   = "functional.cluster",
  E_method         = "gmm"
)

res_single$DOE_score
#> [1] 0.5432022
plot(res_single, type = "bar")
```

![](linear-trajectory_files/figure-gfm/single-doe-1.png)<!-- -->

``` r
plot(res_single, type = "radar")
```

![](linear-trajectory_files/figure-gfm/single-doe-2.png)<!-- -->

``` r
plot(res_single, type = "heatmap")
```

![](linear-trajectory_files/figure-gfm/single-doe-3.png)<!-- -->

## 5. Compare all eight methods at once

`compute_multi_doe_linear()` runs the same evaluation across every
method in one call and ranks them. Pseudotime vectors are passed as a
named list; methods with too few non-missing values (here, `TSCAN` has
some `NA`s) are dropped first.

``` r
pt_list <- setNames(lapply(pt_cols, function(col) cd8t@meta.data[[col]]), pt_cols)
pt_list <- pt_list[sapply(pt_list, function(x) sum(!is.na(x)) > 10)]

res <- compute_multi_doe_linear(
  expr_or_seurat   = cd8t,
  pseudotime_list  = pt_list,
  early_markers    = early_markers,
  terminal_markers = terminal_markers,
  cluster_labels   = "functional.cluster",
  E_method         = "gmm"
)
#> Computing DOE metrics for 8 trajectories...
#> ==========================================================
#> 
#> --- Processing trajectory: Slingshot ---
#> 
#> --- Processing trajectory: CytoTRACE ---
#> 
#> --- Processing trajectory: Monocle3 ---
#> 
#> --- Processing trajectory: DPT ---
#> 
#> --- Processing trajectory: SCORPIUS ---
#> 
#> --- Processing trajectory: TSCAN ---
#> 
#> --- Processing trajectory: PAGA-DPT ---
#> 
#> --- Processing trajectory: Palantir ---
#> 
#> ==========================================================
#> Multi-trajectory DOE analysis complete!
#> Best trajectory: TSCAN (DOE score: 0.598)

res$comparison_summary
#>   trajectory   D_early    D_term    D_comp         O O_orientation   E_early
#> 1      TSCAN 0.7070410 0.4508281 0.5789345 0.7297995             + 0.4342105
#> 2   Monocle3 0.7053013 0.3740680 0.5396846 0.6257121             + 0.4365782
#> 3   PAGA-DPT 0.6922741 0.2850180 0.4886460 0.6530323             + 0.3716814
#> 4  CytoTRACE 0.7840763 0.1918649 0.4879706 0.6114504             + 0.4867257
#> 5  Slingshot 0.7189100 0.3226085 0.5207593 0.4716572             + 0.4336283
#> 6   Palantir 0.6088216 0.1992569 0.4040393 0.6275072             + 0.3392330
#> 7   SCORPIUS 0.7148566 0.2755311 0.4951938 0.5421092             + 0.3864307
#> 8        DPT 0.7308194 0.3092885 0.5200540 0.2625789             + 0.4129794
#>      E_term    E_comp DOE_score has_error D_early_rank D_term_rank D_comp_rank
#> 1 0.5526316 0.4863158 0.5983500     FALSE            5           1           1
#> 2 0.4955752 0.4642097 0.5432022     FALSE            6           2           2
#> 3 0.5162242 0.4321877 0.5246220     FALSE            7           5           6
#> 4 0.3156342 0.3829386 0.4941199     FALSE            1           8           7
#> 5 0.5103245 0.4688606 0.4870924     FALSE            3           3           3
#> 6 0.5132743 0.4084882 0.4800115     FALSE            8           7           8
#> 7 0.3008850 0.3383342 0.4585457     FALSE            4           6           5
#> 8 0.4867257 0.4468301 0.4098210     FALSE            2           4           4
#>   O_rank E_early_rank E_term_rank E_comp_rank DOE_score_rank
#> 1      1            3           1           1              1
#> 2      4            2           5           3              2
#> 3      2            7           2           5              3
#> 4      5            1           7           7              4
#> 5      7            4           4           2              5
#> 6      3            8           3           6              6
#> 7      6            6           8           8              7
#> 8      8            5           6           4              8
res$best_trajectory
#> [1] "TSCAN"
plot(res, type = "bar")
```

![](linear-trajectory_files/figure-gfm/multi-doe-1.png)<!-- -->

``` r
plot(res, type = "radar")
```

![](linear-trajectory_files/figure-gfm/multi-doe-2.png)<!-- -->

``` r
plot(res, type = "heatmap")
```

![](linear-trajectory_files/figure-gfm/multi-doe-3.png)<!-- -->

## 6. Would biologically implausible root gives a low DOE score?

A pseudotime trajectory rooted at the wrong starting cell runs
“backwards” relative to the true biology. To check whether DOE catches
that, we re-run Monocle3 from two different root cells — the
`CD8.NaiveLike` centroid (the biologically correct start) and the
`CD8.TEX` centroid (the wrong end of the trajectory) — and store both as
new `meta.data` columns alongside the eight pre-computed methods. This
needs `monocle3` and `SingleCellExperiment`, neither of which is a
BioTrajX dependency:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")

# monocle3 is distributed from GitHub, not CRAN/Bioconductor
BiocManager::install(c("BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma",
                        "lme4", "S4Vectors", "SummarizedExperiment", "batchelor",
                        "HDF5Array", "terra", "ggrastr"))
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("cole-trapnell-lab/monocle3")
```

If `install_github()` fails with something like
`namespace 'rlang' ... is already loaded, but >= ... is required`,
restart R first — an old `rlang` already loaded in your session can’t be
swapped out mid-session. After restarting, run
`install.packages("rlang")` *before* loading any other package, then
retry the `monocle3` install.

``` r
# Monocle3 is run on the top 2000 highly-variable genes, not the full gene
# set — this matches standard practice and is what the manuscript's own
# root-robustness analysis uses.
cd8t <- FindVariableFeatures(cd8t, nfeatures = 2000, verbose = FALSE)
hvg  <- VariableFeatures(cd8t)

centroid_cell <- function(obj, cluster, label_col = "functional.cluster", n_pcs = 20) {
  pca   <- Embeddings(obj, "pca")[, 1:n_pcs]
  idx   <- which(obj@meta.data[[label_col]] == cluster)
  ctr   <- colMeans(pca[idx, , drop = FALSE])
  dists <- rowSums(sweep(pca[idx, , drop = FALSE], 2, ctr)^2)
  colnames(obj)[idx[which.min(dists)]]
}

root_naive <- centroid_cell(cd8t, "CD8.NaiveLike")
root_tex   <- centroid_cell(cd8t, "CD8.TEX")

run_monocle3_from_root <- function(seurat, start_cell, genes) {
  expr <- as.matrix(GetAssayData(seurat, assay = "RNA", layer = "data"))[genes, ]
  cds  <- monocle3::new_cell_data_set(
    Matrix::Matrix(expr, sparse = TRUE),
    cell_metadata = data.frame(cell_id = colnames(expr), row.names = colnames(expr)),
    gene_metadata = data.frame(gene_short_name = rownames(expr), row.names = rownames(expr))
  )
  cds <- monocle3::preprocess_cds(cds, num_dim = min(30L, ncol(expr) - 1L))
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
  cds <- monocle3::cluster_cells(cds)
  cds <- monocle3::learn_graph(cds, use_partition = FALSE)
  cds <- monocle3::order_cells(cds, root_cells = start_cell)

  pt <- monocle3::pseudotime(cds)
  pt[is.infinite(pt)] <- NA_real_
  rng <- range(pt, na.rm = TRUE)
  (pt - rng[1]) / (rng[2] - rng[1])
}

pt_root_naive <- run_monocle3_from_root(cd8t, root_naive, hvg)
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
pt_root_tex   <- run_monocle3_from_root(cd8t, root_tex, hvg)
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%

cd8t <- AddMetaData(cd8t, pt_root_naive[colnames(cd8t)], col.name = "Monocle3_root_naive")
cd8t <- AddMetaData(cd8t, pt_root_tex[colnames(cd8t)],   col.name = "Monocle3_root_tex")
```

``` r
res_root_naive <- compute_single_doe_linear(
  expr_or_seurat = cd8t, pseudotime = cd8t$Monocle3_root_naive,
  early_markers = early_markers, terminal_markers = terminal_markers,
  cluster_labels = "functional.cluster", E_method = "gmm"
)
res_root_tex <- compute_single_doe_linear(
  expr_or_seurat = cd8t, pseudotime = cd8t$Monocle3_root_tex,
  early_markers = early_markers, terminal_markers = terminal_markers,
  cluster_labels = "functional.cluster", E_method = "gmm"
)

data.frame(
  root      = c("CD8.NaiveLike centroid", "CD8.TEX centroid"),
  DOE_score = c(res_root_naive$DOE_score, res_root_tex$DOE_score)
)
#>                     root DOE_score
#> 1 CD8.NaiveLike centroid 0.5432022
#> 2       CD8.TEX centroid 0.1983723
```

Rooting Monocle3 at the exhausted end rather than the naive end lowers
the DOE score substantially — the same signal that would flag a
genuinely mis-rooted TI run on real data.
