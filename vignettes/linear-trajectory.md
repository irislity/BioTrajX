Evaluating Pseudotime Methods with BioTrajX (Linear Datasets)
================

This article walks through applying BioTrajX’s DOE metrics to a real
CD8⁺ T cell exhaustion dataset (Naive → Terminal exhaustion), comparing
multiple pseudotime methods using the D, O, and E metrics.

## Supported trajectory inference methods

BioTrajX evaluates pseudotime from any TI method. The seven methods
tested in this package are listed below. You are not required to run all
of them — pass whichever pseudotime vectors you have to
`compute_multi_doe_linear()`.

| Method | Type | Language | Notes |
|----|----|----|----|
| Slingshot | MST + principal curves | R | Recommended for linear trajectories |
| CytoTRACE | Gene count structure | R | No root cell needed |
| CytoTRACE2 | Deep learning | R | Requires raw counts |
| DPT | Diffusion pseudotime | R | Requires a root cell |
| SCORPIUS | Dimensionality reduction | R | Designed for linear trajectories |
| PAGA | Graph abstraction + DPT | Python | Requires `scanpy` conda env |
| CellRank | Fate probabilities | Python | Requires `cellrank` conda env |

To run all methods at once, source
`manuscript/scripts/real/run_ti_methods.R`, call `run_all_ti_methods()`,
and add the results to your Seurat object:

``` r
source("manuscript/scripts/real/run_ti_methods.R")

# Identify a root cell (lowest naive marker expression)
expr       <- as.matrix(GetAssayData(cd8t, layer = "data"))
start_cell <- colnames(expr)[which.min(colMeans(expr[early, , drop = FALSE]))]

ti_results <- run_all_ti_methods(expr, start_cell = start_cell)

# Add all pseudotime vectors as metadata columns
cd8t <- AddMetaData(cd8t, ti_results)
```

## 1. Load the CD8⁺ T cell dataset

Download the Seurat object manually from the following link:

<a href="https://transfer.it/t/PFoL2d0hBZAT"
class="uri"><strong>https://transfer.it/t/PFoL2d0hBZAT</strong></a>

Then load it into your R session:

``` r
library(BioTrajX)
library(Seurat)

cd8t <- readRDS("data/cd8t.rds")
```

## 2. Visualize pseudotime and cluster annotations

``` r
DimPlot(cd8t, group.by = "functional.cluster")
FeaturePlot(cd8t, features = "monocle3_pseudotime")
FeaturePlot(cd8t, features = "slingshot_pseudotime")
FeaturePlot(cd8t, features = "cytotrace_score")
```

<embed src="figures/linear/umap_cluster.pdf" width="100%" type="application/pdf" />

<embed src="figures/linear/umap_slingshot.pdf" width="100%" type="application/pdf" />

<embed src="figures/linear/umap_monocle3.pdf" width="100%" type="application/pdf" />

<embed src="figures/linear/umap_cytotrace.pdf" width="100%" type="application/pdf" />

## 3. Prepare expression and marker genes

``` r
expr <- as.matrix(GetAssayData(cd8t, layer = "data"))

# Early-like markers: selected by strongest negative Spearman correlation with
# slingshot pseudotime among biologically interpretable genes (ribosomal/housekeeping
# excluded). These genes monotonically decrease along the CD8+ T cell exhaustion
# trajectory, analogous to early genes in the simulation vignette.
# ModuleScore: CD8.NaiveLike = +0.38, all exhausted clusters negative (to -0.65).
# Refs: Wherry & Kurachi 2015 (Nat Rev Immunol); Miller et al. 2019 (Nat Immunol)
early <- c(
  "LTB", "LDHB", "LEF1", "KLF2", "MAL",
  "SELL", "IL7R", "CCR7", "TCF7"
)

# Terminal markers: selected by strongest positive Spearman correlation with
# slingshot pseudotime — effector/exhaustion genes that monotonically increase
# along the trajectory (score: NaiveLike = -0.54, TEX = +1.30).
# Refs: Sade-Feldman et al. 2018 (Cell); Guo et al. 2018 (Nat Med)
term <- c(
  "SRGN", "GZMK", "CCL5", "NKG7", "RGS1",
  "CST7", "GZMA", "CXCL13", "TIGIT", "LAG3"
)

cluster <- cd8t$functional.cluster
```

``` r
early <- c("CCR7", "TCF7", "IL7R")
term  <- c("PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT", "CTLA4", "CXCL13", "ENTPD1", "ITGAE")

cd8t <- AddModuleScore(cd8t,
  features = early,
  name = 'early_markers'
)

cd8t <- AddModuleScore(cd8t,
  features = term,
  name = 'term_markers'
)
```

``` r
FeaturePlot(cd8t, features = "early_markers1") + ggplot2::ggtitle("Early Marker Expression")
FeaturePlot(cd8t, features = "term_markers1") + ggplot2::ggtitle("Terminal Marker Expression")
```

### Alternative: retrieve markers from MSigDB

The vectors above were curated manually from the literature. As an
alternative,\
`get_markers_msigdb()` can pull the same biological programs directly
from MSigDB\
C2 (curated gene sets). The result is a `marker_set` whose `$early` and\
`$terminal` slots are drop-in replacements for the character vectors
above.

``` r
# install.packages("msigdbr")
ms <- get_markers_msigdb(
  early      = "GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_UP",
  terminal   = "GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_UP",
  collection = "C7"   # immunologic signatures
)

# equivalent to passing `early` and `term` by hand
res_single <- compute_single_doe_linear(
  cd8t,
  cd8t$slingshot_pseudotime,
  early_markers    = ms$early,
  terminal_markers = ms$terminal
)
```

You can also search by keyword when you do not know the exact set name:

``` r
ms <- get_markers_msigdb(
  early      = "NAIVE.*CD8",
  terminal   = "EXHAUST.*CD8",
  collection = "C7",
  exact      = FALSE   # treat strings as grepl() patterns
)
```

## 4. Evaluate a single pseudotime trajectory

``` r
res_single <- compute_single_doe_linear(
  cd8t,
  cd8t$slingshot_pseudotime,
  early_markers    = early,
  terminal_markers = term
)
```

## 5. Visualize single-trajectory DOE results

``` r
plot(res_single, type = "bar")
plot(res_single, type = "radar")
plot(res_single, type = "heatmap")
```

<embed src="figures/linear/single_doe_bar.pdf" width="100%" type="application/pdf" />

<embed src="figures/linear/single_doe_radar.pdf" width="100%" type="application/pdf" />

<embed src="figures/linear/single_doe_heatmap.pdf" width="100%" type="application/pdf" />

## 6. Evaluate multiple pseudotime methods simultaneously

``` r
pseudotime_methods <- list(
  Slingshot  = cd8t$Slingshot,
  CytoTRACE  = cd8t$CytoTRACE,
  CytoTRACE2 = cd8t$CytoTRACE2,
  Monocle3   = cd8t$Monocle3,
  DPT        = cd8t$DPT,
  SCORPIUS   = cd8t$SCORPIUS,
  PAGA       = cd8t$PAGA,
  CellRank   = cd8t$CellRank
)

res <- compute_multi_doe_linear(
  expr,
  pseudotime_list  = pseudotime_methods,
  early_markers    = early,
  terminal_markers = term
)
```

## 7. Visualize multi-trajectory DOE results

``` r
plot(res, type = "bar")
plot(res, type = "radar")
plot(res, type = "heatmap")
```

<embed src="figures/linear/doe_bar.pdf" width="100%" type="application/pdf" />

<embed src="figures/linear/doe_radar.pdf" width="100%" type="application/pdf" />

<embed src="figures/linear/doe_heatmap.pdf" width="100%" type="application/pdf" />
