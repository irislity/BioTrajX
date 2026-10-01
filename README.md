# BioTrajX

BioTrajX is an easy-to-use, lightweight R package for evaluating whether single-cell trajectories are biologically plausible. It is independent of the trajectory inference method and is designed to compare alternative pseudotime orderings, branch-specific trajectories, and root-cell choices before downstream analysis.

BioTrajX provides Directionality, Order, and Endpoint (DOE) metrics for quantitatively evaluating single-cell pseudotime orderings by different trajectory inference methods or produced by different root cell against prior biological knowledge. The package includes wrappers for linear and branched trajectories, along with visualization tools for summarizing and comparing results across competing trajectory inference methods.

![](man/figures/overview.jpeg)

## Installation

Install the development version directly from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("irislity/BioTrajX")
```

## DOE metrics



The package implements three complementary metrics:

- **Directionality (D):** do early and terminal marker programs change in the
  expected direction along pseudotime?
- **Order consistency (O):** do marker programs change coherently and
  monotonically? The sign reports whether trajectory orientation agrees with
  the expected biology.
- **Endpoint validity (E):** are expected early and terminal cells enriched at
  the corresponding ends of the trajectory?

The composite `DOE_score` is the arithmetic mean of available D, O, and E
components. A high score is evidence that the trajectory agrees with the
biological prior encoded by the marker sets; it is not a replacement for
visual inspection or independent validation.

![](man/figures/DOE.jpg)
![](man/figures/doe_schematic.png)


## Getting Started 


Expression data should be log-normalized with genes in rows and cells in
columns as count matrix. Pseudotime vectors can come from any TI method and should be named by
cell.

```r
library(BioTrajX)
```

### Linear trajectories

```r
# compute DOE metrics for a single trajectory
result <- compute_single_doe_linear(expr, pseudotime, early_markers, terminal_markers)

# compare multiple trajectories
comparison <- compute_multi_doe_linear(expr, pseudotime_list, early_markers, terminal_markers)

# visualization of DOE scores
plot(comparison, type = "bar")
plot(comparison, type = "radar")
plot(comparison, type = "heatmap")
```

### Branched trajectories

For a branched trajectory, provide one early and one terminal marker set per
lineage and, when needed, restrict each branch to its own cell-state labels:

```r
# compare multiple trajectories
branched <- compute_multi_doe_branched(
  expr_or_seurat        = expr,
  pseudotime_list       = pseudotime_list,
  early_markers_list    = list(erythroid = stem_markers,
                               b_cell = stem_markers),
  terminal_markers_list = list(erythroid = erythroid_markers,
                               b_cell = b_cell_markers),
  cluster_labels        = cell_state,
  branch_filters        = branch_filters
)
# visualization of DOE scores
plot(branched, scope = "branch", type = "heatmap")
plot(branched, scope = "overall", type = "bar")
plot(branched, scope = "branch", type = "radar", branch_mode = "facet")
```

### Marker Sets
Marker sets can be supplied manually or retrieved with
`get_markers_msigdb()` and `get_markers_cellmarker()`. Use `filter_markers()`
to retain genes detected in the data before scoring.


## Tutorials

- [Linear trajectories](vignettes/linear-trajectory.md): compare pseudotime
  methods on CD8+ T-cell exhaustion, inspect marker trends, and compare roots.
- [Branched trajectories](vignettes/branched-trajectory.md): evaluate
  Stem-to-Erythroid and Stem-to-B-cell trajectories with lineage-specific
  marker sets and scores.
- [Ground-truth validation](vignettes/ground-truth-validation.md): test
  whether if the ranks of DOE are aligned known experimental time in an LCMV
  infection time course.



## Citation and license

If you use BioTrajX, please cite the accompanying manuscript. BioTrajX is
released under the MIT license.



