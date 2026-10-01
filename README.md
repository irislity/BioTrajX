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


![](man/figures/DOE.jpg)
![](man/figures/doe_schematic.png)


## Getting Started 


Expression data should be log-normalized with genes in rows and cells in
columns as count matrix. Pseudotime vectors can come from any TI method and should be named by
cell.

```r
library(BioTrajX)
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



