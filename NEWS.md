# BioTrajX 0.1.0

First public release.

## Core metrics

* `metrics_d()`, `metrics_o()`, `metrics_e()` implement the Directionality,
  Order, and Endpoint components of the DOE framework for scoring a
  pseudotime trajectory against early/terminal marker programs and (for
  Endpoint) ground-truth or GMM-derived cell states.
* `compute_single_doe_linear()` / `compute_multi_doe_linear()` wrap D/O/E
  into a single call for one or many pseudotime vectors over a linear
  trajectory, with ranking and summary tables.
* `compute_single_doe_branched()` / `compute_multi_doe_branched()` extend
  this to branched trajectories, scoring each branch separately and
  combining them into overall summaries.
* `plot_metrics_d()`, `plot_metrics_o()`, `plot_metrics_e()`, and
  `plot()` methods for `doe_results`, `multi_doe_results`,
  `single_doe_branched`, and `multi_doe_branched` provide bar, radar, and
  heatmap visualizations.
* `reverse_pseudotime()` helper for testing orientation sensitivity.

## Marker utilities

* `get_markers_msigdb()`, `get_markers_go()`, and `get_markers_cellmarker()`
  pull early/terminal marker gene sets from MSigDB, GO, and CellMarker.
* `filter_markers()` filters any `marker_set` (from the helpers above, or
  built by hand) down to genes detected at a minimum rate in the data.
* `markers_to_list()` converts a `marker_set` or data frame into the
  plain list form the DOE functions expect.

## Tutorials

* Three vignettes — linear trajectories, branched trajectories, and
  ground-truth validation — walk through applying DOE metrics to real
  CD8+ T-cell exhaustion, stem-to-lineage, and LCMV time-course datasets.
