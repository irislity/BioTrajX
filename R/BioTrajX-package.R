#' BioTrajX: Trajectory Evaluation Metrics for Single-Cell Data
#'
#' BioTrajX provides tools to assess pseudotime and trajectory inference outputs
#' using biologically grounded expectations. It implements the DOPE metrics
#' (Directionality, Order Consistency, Program Coherence, Endpoint Validity),
#' includes plotting utilities for visualization.
#'
#' @description
#' Compute all four DOPE metrics (Directionality, Order, Program coherence, Endpoints validity)
#' for multiple pseudotime trajectories in a single function call. This function loops through
#' each provided pseudotime and computes DOPE metrics, allowing comparison across different
#' trajectory inference methods or parameters. Supports linear or branched trajectories; for
#' branched analyses, cells can be subset using cluster labels before metric computation.
#'
#' @section Main Functions:
#' -
#' -
#' -
#' -
#'
#' @section Individual DOPE Functions:
#' - `metrics_D()` : Directionality
#' - `metrics_O()` : Order consistency
#' - `metrics_P()` : Program coherence
#' - `metrics_E()` : Endpoint validity
#'
#'
#' @docType _PACKAGE
#' @name BioTrajX
#' @keywords internal
NULL
