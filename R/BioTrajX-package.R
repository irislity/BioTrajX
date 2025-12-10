#' BioTrajX: Trajectory Evaluation Metrics for Single-Cell Data
#'
#' BioTrajX provides tools to assess pseudotime and trajectory inference outputs
#' using biologically grounded expectations. It implements the DOE metrics
#' (Directionality, Order Consistency, Endpoint Validity) and
#' includes plotting utilities for visualization.
#'
#' @description
#' Compute the core DOE metrics (Directionality, Order, Endpoints validity)
#' across multiple pseudotime trajectories. Supports linear and branched
#' trajectories.
#'
#' @section Core Metrics:
#' \itemize{
#'   \item \code{metrics_D()} Directionality
#'   \item \code{metrics_E()} Endpoint validity
#'   \item \code{metrics_O()} Order consistency
#' }
#'
#'
#' @section DOE Computation Wrappers:
#' \itemize{
#'   \item \code{compute_single_DOE_linear()}
#'   \item \code{compute_single_DOE_branched()}
#'   \item \code{compute_multi_DOE_linear()}
#'   \item \code{compute_multi_DOE_branched()}
#' }
#'
#'
#' @section Additional Utilities:
#' \itemize{
#'   \item \code{reverse_pseudotime()}
#'   \item \code{plot_metrics_D()}
#'   \item \code{plot_metrics_E()}
#'   \item \code{plot_metrics_O()}
#' }
#'
#'
#' @docType package
#' @name BioTrajX
#' @keywords internal
"_PACKAGE"