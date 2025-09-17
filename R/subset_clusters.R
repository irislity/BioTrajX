#' Subset cells by clusters (Seurat or matrix only)
#'
#' @param expr_or_seurat Seurat object or expression matrix (genes x cells).
#' @param clusters Either:
#'   - character scalar: name of Seurat meta.data column with cluster labels, or
#'   - vector: cluster label per cell (same length/order as columns).
#' @param include Character vector of cluster names to keep (optional).
#' @param exclude Character vector of cluster names to drop (optional).
#' @param min_cells Integer; require at least this many cells after subsetting.
#' @param drop_unused_levels Logical; drop unused factor levels in the returned labels.
#' @param assay For Seurat: assay to set as DefaultAssay after subsetting (optional).
#'
#' @return list with:
#'   - obj: subsetted object of same class as input
#'   - idx: integer indices of kept cells (relative to original)
#'   - cluster_labels: labels of kept cells (character vector)
subset_by_clusters <- function(expr_or_seurat,
                               clusters,
                               include = NULL,
                               exclude = NULL,
                               min_cells = 1,
                               drop_unused_levels = TRUE,
                               assay = NULL) {
  stopifnot(min_cells >= 0)

  is_mat_like <- function(x) is.matrix(x) || inherits(x, "dgCMatrix") || is.data.frame(x)

  # ---- cell names ----
  get_cells <- function(x) {
    if (inherits(x, "Seurat")) {
      Seurat::Cells(x)
    } else if (is_mat_like(x)) {
      colnames(x)
    } else {
      stop("Unsupported class: ", paste(class(x), collapse = ", "),
           ". Only Seurat or matrix/data.frame/dgCMatrix are supported.")
    }
  }

  # ---- labels ----
  get_clusters <- function(x, clusters) {
    n_cells <- length(get_cells(x))
    if (inherits(x, "Seurat") && is.character(clusters) && length(clusters) == 1) {
      if (!clusters %in% colnames(x[[]])) {
        stop("Cluster column '", clusters, "' not found in Seurat meta.data.")
      }
      cl <- x[[clusters]][, 1, drop = TRUE]
    } else {
      # treat as vector
      cl <- clusters
      if (length(cl) != n_cells) {
        stop("Length of 'clusters' (", length(cl), ") must match number of cells (", n_cells, ").")
      }
      # if named, align to cell order
      if (!is.null(names(cl))) {
        cn <- get_cells(x)
        if (!all(cn %in% names(cl))) {
          stop("Named 'clusters' vector must contain all cell names.")
        }
        cl <- cl[cn]
      }
    }
    as.character(cl)
  }

  cells <- get_cells(expr_or_seurat)
  clvec <- get_clusters(expr_or_seurat, clusters)

  if (drop_unused_levels) clvec <- as.character(factor(clvec))

  mask <- rep(TRUE, length(clvec))
  if (!is.null(include)) mask <- clvec %in% include
  if (!is.null(exclude)) mask <- mask & !(clvec %in% exclude)

  idx <- which(mask)
  if (length(idx) < min_cells) {
    stop("Subsetting left ", length(idx), " cells (< min_cells = ", min_cells, "). ",
         "Adjust 'include'/'exclude' or lower 'min_cells'.")
  }

  kept_labels <- clvec[idx]

  # ---- apply subset ----
  if (inherits(expr_or_seurat, "Seurat")) {
    obj_sub <- subset(expr_or_seurat, cells = cells[idx])
    if (!is.null(assay) && assay %in% names(obj_sub@assays)) {
      Seurat::DefaultAssay(obj_sub) <- assay
    }
  } else if (is_mat_like(expr_or_seurat)) {
    obj_sub <- expr_or_seurat[, idx, drop = FALSE]

  } else {
    stop("Unsupported class after checks—this should not happen.")
    }

  list(obj = obj_sub, idx = idx, cluster_labels = kept_labels)
  }

