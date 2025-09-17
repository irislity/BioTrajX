# Branched DOPE

# Utilities ----
`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x

.stop_if_missing <- function(lst, required_names, label) {
  miss <- setdiff(required_names, names(lst))
  if (length(miss)) stop("Missing names in ", label, ": ", paste(miss, collapse = ", "))
}

# Expect a named list like:
# naive_markers_list   = list(branchA = c("..."), branchB = c("..."), ...)
# terminal_markers_list= list(branchA = c("..."), branchB = c("..."), ...)
# pathways_list        = list(branchA = c("..."), branchB = c("..."), ...) or NULL if using gene_sets_list
# gene_sets_list       = list(branchA = list(mod1=c(...), ...), branchB=..., ...) or NULL
# branch_filters       = named list controlling which cells belong to the branch (optional):
#   branch_filters = list(
#     branchA = list(include = c("Stem","ProgA"), exclude = NULL, min_cells = 20),
#     branchB = list(include = c("ProgB","TermB"))
#   )
#
# If branch_filters[[b]] is NULL, all cells are used (no extra subsetting).
# If you pass cluster_labels (meta column name or vector), subsetting is by clusters.

# ------- Core: single trajectory across multiple branches -------
#' Compute DOPE for one pseudotime with multiple branches
#' @param expr_or_seurat matrix/Seurat
#' @param pseudotime numeric vector (named by cell or same length as ncol(expr))
#' @param naive_markers_list named list: branch -> character vector (markers)
#' @param terminal_markers_list named list: branch -> character vector (markers)
#' @param pathways_list named list: branch -> character vector (KEGG queries) OR NULL when using gene_sets_list
#' @param gene_sets_list named list: branch -> named list of modules (ID space = id_type) OR NULL to use KEGG
#' @param cluster_labels Seurat meta column name OR per-cell vector (required for branch subsetting)
#' @param branch_filters named list: branch -> list(include=..., exclude=..., min_cells=..., drop_unused_levels=TRUE)
#' @param E_method,id_type,species,assay,slot,min_remaining,min_fraction,min_genes_per_module,tol,plot_E,verbose see previous
#' @return list with per-branch DOPE + aggregate
compute_DOPE_single_branched <- function(expr_or_seurat,
                                         pseudotime,
                                         naive_markers_list,
                                         terminal_markers_list,
                                         pathways_list = NULL,
                                         gene_sets_list = NULL,
                                         cluster_labels = NULL,
                                         branch_filters = NULL,
                                         E_method = c("gmm","clusters","combined"),
                                         id_type = c("SYMBOL","ENSEMBL","ENTREZID"),
                                         species = "Homo sapiens",
                                         assay = NULL,
                                         slot = "data",
                                         min_remaining = 10,
                                         min_fraction  = 0.20,
                                         min_genes_per_module = 3,
                                         tol = 1e-8,
                                         plot_E = FALSE,
                                         verbose = TRUE) {
  E_method <- match.arg(E_method)
  id_type  <- match.arg(id_type)

  # sanity on branches
  branch_names <- union(names(naive_markers_list), names(terminal_markers_list))
  if (is.null(branch_names) || any(branch_names == "")) stop("Branch lists must be named.")
  .stop_if_missing(naive_markers_list,   branch_names, "naive_markers_list")   # ensures all present
  .stop_if_missing(terminal_markers_list,branch_names, "terminal_markers_list")
  if (!is.null(pathways_list))  .stop_if_missing(pathways_list,  branch_names, "pathways_list")
  if (!is.null(gene_sets_list)) .stop_if_missing(gene_sets_list, branch_names, "gene_sets_list")

  # expression & cell order
  expr <- if (inherits(expr_or_seurat, "Seurat")) {
    if (is.null(assay)) assay <- Seurat::DefaultAssay(expr_or_seurat)
    as.matrix(Seurat::GetAssayData(expr_or_seurat, assay = assay, slot = slot))
  } else {
    as.matrix(expr_or_seurat)
  }
  cells_all <- if (inherits(expr_or_seurat, "Seurat")) Seurat::Cells(expr_or_seurat) else colnames(expr)

  # align pseudotime & clusters to cell order
  .align_to_cells <- function(x, target_cells, what="vector"){
    if (!is.null(names(x))) {
      if (!all(target_cells %in% names(x))) stop("Named ", what, " missing some cells.")
      return(x[target_cells])
    } else {
      if (length(x) != length(target_cells))
        stop("Length of ", what, " (", length(x), ") must match #cells (", length(target_cells), ").")
      return(x)
    }
  }
  pt_all <- .align_to_cells(pseudotime, cells_all, "pseudotime")

  if (!is.null(cluster_labels)) {
    if (inherits(expr_or_seurat, "Seurat") && is.character(cluster_labels) && length(cluster_labels) == 1) {
      if (!cluster_labels %in% colnames(expr_or_seurat[[]]))
        stop("Cluster column '", cluster_labels, "' not found in Seurat meta.data.")
      cl_all <- expr_or_seurat[[cluster_labels]][,1,drop=TRUE]
      names(cl_all) <- cells_all
    } else {
      cl_all <- .align_to_cells(cluster_labels, cells_all, "cluster_labels")
    }
  } else {
    cl_all <- NULL
  }

  # branch subsetting helper (reuses your subset_by_clusters idea but inline & minimal)
  subset_for_branch <- function(include=NULL, exclude=NULL, min_cells=10, drop_unused_levels=TRUE) {
    mask <- rep(TRUE, length(cells_all))
    if (!is.null(cl_all)) {
      if (!is.null(include)) mask <- cl_all %in% include
      if (!is.null(exclude)) mask <- mask & !(cl_all %in% exclude)
    }
    idx <- which(mask)
    if (length(idx) < min_cells) stop("Branch subsetting left ", length(idx), " cells (< ", min_cells, ").")
    list(
      expr = expr[, idx, drop=FALSE],
      pt   = pt_all[idx],
      cl   = if (is.null(cl_all)) NULL else as.character(if (drop_unused_levels) factor(cl_all[idx]) else cl_all[idx]),
      cells= cells_all[idx]
    )
  }

  # per-branch compute using your existing metrics_* functions
  run_branch <- function(b) {
    if (verbose) message(">>> Branch: ", b)
    filt <- branch_filters[[b]] %||% list()
    include <- filt$include %||% NULL
    exclude <- filt$exclude %||% NULL
    minc    <- filt$min_cells %||% 10
    drop_lv <- filt$drop_unused_levels %||% TRUE

    sub <- subset_for_branch(include, exclude, minc, drop_lv)

    # pull branch-specific definitions
    naive_markers  <- naive_markers_list[[b]]
    terminal_markers <- terminal_markers_list[[b]]

    gene_sets <- if (!is.null(gene_sets_list)) gene_sets_list[[b]] else NULL
    pathways  <- if (!is.null(pathways_list))  pathways_list[[b]]  else NULL

    # D
    D_res <- tryCatch({
      metrics_D(sub$expr, naive_markers, terminal_markers, sub$pt)
    }, error = function(e) list(D_naive=NA_real_, D_term=NA_real_, error=e$message))

    # O
    O_res <- tryCatch({
      metrics_O(sub$expr, naive_markers, terminal_markers, sub$pt, tol = tol)
    }, error = function(e) list(O=NA_real_, error=e$message))

    # P
    P_res <- tryCatch({
      metrics_P(sub$expr,
                sub$pt,
                pathways = pathways,
                gene_sets = gene_sets,
                naive_markers = naive_markers,
                terminal_markers = terminal_markers,
                id_type = id_type,
                species = species,
                assay = NULL, slot = "data",
                min_remaining = min_remaining,
                min_fraction = min_fraction,
                min_genes_per_module = min_genes_per_module,
                verbose = verbose)
    }, error = function(e) list(P=NA_real_, error=e$message))

    # E
    naive_scores <- colMeans(sub$expr[intersect(naive_markers, rownames(sub$expr)), , drop=FALSE])
    term_scores  <- colMeans(sub$expr[intersect(terminal_markers, rownames(sub$expr)), , drop=FALSE])
    E_res <- tryCatch({
      metrics_E(sub$pt,
                cluster_labels = sub$cl,
                naive_marker_scores = naive_scores,
                terminal_marker_scores = term_scores,
                naive_clusters = NULL, terminal_clusters = NULL,   # optionally supply per-branch cluster anchors if you want
                method = E_method, plot = FALSE)
    }, error = function(e) list(E_naive=NA_real_, E_term=NA_real_, E_comp=NA_real_, error=e$message))

    comp_vec <- c(D_res$D_naive %||% NA_real_,
                  D_res$D_term  %||% NA_real_,
                  O_res$O       %||% NA_real_,
                  P_res$P       %||% NA_real_,
                  E_res$E_comp  %||% NA_real_)
    score <- if (all(is.na(comp_vec))) NA_real_ else mean(comp_vec, na.rm = TRUE)

    list(
      branch = b,
      n_cells = ncol(sub$expr),
      D = D_res, O = O_res, P = P_res, E = E_res,
      DOPE_score = score
    )
  }

  branch_results <- lapply(branch_names, run_branch)
  names(branch_results) <- branch_names

  # aggregate score across branches (weighted by cells by default)
  w <- vapply(branch_results, function(x) x$n_cells, numeric(1))
  s <- vapply(branch_results, function(x) x$DOPE_score, numeric(1))
  agg <- if (all(is.na(s))) NA_real_ else stats::weighted.mean(s, w, na.rm = TRUE)

  list(
    branches = branch_results,
    aggregate_DOPE = agg,
    weights = w,
    params = list(E_method=E_method, id_type=id_type, species=species,
                  min_remaining=min_remaining, min_fraction=min_fraction,
                  min_genes_per_module=min_genes_per_module, tol=tol)
  )
}

# ------- Multi-trajectory (branched) -------
#' Compute DOPE across multiple pseudotime trajectories with branches
#' @param pseudotime_list named list: trajectory -> pseudotime vector
#' @return list with per-trajectory results and comparison table
compute_multi_DOPE_branched <- function(expr_or_seurat,
                                        pseudotime_list,
                                        naive_markers_list,
                                        terminal_markers_list,
                                        pathways_list = NULL,
                                        gene_sets_list = NULL,
                                        cluster_labels = NULL,
                                        branch_filters = NULL,
                                        E_method = c("gmm","clusters","combined"),
                                        id_type = c("SYMBOL","ENSEMBL","ENTREZID"),
                                        species = "Homo sapiens",
                                        assay = NULL,
                                        slot  = "data",
                                        min_remaining = 10,
                                        min_fraction  = 0.20,
                                        min_genes_per_module = 3,
                                        tol = 1e-8,
                                        plot_E = FALSE,
                                        verbose = TRUE,
                                        parallel = FALSE,
                                        n_cores = NULL) {
  E_method <- match.arg(E_method)
  id_type  <- match.arg(id_type)

  if (!is.list(pseudotime_list)) stop("pseudotime_list must be a named list.")
  if (is.null(names(pseudotime_list))) names(pseudotime_list) <- paste0("trajectory_", seq_along(pseudotime_list))

  run_one <- function(name, pt) {
    if (verbose) message("\n=== Trajectory: ", name, " ===")
    out <- compute_DOPE_single_branched(
      expr_or_seurat = expr_or_seurat,
      pseudotime = pt,
      naive_markers_list = naive_markers_list,
      terminal_markers_list = terminal_markers_list,
      pathways_list = pathways_list,
      gene_sets_list = gene_sets_list,
      cluster_labels = cluster_labels,
      branch_filters = branch_filters,
      E_method = E_method,
      id_type = id_type,
      species = species,
      assay = assay,
      slot = slot,
      min_remaining = min_remaining,
      min_fraction = min_fraction,
      min_genes_per_module = min_genes_per_module,
      tol = tol,
      plot_E = plot_E,
      verbose = verbose
    )
    out$trajectory <- name
    out
  }

  if (parallel && length(pseudotime_list) > 1 && requireNamespace("parallel", quietly = TRUE)) {
    if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    res_list <- parallel::parLapply(cl, seq_along(pseudotime_list), function(i) {
      name <- names(pseudotime_list)[i]; pt <- pseudotime_list[[i]]
      run_one(name, pt)
    })
    names(res_list) <- names(pseudotime_list)
  } else {
    res_list <- Map(run_one, names(pseudotime_list), pseudotime_list)
    names(res_list) <- names(pseudotime_list)
  }

  # comparison table
  comp <- data.frame(
    trajectory = names(res_list),
    aggregate_DOPE = vapply(res_list, function(x) x$aggregate_DOPE, numeric(1)),
    stringsAsFactors = FALSE
  )
  comp <- comp[order(comp$aggregate_DOPE, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  best <- if (all(is.na(comp$aggregate_DOPE))) NA_character_ else comp$trajectory[1]

  list(
    results = res_list,
    comparison = comp,
    best_trajectory = best,
    args = list(E_method=E_method, id_type=id_type, species=species, parallel=parallel, n_cores=n_cores)
  )
}

# =========================
# Example usage
# =========================
# branches <- c("Stem_to_Ery","Stem_to_Meg")
# naive_markers_list <- list(
#   Stem_to_Ery = c("CD34","KIT","GATA2"),
#   Stem_to_Meg = c("CD34","KIT","GATA2")
# )
# terminal_markers_list <- list(
#   Stem_to_Ery = c("GYPA","HBB","HBA1","TFRC","KLF1","ALAS2","EPOR"),
#   Stem_to_Meg = c("ITGA2B","GP9","PF4") # example megakaryocyte markers
# )
# pathways_list <- list(
#   Stem_to_Ery = c("KEGG_HEMATOPOIETIC_CELL_LINEAGE","KEGG_JAK_STAT_SIGNALING_PATHWAY"),
#   Stem_to_Meg = c("KEGG_HEMATOPOIETIC_CELL_LINEAGE")
# )
# branch_filters <- list(
#   Stem_to_Ery = list(include = c("Stem_Progenitors","Erythroid_progenitors_Erythroblasts"), min_cells = 50),
#   Stem_to_Meg = list(include = c("Stem_Progenitors","Megakaryocyte_progenitors"), min_cells = 50)
# )
# res_multi <- compute_multi_DOPE_branched(
#   expr_or_seurat = seu,
#   pseudotime_list = list(mono = seu$pt_mono, sling = seu$pt_sling),
#   naive_markers_list = naive_markers_list,
#   terminal_markers_list = terminal_markers_list,
#   pathways_list = pathways_list,
#   cluster_labels = "Phenotype",
#   branch_filters = branch_filters,
#   id_type = "SYMBOL",
#   species = "Mus musculus",
#   verbose = TRUE
# )
