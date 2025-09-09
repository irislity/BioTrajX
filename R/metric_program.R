#' Program coherence (P) via GAM(R^2) of module scores ~ s(pseudotime)
#'
#' @description
#' Computes DOPE's Program coherence (P) by:
#' 1) Assembling pathway gene sets (either from KEGG via msigdbr or user-supplied),
#' 2) Removing any overlap with naïve/terminal markers,
#' 3) Computing per-cell mean-z module scores, and
#' 4) Fitting mgcv::gam(score ~ s(pseudotime)) per module and averaging R^2.
#'
#' @param expr_or_seurat Matrix (genes x cells) **or** a Seurat object.
#' @param pseudotime Numeric vector of length = #cells; name by cell if possible.
#' @param pathways Character vector of KEGG queries (e.g., "hsa04660" or name fragments).
#'                 Ignored if \code{gene_sets} is provided.
#' @param gene_sets Optional named list of gene vectors (ID space given by \code{id_type}).
#'                  If supplied, skips KEGG lookup.
#' @param naive_markers,terminal_markers Character vectors in the same ID space as \code{id_type}.
#' @param id_type One of "SYMBOL","ENSEMBL","ENTREZID". Default "SYMBOL".
#' @param species Species label for msigdbr (default "Homo sapiens").
#' @param assay,slot If a Seurat object is provided, which assay/slot to pull (default: current assay, slot "data").
#' @param min_remaining Warn if filtered pathway < this many genes (default 10).
#' @param min_fraction Warn if filtered pathway < this fraction of original (default 0.20).
#' @param min_genes_per_module Drop modules with < this many present genes (default 3).
#' @param verbose Logical; print progress messages.
#'
#' @return A list with:
#' \item{P}{Scalar coherence score in [0,1] (mean R^2 across modules).}
#' \item{R2}{Named numeric vector of per-module R^2 (or dev.expl fallback).}
#' \item{modules_used}{Character vector of module names kept.}
#' \item{modules_dropped}{Character vector of modules dropped (< min_genes_per_module or all-NA).}
#' \item{sets_used}{Named list of genes per module actually used (post-overlap filter & present in expr).}
#' \item{summary}{A tibble summarizing modules, sizes, and overlaps (if KEGG used).}
#'
#' @export
metrics_P <- function(expr_or_seurat,
                      pseudotime,
                      pathways = NULL,
                      gene_sets = NULL,
                      naive_markers = character(0),
                      terminal_markers = character(0),
                      id_type = c("SYMBOL","ENSEMBL","ENTREZID"),
                      species = "Homo sapiens",
                      assay = NULL,
                      slot  = "data",
                      min_remaining = 10,
                      min_fraction  = 0.20,
                      min_genes_per_module = 3,
                      verbose = TRUE) {

  id_type <- match.arg(id_type)

  ## ---------- helpers (internal) ----------
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  .as_matrix <- function(x) {
    if (inherits(x, "Seurat")) {
      if (!requireNamespace("Seurat", quietly = TRUE))
        stop("Seurat not installed but a Seurat object was supplied.")
      a <- if (is.null(assay)) Seurat::DefaultAssay(x) else assay
      Seurat::GetAssayData(x, assay = a, slot = slot)
    } else {
      as.matrix(x)
    }
  }

  .map_ids <- function(symbols, to) {
    if (to == "SYMBOL") return(unique(symbols))
    if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
        !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("ID mapping requested but AnnotationDbi/org.Hs.eg.db not available.")
    }
    res <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(symbols),
                                 keytype = "SYMBOL", columns = to)
    unique(stats::na.omit(res[[to]]))
  }

  .get_kegg_symbols <- function(query, species) {
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      stop("msigdbr not installed; cannot fetch KEGG sets. Supply `gene_sets` instead.")
    }
    tbl <- msigdbr::msigdbr(species = species, category = "C2", subcategory = "CP:KEGG")
    if (grepl("^hsa\\d+$", tolower(query))) {
      # exact KEGG ID match
      tbl <- tbl[tolower(tbl$gs_exact_source) == tolower(query), , drop = FALSE]
    } else {
      # name fragment (case-insensitive)
      tbl <- tbl[grepl(toupper(query), tbl$gs_name, fixed = TRUE), , drop = FALSE]
    }
    unique(tbl$gene_symbol)
  }

  .remove_overlaps <- function(genes, naive, terminal) {
    setdiff(unique(genes), unique(c(naive, terminal)))
  }

  .z_by_gene <- function(mat) {
    mu <- Matrix::rowMeans(mat)
    sd <- sqrt(Matrix::rowMeans((mat - mu)^2))
    sd[!is.finite(sd) | sd == 0] <- 1
    (mat - mu) / sd
  }

  .module_scores <- function(mat_z, gene_sets_named) {
    mod_names <- names(gene_sets_named)
    if (is.null(mod_names)) mod_names <- paste0("module_", seq_along(gene_sets_named))
    sc <- vapply(seq_along(gene_sets_named), function(j) {
      g <- intersect(gene_sets_named[[j]], rownames(mat_z))
      if (length(g) == 0L) return(rep(NA_real_, ncol(mat_z)))
      Matrix::colMeans(mat_z[g, , drop = FALSE])
    }, FUN.VALUE = numeric(ncol(mat_z)))
    colnames(sc) <- mod_names
    t(sc) # modules x cells
  }

  .fit_gams <- function(Z_mod_by_cell, tvec) {
    if (!requireNamespace("mgcv", quietly = TRUE))
      stop("mgcv not installed; required for GAM fitting.")
    lapply(seq_len(nrow(Z_mod_by_cell)), function(i) {
      z <- Z_mod_by_cell[i, ]
      if (all(is.na(z))) return(NULL)
      mgcv::gam(z ~ s(tvec), method = "REML")
    })
  }

  .extract_R2 <- function(fits) {
    r2 <- vapply(fits, function(f) {
      if (is.null(f)) return(NA_real_)
      s <- summary(f)
      if (!is.null(s$r.sq)) s$r.sq else if (!is.null(s$dev.expl)) s$dev.expl else NA_real_
    }, numeric(1))
    r2
  }

  ## ---------- 1) get expression matrix & align ----------
  mat <- .as_matrix(expr_or_seurat)              # genes x cells
  if (is.null(colnames(mat))) colnames(mat) <- paste0("cell_", seq_len(ncol(mat)))
  # align pseudotime to colnames if named
  if (!is.null(names(pseudotime))) {
    common <- intersect(colnames(mat), names(pseudotime))
    if (length(common) == 0L)
      stop("No overlap between colnames(expr) and names(pseudotime).")
    mat <- mat[, common, drop = FALSE]
    tvec <- as.numeric(pseudotime[common])
  } else {
    if (length(pseudotime) != ncol(mat))
      stop("Length of `pseudotime` must equal ncol(expr).")
    tvec <- as.numeric(pseudotime)
  }

  ## ---------- 2) build gene sets ----------
  used_sets <- list()
  sets_source <- if (!is.null(gene_sets)) "user" else "kegg"

  if (sets_source == "user") {
    # assume gene_sets already in id_type
    used_sets <- gene_sets
  } else {
    if (is.null(pathways) || length(pathways) == 0L)
      stop("Provide either `gene_sets` or a non-empty `pathways` vector.")
    for (q in pathways) {
      syms <- .get_kegg_symbols(q, species = species)
      if (length(syms) == 0L) {
        .msg("No KEGG genes for '%s' (species=%s); skipping.", q, species)
        next
      }
      ids  <- .map_ids(syms, to = id_type)
      ids2 <- .remove_overlaps(ids, naive_markers, terminal_markers)
      # warn if filtered small
      if (length(ids) > 0) {
        frac <- length(ids2) / length(ids)
        if (length(ids2) < min_remaining || frac < min_fraction) {
          .msg("Pathway %s small after filtering: %d/%d (%.1f%%) kept.",
               q, length(ids2), length(ids), 100*frac)
        }
      }
      used_sets[[q]] <- ids2
    }
  }

  if (length(used_sets) == 0L)
    stop("No usable gene sets (all empty after filtering?).")

  ## ---------- 3) drop tiny modules & compute module scores ----------
  # keep only modules with >= min_genes_per_module present in mat
  present_counts <- vapply(used_sets, function(g) length(intersect(g, rownames(mat))), integer(1))
  kept <- names(present_counts)[present_counts >= min_genes_per_module]
  dropped <- names(present_counts)[present_counts <  min_genes_per_module]

  if (length(kept) == 0L) {
    warning("All modules dropped (< min_genes_per_module present). Returning NA.")
    return(list(
      P = NA_real_, R2 = setNames(numeric(0), character(0)),
      modules_used = character(0), modules_dropped = names(used_sets),
      sets_used = used_sets, summary = utils::stack(list(message = "All modules dropped"))
    ))
  }

  used_sets_kept <- used_sets[kept]
  mat_z <- .z_by_gene(mat)
  # modules x cells
  Z_mod_by_cell <- .module_scores(mat_z, used_sets_kept)

  ## ---------- 4) fit GAM per module ----------
  fits <- .fit_gams(Z_mod_by_cell, tvec)
  names(fits) <- rownames(Z_mod_by_cell)

  R2 <- .extract_R2(fits)
  names(R2) <- rownames(Z_mod_by_cell)

  # overall P = mean R^2
  P <- mean(R2, na.rm = TRUE)

  ## ---------- 5) summary ----------
  smry <- tibble::tibble(
    module          = names(used_sets),
    n_genes_declared= vapply(used_sets, length, integer(1)),
    n_genes_present = present_counts,
    kept            = names(used_sets) %in% kept
  )

  list(
    P = P,
    R2 = R2,
    modules_used = kept,
    modules_dropped = dropped,
    sets_used = used_sets_kept,
    summary = smry
  )
}
