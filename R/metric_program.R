#' P metric: Program coherence along pseudotime (species-aware, fixed for KEGG mapping)
#'
#' @description
#' Compute the Program (P) metric by fitting GAMs of module (gene set) z-scores
#' against pseudotime and aggregating explained variance (R^2) across modules.
#' Gene sets can come from user input (`gene_sets`) or from MSigDB KEGG via
#' `msigdbr` given `pathways`. Supports Seurat or matrix input. Species-aware
#' ID mapping via the appropriate Bioconductor OrgDb.
#'
#' @param expr_or_seurat Matrix (genes x cells) **or** a Seurat object.
#' @param pseudotime Numeric vector (length = #cells).
#' @param pathways Character vector of KEGG queries (e.g. "mmu04640").
#' @param gene_sets Named list of character vectors (modules).
#' @param naive_markers,terminal_markers Character vectors of marker IDs to exclude.
#' @param id_type One of c("SYMBOL","ENSEMBL","ENTREZID").
#' @param species Species common name (e.g., "Mus musculus").
#' @param assay,slot Seurat assay/slot.
#' @param min_remaining,min_fraction,min_genes_per_module Filtering controls.
#' @param verbose Logical.
#'
#' @return List with P metric results.
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
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  ## ---------- expression as matrix ----------
  .as_matrix <- function(x) {
    if (inherits(x, "Seurat")) {
      a <- if (is.null(assay)) Seurat::DefaultAssay(x) else assay
      Seurat::GetAssayData(x, assay = a, slot = slot)
    } else {
      m <- as.matrix(x)
      if (is.null(rownames(m))) stop("Expression matrix must have gene rownames.")
      m
    }
  }

  ## ---------- species → OrgDb resolver ----------
  .org_pkg_for_species <- function(species) {
    sp <- tolower(trimws(species))
    aliases <- list(
      "homo sapiens"     = "org.Hs.eg.db",
      "mus musculus"     = "org.Mm.eg.db",
      "rattus norvegicus"= "org.Rn.eg.db"
    )
    aliases[[sp]]
  }

  ## ---------- ID mapping helper ----------
  .map_ids <- function(symbols, to, species) {
    sy <- unique(symbols)
    if (to == "SYMBOL") return(sy)

    pkg <- .org_pkg_for_species(species)
    if (is.null(pkg)) stop(sprintf("No OrgDb mapping for %s", species))
    OrgDb <- getExportedValue(pkg, pkg)

    res <- AnnotationDbi::select(OrgDb,
                                 keys = sy,
                                 keytype = "SYMBOL",
                                 columns = to)
    unique(stats::na.omit(res[[to]]))
  }

  ## ---------- KEGG from msigdbr ----------
  .get_kegg_ids <- function(query, species, id_type) {
    tbl <- msigdbr::msigdbr(species = species,
                            category = "C2",
                            subcategory = "CP:KEGG")
    if (nrow(tbl) == 0L) {
      stop(sprintf("No KEGG sets found for species '%s'.", species))
    }

    q <- toupper(query)

    # Case 1: numeric KEGG code (e.g. "04640")
    if (grepl("^\\d+$", q)) {
      code <- paste0("HSA", q)
      tbl <- tbl[tolower(tbl$gs_exact_source) == tolower(code), , drop = FALSE]

      # Case 2: organism-prefixed KEGG code (mmu04640, hsa04640, dre04640...)
    } else if (grepl("^[A-Z]{3}\\d+$", q)) {
      code <- sub("^[A-Z]{3}", "HSA", q)  # force to HSA#####
      tbl <- tbl[tolower(tbl$gs_exact_source) == tolower(code), , drop = FALSE]

      # Case 3: KEGG_ style (KEGG_HEMATOPOIETIC_CELL_LINEAGE)
    } else if (grepl("^KEGG_", q)) {
      tbl <- tbl[tbl$gs_name == q, , drop = FALSE]

      # Case 4: Plain text pathway name (e.g. "Hematopoietic cell lineage")
    } else {
      q <- gsub("[^A-Z0-9]", "_", q)   # normalize spaces → underscores
      q <- paste0("KEGG_", q)
      tbl <- tbl[grepl(q, tbl$gs_name, fixed = TRUE), , drop = FALSE]
    }

    if (nrow(tbl) == 0L) {
      return(character(0))
    }

    syms <- unique(tbl$gene_symbol)
    .map_ids(syms, to = id_type, species = species)
  }


  ## ---------- utilities ----------
  .remove_overlaps <- function(genes, naive, terminal) {
    setdiff(unique(genes), unique(c(naive, terminal)))
  }
  .row_z <- function(mat) {
    mu <- rowMeans(mat)
    sd <- sqrt(rowMeans((mat - mu)^2))
    sd[!is.finite(sd) | sd == 0] <- 1
    sweep(mat, 1, mu, "-") / sd
  }
  .module_scores <- function(mat_z, gene_sets_named) {
    mods <- names(gene_sets_named)
    sc <- vapply(seq_along(gene_sets_named), function(j) {
      g <- intersect(gene_sets_named[[j]], rownames(mat_z))
      if (length(g) == 0L) return(rep(NA_real_, ncol(mat_z)))
      colMeans(mat_z[g, , drop = FALSE])
    }, FUN.VALUE = numeric(ncol(mat_z)))
    colnames(sc) <- mods
    t(sc)
  }
  .fit_gams <- function(Z_mod_by_cell, tvec) {
    lapply(seq_len(nrow(Z_mod_by_cell)), function(i) {
      z <- Z_mod_by_cell[i, ]
      if (all(is.na(z))) return(NULL)
      mgcv::gam(z ~ s(tvec), method = "REML")
    })
  }
  .extract_R2 <- function(fits) {
    vapply(fits, function(f) {
      if (is.null(f)) return(NA_real_)
      s <- summary(f)
      if (!is.null(s$r.sq)) s$r.sq else if (!is.null(s$dev.expl)) s$dev.expl else NA_real_
    }, numeric(1))
  }

  ## ---------- 1) load expression ----------
  mat <- .as_matrix(expr_or_seurat)
  if (is.null(colnames(mat))) colnames(mat) <- paste0("cell_", seq_len(ncol(mat)))

  if (!is.null(names(pseudotime))) {
    common <- intersect(colnames(mat), names(pseudotime))
    mat <- mat[, common, drop = FALSE]
    tvec <- as.numeric(pseudotime[common])
  } else {
    tvec <- as.numeric(pseudotime)
  }

  ## ---------- 2) build gene sets ----------
  used_sets <- list()
  if (!is.null(gene_sets)) {
    used_sets <- gene_sets
  } else {
    for (q in pathways) {
      ids <- .get_kegg_ids(q, species, id_type)
      ids2 <- .remove_overlaps(ids, naive_markers, terminal_markers)
      used_sets[[q]] <- ids2
      if (length(ids2) == 0L) .msg("No KEGG genes for '%s' (species=%s)", q, species)
    }
  }
  if (length(used_sets) == 0L) stop("No usable gene sets (all empty).")

  ## ---------- 3) drop tiny modules ----------
  present_counts <- vapply(used_sets, function(g)
    length(intersect(g, rownames(mat))), integer(1))
  kept <- names(present_counts)[present_counts >= min_genes_per_module]
  used_sets_kept <- used_sets[kept]

  ## ---------- 4) scores + GAM ----------
  mat_z <- .row_z(mat)
  Z_mod_by_cell <- .module_scores(mat_z, used_sets_kept)
  fits <- .fit_gams(Z_mod_by_cell, tvec)
  R2 <- .extract_R2(fits)
  names(R2) <- rownames(Z_mod_by_cell)

  P <- mean(R2, na.rm = TRUE)

  smry <- tibble::tibble(
    module = names(used_sets),
    n_genes_declared = vapply(used_sets, length, integer(1)),
    n_genes_present  = present_counts,
    kept             = names(used_sets) %in% kept
  )

  list(P = P, R2 = R2, modules_used = kept,
       modules_dropped = setdiff(names(used_sets), kept),
       sets_used = used_sets_kept, summary = smry)
}
