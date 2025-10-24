


.metric_P_core <- function(module_score, pseudotime,
                           threshold = c("quantile","mean_sd"),
                           q = 0.8, alpha = 0.5,
                           soft_sigma = c("mad","sd"), k_sigma = 1.0,
                           M = 500, n_bins = 10, seed = 1L) {
  stopifnot(length(module_score) == length(pseudotime))
  set.seed(seed)
  N <- length(module_score)
  o <- order(pseudotime)
  s <- module_score[o]
  t <- pseudotime[o]

  threshold  <- match.arg(threshold)
  soft_sigma <- match.arg(soft_sigma)

  # (optional) safety: accept q in 0–1 or 1–100
  if (q > 1) q <- q/100
  if (q <= 0 || q >= 1) stop("'q' must be in (0,1).")

  if (threshold == "quantile") {
    tau <- as.numeric(quantile(s, q, na.rm = TRUE))
  } else {
    tau <- mean(s, na.rm = TRUE) + alpha*sd(s, na.rm = TRUE)
  }

  y <- ifelse(s >= tau, 1L, -1L)
  sig <- if (soft_sigma == "mad") mad(s) else sd(s)
  if (!is.finite(sig) || sig == 0) sig <- 1e-6
  z <- tanh(k_sigma * (s - tau) / sig)

  e_hard <- ifelse(y[-N] == 1 & y[-1] == 1, 1, -1)
  C_hard <- mean(e_hard)
  e_soft <- z[-N] * z[-1]
  C_soft <- mean(e_soft)

  if (any(y == 1)) var_pos <- stats::var(t[y == 1]) else var_pos <- stats::var(t) * 1e6
  var_all <- stats::var(t)
  T_var <- 1 - (var_pos / var_all)
  T_var <- max(0, min(1, T_var))

  runs <- rle(y)
  idx <- cumsum(runs$lengths)
  start_idx <- c(1, head(idx, -1) + 1)
  end_idx <- idx
  pos_runs <- which(runs$values == 1)
  if (length(pos_runs) > 0) {
    spans <- mapply(function(si, ei) t[ei] - t[si], start_idx[pos_runs], end_idx[pos_runs])
    sizes <- runs$lengths[pos_runs]
    Dbar <- sum(sizes * spans) / sum(sizes)
  } else {
    Dbar <- max(t) - min(t)
  }

  Tall <- max(t) - min(t)
  if (Tall == 0) Tall <- 1e-6
  T_run <- 1 - (Dbar / Tall)
  T_run <- max(0, min(1, T_run))
  T_peak <- (T_var + T_run)/2

  bin_id <- cut(t,
                breaks = quantile(t, probs = seq(0,1,length.out=n_bins+1),
                                  na.rm = TRUE, type = 7),
                include.lowest = TRUE, labels = FALSE)

  sim_stats <- replicate(M, {
    s_sim <- s
    for (b in seq_len(n_bins)) {
      ii <- which(bin_id == b)
      if (length(ii) > 1) s_sim[ii] <- sample(s_sim[ii], length(ii), replace = FALSE)
    }
    y_sim <- ifelse(s_sim >= tau, 1L, -1L)
    z_sim <- tanh(k_sigma * (s_sim - tau) / sig)
    C_hard_sim <- mean(ifelse(y_sim[-N] == 1 & y_sim[-1] == 1, 1, -1))
    C_soft_sim <- mean(z_sim[-N] * z_sim[-1])
    if (any(y_sim == 1)) var_pos_sim <- stats::var(t[y_sim == 1]) else var_pos_sim <- var_all * 1e6
    T_var_sim <- 1 - (var_pos_sim / var_all)
    T_var_sim <- max(0, min(1, T_var_sim))

    runs_sim <- rle(y_sim)
    idx_sim <- cumsum(runs_sim$lengths)
    start_idx_sim <- c(1, head(idx_sim, -1) + 1)
    end_idx_sim <- idx_sim
    pos_runs_sim <- which(runs_sim$values == 1)
    if (length(pos_runs_sim) > 0) {
      spans_sim <- mapply(function(si, ei) t[ei] - t[si], start_idx_sim[pos_runs_sim], end_idx_sim[pos_runs_sim])
      sizes_sim <- runs_sim$lengths[pos_runs_sim]
      Dbar_sim <- sum(sizes_sim * spans_sim) / sum(sizes_sim)
    } else {
      Dbar_sim <- Tall
    }
    T_run_sim <- 1 - (Dbar_sim / Tall)
    T_run_sim <- max(0, min(1, T_run_sim))
    T_peak_sim <- (T_var_sim + T_run_sim)/2
    c(C_hard_sim = C_hard_sim, C_soft_sim = C_soft_sim, T_peak_sim = T_peak_sim)
  })

  C_hard_null <- sim_stats["C_hard_sim",]
  C_soft_null <- sim_stats["C_soft_sim",]
  T_peak_null <- sim_stats["T_peak_sim",]

  Pc_hard <- mean(C_hard_null <= C_hard)
  Pc_soft <- mean(C_soft_null <= C_soft)
  Pt      <- mean(T_peak_null <= T_peak)

  p_C_hard <- mean(C_hard_null >= C_hard)
  p_C_soft <- mean(C_soft_null >= C_soft)
  p_T      <- mean(T_peak_null >= T_peak)

  P_final <- mean(c(Pc_soft, Pt))

  list(
    threshold = tau,
    C_hard = C_hard, C_soft = C_soft, T_var = T_var, T_run = T_run, T_peak = T_peak,
    null = list(C_hard = C_hard_null, C_soft = C_soft_null, T_peak = T_peak_null),
    percentiles = c(Pc_hard = Pc_hard, Pc_soft = Pc_soft, Pt = Pt),
    p_values = c(p_C_hard = p_C_hard, p_C_soft = p_C_soft, p_T = p_T),
    P = P_final
  )
}


#=====
.get_kegg_module_scores <- function(expr_or_seurat,
                                    pathways,
                                    naive_markers = character(0),
                                    terminal_markers = character(0),
                                    id_type = c("SYMBOL","ENSEMBL","ENTREZID"),
                                    species = "Homo sapiens",
                                    assay = NULL,
                                    slot = "data") {
  id_type <- match.arg(id_type)

  ## --- Expression as matrix ---
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

  ## --- Species → OrgDb resolver ---
  .org_pkg_for_species <- function(species) {
    sp <- tolower(trimws(species))
    aliases <- list(
      "homo sapiens"     = "org.Hs.eg.db",
      "mus musculus"     = "org.Mm.eg.db",
      "rattus norvegicus"= "org.Rn.eg.db"
    )
    aliases[[sp]]
  }

  ## --- ID mapping ---
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

  ## --- Get KEGG pathway genes ---
  .get_kegg_ids <- function(query, species, id_type) {
    tbl <- msigdbr::msigdbr(species = species,
                            category = "C2",
                            subcategory = "CP:KEGG")
    if (nrow(tbl) == 0L) stop(sprintf("No KEGG sets for '%s'.", species))
    q <- toupper(query)
    if (grepl("^\\d+$", q)) {
      code <- paste0("HSA", q)
      tbl <- tbl[tolower(tbl$gs_exact_source) == tolower(code), ]
    } else if (grepl("^[A-Z]{3}\\d+$", q)) {
      code <- sub("^[A-Z]{3}", "HSA", q)
      tbl <- tbl[tolower(tbl$gs_exact_source) == tolower(code), ]
    } else if (grepl("^KEGG_", q)) {
      tbl <- tbl[tbl$gs_name == q, ]
    } else {
      q <- gsub("[^A-Z0-9]", "_", q)
      q <- paste0("KEGG_", q)
      tbl <- tbl[grepl(q, tbl$gs_name, fixed = TRUE), ]
    }
    if (nrow(tbl) == 0L) return(character(0))
    syms <- unique(tbl$gene_symbol)
    .map_ids(syms, to = id_type, species = species)
  }

  ## --- Helpers ---
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

  ## --- main ---
  mat <- .as_matrix(expr_or_seurat)
  used_sets <- list()
  for (q in pathways) {
    ids <- .get_kegg_ids(q, species, id_type)
    ids2 <- .remove_overlaps(ids, naive_markers, terminal_markers)
    used_sets[[q]] <- ids2
  }

  mat_z <- .row_z(mat)
  scores <- .module_scores(mat_z, used_sets)
  list(scores = scores, gene_sets = used_sets)
}

# ===

metrics_P <- function(expr_or_seurat,
                      pseudotime,
                      pathways,
                      naive_markers = character(0),
                      terminal_markers = character(0),
                      id_type = c("SYMBOL","ENSEMBL","ENTREZID"),
                      species = "Homo sapiens",
                      assay = NULL,
                      slot = "data",
                      verbose = TRUE,
                      ...) {
  id_type <- match.arg(id_type)
  .msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))

  ## 1. Get KEGG module scores
  res <- .get_kegg_module_scores(expr_or_seurat = expr_or_seurat,
                                 pathways       = pathways,
                                 naive_markers  = naive_markers,
                                 terminal_markers = terminal_markers,
                                 id_type        = id_type,
                                 species        = species,
                                 assay          = assay,
                                 slot           = slot)

  module_scores <- res$scores
  gene_sets     <- res$gene_sets

  ## 2. Align pseudotime
  if (!is.null(names(pseudotime))) {
    common <- intersect(colnames(module_scores), names(pseudotime))
    pseudotime <- pseudotime[common]
    module_scores <- module_scores[, common, drop = FALSE]
  }

  ## 3. Compute .metric_P_core() for each pathway
  P_results <- list()
  for (m in rownames(module_scores)) {
    sc <- as.numeric(module_scores[m, ])
    if (all(is.na(sc))) {
      .msg("Skipping %s (no valid scores)", m)
      next
    }
    .msg("Computing P metric for %s ...", m)
    P_results[[m]] <- .metric_P_core(sc, pseudotime, ...)
  }


    ## 4) Minimal summary (only requested fields)
    summary_tbl <- tibble::tibble(
      pathway  = names(P_results),
      P_value  = vapply(P_results, function(x) x$P, numeric(1)),
      C_soft   = vapply(P_results, function(x) x$C_soft, numeric(1)),
      T_peak   = vapply(P_results, function(x) x$T_peak, numeric(1)),
      n_genes  = vapply(names(P_results),
                        function(p) length(gene_sets[[p]]), integer(1))
    )

    ## 5) Store "everything else" in a separate list
    keep_in_summary <- c("P", "C_soft", "T_peak")
    details <- lapply(P_results, function(x) {
      x[setdiff(names(x), keep_in_summary)]
    })

    ## 6) Overall P (mean of per-pathway P)
    P_overall <- mean(summary_tbl$P_value, na.rm = TRUE)

    list(
      P_overall    = P_overall,
      summary      = summary_tbl,   # only P_value, C_soft, T_peak, n_genes
      details      = details,       # all other fields per pathway
      module_scores = module_scores,
      gene_sets     = gene_sets
    )
  }










#' Plot Program Coherence (P) Metric Results
plot_P_module <- function(module_score, pseudotime, core, main = NULL, span = 0.5,
                          pch = 16, cex = 0.7) {
  stopifnot(length(module_score) == length(pseudotime))
  ord <- order(pseudotime)
  t <- as.numeric(pseudotime[ord])
  s <- as.numeric(module_score[ord])

  tau <- core$threshold
  y   <- ifelse(s >= tau, 1L, -1L)
  cols <- ifelse(y == 1L, "red", "grey70")

  plot(t, s, col = cols, pch = pch, cex = cex,
       xlab = "Pseudotime", ylab = "Module score (z)",
       main = main)

  # smooth curve (base R loess)
  fit <- try(stats::loess(s ~ t, span = span), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    tt <- seq(min(t), max(t), length.out = 200)
    lines(tt, predict(fit, newdata = data.frame(t = tt)), lwd = 2)
  }

  # threshold line
  abline(h = tau, col = "red", lty = 2)

  legend("topleft", inset = 0.01,
         legend = c("+1 (>= τ)", "-1 (< τ)"),
         col = c("red", "grey70"), pch = 16, bty = "n")
}
# plot_P_module(P$module_scores[1,], pseudotime = seu$monocle3_pseudotime, core = P$details$hsa04660, main = "hsa04660")

