# =============================================================================
# run_ti_methods.R
#
# Run eight trajectory inference (TI) methods (Slingshot, CytoTRACE, DPT,
# SCORPIUS, TSCAN, Monocle3, PAGA-DPT, Palantir) and return a 1D
# pseudotime vector (named numeric [0, 1]) per method.
#
# R methods  (loaded via CRAN / Bioconductor):
#   Slingshot · CytoTRACE · DPT (destiny) · SCORPIUS · TSCAN · Monocle3
#
# Python methods (loaded via reticulate — see SETUP below):
#   PAGA-DPT · Palantir
#
# INPUT
#   expr        genes × cells log-normalised expression matrix
#               (rownames = genes, colnames = cells)
#   start_cell  name of one root cell (used by directed methods)
#   n_pcs       number of PCs to compute for dimension reduction
#
# OUTPUT
#   run_all_ti_methods() returns a data.frame:
#     columns = methods, rows = cells, values = pseudotime [0, 1]
#     NA where a method failed or requires unavailable dependencies
#
# PYTHON SETUP
#   conda create -n ti python=3.10
#   conda activate ti
#   pip install scanpy palantir
#   # In R: reticulate::use_condaenv("ti")
#
# USAGE
#   source("manuscript/run_ti_methods.R")
#   f   <- make_linear_fixture(n_cells = 300)          # BioTrajX fixture
#   res <- run_all_ti_methods(f$expr, start_cell = names(f$pseudotime)[1])
#   cor(res, f$pseudotime, use = "pairwise.complete.obs", method = "spearman")
# =============================================================================

# ── Threading: must be set before any BLAS/OMP library loads ─────────────────
# Prevents segfault when Python (via reticulate) and R share an OMP runtime.
Sys.setenv(OMP_NUM_THREADS       = "1",
           OPENBLAS_NUM_THREADS  = "1",
           MKL_NUM_THREADS       = "1",
           VECLIB_MAXIMUM_THREADS= "1")

# ── Repo root & fixture loader ────────────────────────────────────────────────
.repo_root <- tryCatch({
  src <- rstudioapi::getSourceEditorContext()$path
  normalizePath(file.path(dirname(src), ".."))
}, error = function(e) getwd())

source(file.path(.repo_root, "manuscript", "scripts", "synthetic", "helper-fixtures.R"))

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Rescale a numeric vector to [0, 1]; returns NA vector on failure.
.rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(setNames(rep(NA_real_, length(x)), names(x)))
  (x - rng[1]) / diff(rng)
}

#' PCA on a genes × cells matrix; returns cells × PCs matrix.
#' Drops zero/low-variance genes and scales to avoid singular matrices.
.pca <- function(expr, n_pcs = 30) {
  gene_sd  <- apply(expr, 1, sd)
  expr     <- expr[gene_sd > 1e-6, , drop = FALSE]
  n_pcs    <- min(n_pcs, nrow(expr) - 1L, ncol(expr) - 1L)
  prcomp(t(expr), center = TRUE, scale. = TRUE, rank. = n_pcs)$x
}

#' Return NA pseudotime with a message when a package is missing.
.missing_pkg <- function(pkg, cells) {
  message(sprintf("  [skip] %s not installed — install with: %s",
                  pkg, .install_hint(pkg)))
  setNames(rep(NA_real_, length(cells)), cells)
}

.install_hint <- function(pkg) {
  bioc <- c("slingshot", "SingleCellExperiment", "destiny", "BiocGenerics", "TSCAN")
  if (pkg %in% bioc)
    sprintf('BiocManager::install("%s")', pkg)
  else if (pkg == "CytoTRACE")
    'devtools::install_github("gunsagargulati/CytoTRACE")'
  else if (pkg == "SCORPIUS")
    'remotes::install_github("rcannood/SCORPIUS")'
  else
    sprintf('install.packages("%s")', pkg)
}

# =============================================================================
# R-NATIVE METHODS
# =============================================================================

# ── 1. Slingshot ──────────────────────────────────────────────────────────────
#' Slingshot pseudotime (Bioconductor)
#'
#' Fits a minimum spanning tree and principal curves through a PCA embedding.
#' For 1D output we take the pseudotime of the first (or only) lineage and
#' average across lineages for cells assigned to multiple paths.
#'
#' @param expr        genes × cells log-normalised matrix
#' @param start_cell  name of the root cell; its PCA cluster becomes start.clus
#' @param n_pcs       PCs to use for the embedding (default 20)
#' @param precomp_pca optional cells × PCs matrix to use instead of recomputing
#' @return Named numeric [0, 1] per cell; NA for cells with no assignment
run_slingshot <- function(expr, start_cell = NULL, n_pcs = 20, precomp_pca = NULL) {
  if (!requireNamespace("slingshot", quietly = TRUE))
    return(.missing_pkg("slingshot", colnames(expr)))
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
    return(.missing_pkg("SingleCellExperiment", colnames(expr)))

  message("  [Slingshot] computing PCA + fitting curves ...")
  # Use at most 10 PCs — more causes near-singular covariance in curve fitting
  pca_raw <- if (!is.null(precomp_pca)) precomp_pca else .pca(expr, n_pcs)
  pca     <- pca_raw[, seq_len(min(10L, ncol(pca_raw))), drop = FALSE]
  sce    <- SingleCellExperiment::SingleCellExperiment(
    assays  = list(logcounts = expr),
    reducedDims = list(PCA = pca)
  )

  # k-means clusters for Slingshot's MST — cap at 10 to avoid singularity
  k       <- min(10L, max(3L, round(sqrt(ncol(expr)) / 3)))
  cl      <- stats::kmeans(pca, centers = k, nstart = 25, iter.max = 50)$cluster
  start_c <- if (!is.null(start_cell)) as.character(cl[start_cell]) else NULL

  sds <- tryCatch(
    slingshot::slingshot(sce, clusterLabels = cl,
                         reducedDim = "PCA", start.clus = start_c),
    error = function(e) { message("    Slingshot error: ", e$message); NULL }
  )
  if (is.null(sds)) return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))

  # slingPseudotime: cells × lineages; slingCurveWeights: cells × lineages.
  # Weighted mean collapses multiple lineages to a single 1D pseudotime,
  # giving more weight to the lineage each cell is most confidently assigned to.
  pt_mat  <- slingshot::slingPseudotime(sds)
  wt_mat  <- slingshot::slingCurveWeights(sds)
  # Mask weights where pseudotime is NA, then row-normalise
  wt_mat[is.na(pt_mat)] <- 0
  pt_mat[is.na(pt_mat)] <- 0
  row_sums <- rowSums(wt_mat)
  # Cells with zero total weight (unassigned) get NA
  pt <- ifelse(row_sums == 0, NA_real_,
               rowSums(pt_mat * wt_mat) / row_sums)
  .rescale01(setNames(pt, colnames(expr)))
}

# ── 2. CytoTRACE ──────────────────────────────────────────────────────────────
#' CytoTRACE v1 developmental potential score
#'
#' Scores cells by gene count structure; higher score = more primitive (naive).
#' We invert so 0 = most naive → 1 = most differentiated (terminal).
#'
#' CytoTRACE is sensitive to the number of genes detected per cell, so passing
#' the full gene matrix (not an HVG subset) is recommended.  Supply
#' fullgene_expr (all genes × cells, log-normalised) via run_all_ti_methods();
#' falls back to expr (which may be HVG-subsetted) if not provided.
#'
#' Install: devtools::install_github("gunsagargulati/CytoTRACE")
#'
#' @param expr          genes × cells log-normalised matrix (fallback)
#' @param fullgene_expr all-gene log-normalised matrix; used in preference to expr
#' @return Named numeric [0, 1] per cell
run_cytotrace <- function(expr, fullgene_expr = NULL) {
  if (!requireNamespace("CytoTRACE", quietly = TRUE))
    return(.missing_pkg("CytoTRACE", colnames(expr)))

  mat <- if (!is.null(fullgene_expr)) {
    message("  [CytoTRACE] using full-gene matrix (", nrow(fullgene_expr), " genes) ...")
    fullgene_expr
  } else {
    message("  [CytoTRACE] computing developmental potential ...")
    expr
  }

  res <- tryCatch(
    CytoTRACE::CytoTRACE(mat, ncores = 1L),
    error = function(e) { message("    CytoTRACE error: ", e$message); NULL }
  )
  if (is.null(res)) return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))

  # CytoTRACE score: 1 = most primitive; invert so 1 = most differentiated
  raw <- res$CytoTRACE[colnames(expr)]
  .rescale01(1 - raw)
}

# ── 3. Monocle 3 ──────────────────────────────────────────────────────────────
#' Monocle 3 pseudotime via principal graph
#'
#' Builds a cell_data_set, runs UMAP, learns a principal graph, and roots it
#' at the partition containing start_cell. For a single connected trajectory
#' this produces a well-defined 1D pseudotime.
#'
#' @param expr       genes × cells log-normalised matrix
#' @param start_cell name of the root cell
#' @return Named numeric [0, 1] per cell; Inf cells mapped to NA
run_monocle3 <- function(expr, start_cell = NULL) {
  for (pkg in c("monocle3", "SingleCellExperiment", "Matrix")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      return(.missing_pkg(pkg, colnames(expr)))
  }

  message("  [Monocle3] building CDS, learning graph ...")
  gene_meta <- data.frame(gene_short_name = rownames(expr),
                          row.names       = rownames(expr))
  cell_meta <- data.frame(cell_id = colnames(expr),
                          row.names = colnames(expr))

  cds <- tryCatch({
    monocle3::new_cell_data_set(
      Matrix::Matrix(expr, sparse = TRUE),
      cell_metadata = cell_meta,
      gene_metadata = gene_meta
    )
  }, error = function(e) { message("    CDS error: ", e$message); NULL })
  if (is.null(cds)) return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))

  cds <- monocle3::preprocess_cds(cds, num_dim = min(30L, ncol(expr) - 1L))
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
  cds <- monocle3::cluster_cells(cds)
  cds <- monocle3::learn_graph(cds, use_partition = FALSE)

  # Root at the partition/node nearest to start_cell
  if (!is.null(start_cell) && start_cell %in% colnames(cds)) {
    cds <- monocle3::order_cells(cds, root_cells = start_cell)
  } else {
    # Auto-root at the cell with the lowest UMAP1 coordinate
    umap  <- SingleCellExperiment::reducedDim(cds, "UMAP")
    root  <- rownames(umap)[which.min(umap[, 1])]
    cds   <- monocle3::order_cells(cds, root_cells = root)
  }

  pt  <- monocle3::pseudotime(cds)
  pt[is.infinite(pt)] <- NA_real_
  .rescale01(pt)
}

# ── 4. DPT (destiny) ─────────────────────────────────────────────────────────
#' Diffusion Pseudotime via the destiny package
#'
#' Constructs a diffusion map on the expression space and computes DPT from
#' the root cell. DPT is inherently 1D (a single pseudotime axis).
#'
#' @param expr       genes × cells log-normalised matrix
#' @param start_cell name of the root cell (if NULL, uses cell 1)
#' @param n_pcs      number of PCs to use as input to DiffusionMap
#' @return Named numeric [0, 1] per cell
run_dpt <- function(expr, start_cell = NULL, n_pcs = 30, precomp_pca = NULL) {
  if (!requireNamespace("destiny", quietly = TRUE))
    return(.missing_pkg("destiny", colnames(expr)))

  message("  [DPT] computing diffusion map ...")
  pca <- if (!is.null(precomp_pca)) precomp_pca else .pca(expr, n_pcs)
  root_idx <- if (!is.null(start_cell) && start_cell %in% rownames(pca))
    which(rownames(pca) == start_cell) else 1L

  dm <- tryCatch(
    destiny::DiffusionMap(pca, verbose = FALSE),
    error = function(e) { message("    DiffusionMap error: ", e$message); NULL }
  )
  if (is.null(dm)) return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))

  dpt_obj <- tryCatch(
    destiny::DPT(dm, tips = root_idx),
    error = function(e) { message("    DPT error: ", e$message); NULL }
  )
  if (is.null(dpt_obj)) return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))

  pt <- dpt_obj$dpt
  .rescale01(setNames(pt, colnames(expr)))
}

# ── 5. SCORPIUS ───────────────────────────────────────────────────────────────
# ── SCORPIUS ───────────────────────────────────────────────────────────────────
#' SCORPIUS pseudotime (linear trajectories)
#'
#' SCORPIUS is designed for 1D linear trajectories and returns a pseudotime
#' vector in [0, 1] directly. It works on a cells × genes matrix.
#'
#' SCORPIUS's principal-curve fit has no native root-cell argument, so its
#' output orientation is arbitrary. If start_cell is supplied, the trajectory
#' is flipped (1 - time) when needed so that start_cell falls in the earlier
#' half, anchoring it to the same root used by the other TI methods.
#'
#' @param expr genes × cells log-normalised matrix
#' @param start_cell name of the root cell (NULL = leave orientation as-is)
#' @return Named numeric [0, 1] per cell
run_scorpius <- function(expr, start_cell = NULL) {
  if (!requireNamespace("SCORPIUS", quietly = TRUE))
    return(.missing_pkg("SCORPIUS", colnames(expr)))

  message("  [SCORPIUS] inferring 1D trajectory ...")
  # SCORPIUS expects cells × genes
  expr_T  <- t(expr)
  space   <- tryCatch(
    SCORPIUS::reduce_dimensionality(expr_T, dist = "spearman"),
    error = function(e) { message("    SCORPIUS reduce error: ", e$message); NULL }
  )
  if (is.null(space)) return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))

  traj <- tryCatch(
    SCORPIUS::infer_trajectory(space),
    error = function(e) { message("    SCORPIUS infer error: ", e$message); NULL }
  )
  if (is.null(traj)) return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))

  pt <- .rescale01(setNames(traj$time, colnames(expr)))

  if (!is.null(start_cell) && start_cell %in% names(pt) && pt[[start_cell]] > 0.5) {
    message(sprintf("  [SCORPIUS] flipping orientation to anchor root cell %s", start_cell))
    pt <- 1 - pt
  }

  pt
}

# ── 6. TSCAN ──────────────────────────────────────────────────────────────────
#' TSCAN pseudotime (MST-based ordering)
#'
#' TSCAN clusters cells, fits a minimum spanning tree over cluster centers,
#' and orders cells along the longest path through the tree. Like SCORPIUS,
#' it has no native root-cell argument, so cells are ranked by their position
#' along the ordered path and the result is flipped (1 - time) when needed
#' so start_cell falls in the earlier half, anchoring it to the same root
#' used by the other TI methods.
#'
#' Cells that fall on a side branch excluded from the main MST path (rather
#' than the single ordered path TSCANorder() returns) are set to NA.
#'
#' @param expr genes × cells log-normalised matrix
#' @param start_cell name of the root cell (NULL = leave orientation as-is)
#' @return Named numeric [0, 1] per cell; NA for cells off the main path
run_tscan <- function(expr, start_cell = NULL) {
  if (!requireNamespace("TSCAN", quietly = TRUE))
    return(.missing_pkg("TSCAN", colnames(expr)))

  message("  [TSCAN] clustering + MST ordering ...")
  mc <- tryCatch(
    TSCAN::exprmclust(expr, clusternum = 2:9, reduce = TRUE),
    error = function(e) { message("    TSCAN exprmclust error: ", e$message); NULL }
  )
  if (is.null(mc)) return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))

  ord <- tryCatch(
    TSCAN::TSCANorder(mc, orderonly = TRUE),
    error = function(e) { message("    TSCAN order error: ", e$message); NULL }
  )
  if (is.null(ord) || length(ord) == 0)
    return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))

  # TSCANorder() returns only cells on the main path; off-path cells get NA
  pt <- setNames(rep(NA_real_, ncol(expr)), colnames(expr))
  pt[ord] <- seq_along(ord)
  pt <- .rescale01(pt)

  if (!is.null(start_cell) && start_cell %in% names(pt) &&
      !is.na(pt[[start_cell]]) && pt[[start_cell]] > 0.5) {
    message(sprintf("  [TSCAN] flipping orientation to anchor root cell %s", start_cell))
    pt <- 1 - pt
  }

  pt
}


# =============================================================================
# PYTHON METHODS  (via reticulate)
# =============================================================================
# These wrappers require a Python environment with the relevant packages.
# Set .python_env below or export RETICULATE_PYTHON before sourcing this file.
# If the environment is unavailable the method returns NA with a message.
#
# "sc_tutorial" (a large, unpinned tutorial env) had multiple bundled
# libomp.dylib copies (sklearn, torch) colliding at runtime and segfaulting
# scanpy's sc.pp.neighbors() / leiden on this machine. "biotrajx_ti" is a
# minimal conda-forge env (python-igraph, leidenalg, scanpy, pandas,
# scikit-learn, palantir — see environment.yml at the repo root) that avoids
# the conflict.
.python_env <- Sys.getenv("RETICULATE_CONDAENV", unset = "biotrajx_ti")

.ensure_reticulate <- function(condaenv = .python_env) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    message('  [Python] reticulate not installed — install.packages("reticulate")')
    return(FALSE)
  }
  # Set threading env vars in R before any OMP libraries initialise in Python
  Sys.setenv(OMP_NUM_THREADS      = "1",
             OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS      = "1",
             VECLIB_MAXIMUM_THREADS = "1")
  tryCatch(
    reticulate::use_condaenv(condaenv, required = TRUE),
    error = function(e) {
      message(sprintf("  [Python] conda env '%s' not found: %s", condaenv, e$message))
      return(FALSE)
    }
  )
  TRUE
}

# ── Python subprocess helper ───────────────────────────────────────────────────
# Running scanpy inside the same R process causes a segfault on macOS ARM64:
# numpy's LAPACK (scipy) conflicts with R's Accelerate BLAS.
# Solution: write data to temp CSV, run Python in a fully isolated subprocess
# via `conda run`, and read the results back.

.conda_bin <- "/opt/homebrew/Caskroom/miniconda/base/bin/conda"

.run_py_subprocess <- function(py_lines, condaenv = .python_env) {
  script <- tempfile(fileext = ".py")
  writeLines(py_lines, script)
  on.exit(unlink(script), add = TRUE)
  system2(.conda_bin,
          args   = c("run", "--no-capture-output", "-n", condaenv, "python", script),
          stdout = TRUE, stderr = TRUE, wait = TRUE)
}

# ── 7. PAGA-DPT (scanpy) ──────────────────────────────────────────────────────
#' PAGA-guided diffusion pseudotime (scanpy) — runs in a subprocess
#'
#' Runs scanpy's PAGA graph abstraction to improve DPT connectivity, then
#' computes Diffusion Pseudotime anchored at start_cell.
#'
#' Python requirements: scanpy, anndata  (conda env: sc_tutorial)
#'
#' @param expr       genes × cells log-normalised matrix
#' @param start_cell name of the root cell
#' @param precomp_pca optional cells × PCs matrix; passed directly to AnnData
#' @return Named numeric [0, 1] per cell
run_paga <- function(expr, start_cell = NULL, precomp_pca = NULL) {
  pca_mat <- if (!is.null(precomp_pca))
    precomp_pca[colnames(expr), seq_len(min(30L, ncol(precomp_pca))), drop = FALSE]
  else
    .pca(expr, 30L)

  pca_file <- tempfile(fileext = ".csv")
  out_file  <- tempfile(fileext = ".csv")
  on.exit({ unlink(pca_file); unlink(out_file) }, add = TRUE)

  write.csv(as.data.frame(pca_mat), pca_file, row.names = TRUE)

  root_idx <- if (!is.null(start_cell) && start_cell %in% colnames(expr))
    which(colnames(expr) == start_cell) - 1L else 0L

  py <- c(
    "import warnings; warnings.filterwarnings('ignore')",
    "import pandas as pd, numpy as np, scanpy as sc, anndata as ad",
    sprintf("pca_df   = pd.read_csv('%s', index_col=0)", pca_file),
    "cells    = pca_df.index.tolist()",
    "pca      = pca_df.values.astype('float32')",
    sprintf("root_idx = %d", root_idx),
    sprintf("out_file = '%s'", out_file),
    "adata = ad.AnnData(np.zeros((len(cells), 1)))",
    "adata.obs_names = cells",
    "adata.obsm['X_pca'] = pca",
    "sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')",
    # PAGA requires cluster labels (leiden) before it can run
    "sc.tl.leiden(adata, resolution=0.5)",
    "sc.tl.paga(adata, groups='leiden')",
    "adata.uns['iroot'] = int(root_idx)",
    "sc.tl.dpt(adata)",
    "adata.obs[['dpt_pseudotime']].to_csv(out_file)"
  )

  message("  [PAGA-DPT] running PAGA+DPT in Python subprocess ...")
  out <- .run_py_subprocess(py)
  if (!file.exists(out_file) || file.info(out_file)$size == 0) {
    message("  [PAGA-DPT] subprocess failed:\n", paste(tail(out, 20), collapse = "\n"))
    return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))
  }

  pt_df <- read.csv(out_file, row.names = 1, check.names = FALSE)
  pt    <- pt_df[[1]][match(colnames(expr), rownames(pt_df))]
  .rescale01(setNames(as.numeric(pt), colnames(expr)))
}


# ── 8. Palantir ───────────────────────────────────────────────────────────────
#' Palantir pseudotime via multiscale diffusion components
#'
#' Runs diffusion maps on the PCA embedding, projects into a multiscale
#' diffusion space, then computes Palantir pseudotime from start_cell.
#' Palantir requires an explicit root cell (there is no unrooted mode).
#'
#' Python requirements: palantir, scanpy, anndata  (conda env: see PYTHON SETUP)
#'
#' @param expr        genes × cells log-normalised matrix
#' @param start_cell  name of the root cell; falls back to the first cell if NULL
#' @param precomp_pca optional cells × PCs matrix; passed directly to AnnData
#' @return Named numeric [0, 1] per cell
run_palantir <- function(expr, start_cell = NULL, precomp_pca = NULL) {
  pca_mat <- if (!is.null(precomp_pca))
    precomp_pca[colnames(expr), seq_len(min(30L, ncol(precomp_pca))), drop = FALSE]
  else
    .pca(expr, 30L)

  pca_file <- tempfile(fileext = ".csv")
  out_file  <- tempfile(fileext = ".csv")
  on.exit({ unlink(pca_file); unlink(out_file) }, add = TRUE)

  write.csv(as.data.frame(pca_mat), pca_file, row.names = TRUE)

  root_cell <- if (!is.null(start_cell) && start_cell %in% colnames(expr))
    start_cell else colnames(expr)[1]

  py <- c(
    "import warnings; warnings.filterwarnings('ignore')",
    "import pandas as pd, numpy as np, scanpy as sc, anndata as ad, palantir",
    sprintf("pca_df    = pd.read_csv('%s', index_col=0)", pca_file),
    "cells     = pca_df.index.tolist()",
    "pca       = pca_df.values.astype('float32')",
    sprintf("root_cell = '%s'", root_cell),
    sprintf("out_file  = '%s'", out_file),
    "adata = ad.AnnData(np.zeros((len(cells), 1)))",
    "adata.obs_names = cells",
    "adata.obsm['X_pca'] = pca",
    "sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')",
    "palantir.utils.run_diffusion_maps(adata, pca_key='X_pca')",
    "palantir.utils.determine_multiscale_space(adata)",
    "palantir.core.run_palantir(adata, early_cell=root_cell)",
    "adata.obs[['palantir_pseudotime']].to_csv(out_file)"
  )

  message("  [Palantir] running diffusion maps + Palantir in Python subprocess ...")
  out <- .run_py_subprocess(py)
  if (!file.exists(out_file) || file.info(out_file)$size == 0) {
    message("  [Palantir] subprocess failed:\n", paste(tail(out, 20), collapse = "\n"))
    return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))
  }

  pt_df <- read.csv(out_file, row.names = 1, check.names = FALSE)
  pt    <- pt_df[[1]][match(colnames(expr), rownames(pt_df))]
  .rescale01(setNames(as.numeric(pt), colnames(expr)))
}

# =============================================================================
# ORCHESTRATOR
# =============================================================================

#' Run all eight TI methods and return a tidy data frame of pseudotimes
#'
#' Each method is wrapped in tryCatch so a single failure does not abort the
#' run.  Methods are skipped gracefully when their dependencies are absent.
#'
#' @param expr        genes × cells log-normalised matrix
#' @param start_cell  name of the root cell (NULL = auto-detect)
#' @param methods     character vector of method names to run; default = all
#' @param n_pcs       number of PCs for embedding-based methods
#' @return data.frame with one row per cell and one column per method,
#'         values in [0, 1] (NA where unavailable)
run_all_ti_methods <- function(
    expr,
    start_cell         = NULL,
    fullgene_expr      = NULL,   # all-gene log-norm matrix for CytoTRACE; falls back to expr
    precomp_pca        = NULL,   # optional cells × PCs matrix passed to Slingshot / DPT
    methods            = c("Slingshot", "CytoTRACE", "Monocle3", "DPT", "SCORPIUS",
                           "TSCAN", "PAGA-DPT", "Palantir"),
    n_pcs              = 30L
) {
  stopifnot(is.matrix(expr), !is.null(rownames(expr)), !is.null(colnames(expr)))

  cells <- colnames(expr)

  dispatch <- list(
    Slingshot   = function() run_slingshot(expr,   start_cell, n_pcs, precomp_pca),
    CytoTRACE   = function() run_cytotrace(expr, fullgene_expr),
    Monocle3    = function() run_monocle3(expr,    start_cell),
    DPT         = function() run_dpt(expr,         start_cell, n_pcs, precomp_pca),
    SCORPIUS    = function() run_scorpius(expr, start_cell),
    TSCAN       = function() run_tscan(expr,       start_cell),
    "PAGA-DPT"  = function() run_paga(expr,        start_cell, precomp_pca),
    Palantir    = function() run_palantir(expr,    start_cell, precomp_pca)
  )

  results <- vector("list", length(methods))
  names(results) <- methods

  for (m in methods) {
    message(sprintf("\n=== %s ===", m))
    results[[m]] <- tryCatch(
      dispatch[[m]](),
      error = function(e) {
        message(sprintf("  [%s] unexpected error: %s", m, e$message))
        setNames(rep(NA_real_, length(cells)), cells)
      }
    )
    # Align to cell order in expr (defensive)
    results[[m]] <- results[[m]][cells]
  }

  as.data.frame(results, row.names = cells, check.names = FALSE)
}

# =============================================================================
# ENTRY POINT
# =============================================================================

if (!interactive() && sys.nframe() == 0L) {
  message("Building linear fixture (n_cells = 300) ...")
  f <- make_linear_fixture(n_cells = 300)

  res <- run_all_ti_methods(
    expr       = f$expr,
    start_cell = names(f$pseudotime)[1]
  )

  message("\n── Spearman correlation with ground-truth pseudotime ──")
  cors <- sapply(res, function(pt) {
    cor(pt, f$pseudotime, use = "pairwise.complete.obs", method = "spearman")
  })
  print(round(cors, 3))

  out <- file.path(.repo_root, "manuscript", "ti_pseudotime_results.csv")
  write.csv(cbind(ground_truth = f$pseudotime, res), out, row.names = TRUE)
  message("\nResults written to: ", out)
}
