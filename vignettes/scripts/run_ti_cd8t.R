# =============================================================================
# run_ti_cd8t.R
#
# Run all TI methods on the CD8 T-cell dataset and save pseudotimes.
#
# INPUT:   data/cd8t.rds
# OUTPUT:  manuscript/results/linear_cd8t/ti_pseudotime_cd8t.csv
#
# Usage:
#   Rscript manuscript/scripts/real/run_ti_cd8t.R
# =============================================================================

Sys.setenv(OMP_NUM_THREADS        = "1",
           OPENBLAS_NUM_THREADS   = "1",
           MKL_NUM_THREADS        = "1",
           VECLIB_MAXIMUM_THREADS = "1")

suppressPackageStartupMessages(library(Seurat))

repo_root <- here::here()
source(file.path(repo_root, "manuscript", "scripts", "real", "run_ti_methods.R"))

# =============================================================================
# 1. LOAD DATA
# =============================================================================
message("Loading cd8t.rds ...")
seurat <- readRDS(file.path(repo_root, "data", "cd8t.rds"))
message(sprintf("  %d cells × %d genes", ncol(seurat), nrow(seurat)))

# =============================================================================
# 2. PREPARE EXPRESSION MATRICES
# =============================================================================
seurat_pca <- Embeddings(seurat, "pca")

message("Selecting top 2000 highly variable genes ...")
seurat <- FindVariableFeatures(seurat, nfeatures = 2000, verbose = FALSE)
hvg    <- VariableFeatures(seurat)

expr          <- as.matrix(GetAssayData(seurat, layer = "data")[hvg, ])
fullgene_expr <- as.matrix(GetAssayData(seurat, layer = "data"))

message(sprintf("  log-norm HVG:  %d × %d", nrow(expr),          ncol(expr)))
message(sprintf("  log-norm (all):%d × %d", nrow(fullgene_expr), ncol(fullgene_expr)))

# =============================================================================
# 3. ROOT CELL — most primitive cell by CytoTRACE score
# =============================================================================
message("Computing CytoTRACE score for root cell selection ...")
cytotrace_pt <- run_cytotrace(expr, fullgene_expr)
# run_cytotrace() inverts: 0 = primitive, 1 = differentiated — so root = min
start_cell   <- names(which.min(cytotrace_pt))
message(sprintf("  Root cell: %s", start_cell))

# =============================================================================
# 4. RUN TI METHODS
# =============================================================================
set.seed(42)
message("\nRunning TI methods ...")
ti_results <- run_all_ti_methods(
  expr          = expr,
  fullgene_expr = fullgene_expr,
  start_cell    = start_cell,
  precomp_pca   = seurat_pca,
  methods       = c("Slingshot", "CytoTRACE", "Monocle3", "DPT", "SCORPIUS",
                    "TSCAN", "PAGA-DPT", "Palantir")
)

# =============================================================================
# 5. SAVE
# =============================================================================
out_dir <- file.path(repo_root, "manuscript", "results", "linear_cd8t")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "ti_pseudotime_cd8t.csv")
write.csv(ti_results, out, row.names = TRUE)
message("\nSaved: ", out)
message("Done.")
