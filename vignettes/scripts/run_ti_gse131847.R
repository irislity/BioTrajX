# =============================================================================
# run_ti_gse131847.R
#
# Run all TI methods on the GSE131847 LCMV time-course CD8 T-cell dataset
# (naive, d3, d4, d5, d6, d7, d10, d14, d21, d32, d60, d90 post-infection)
# and save pseudotimes.
#
# INPUT:   data/GSE131847_seu.rds
# OUTPUT:  manuscript/results/linear_gse131847/ti_pseudotime_gse131847.csv
#
# Usage:
#   Rscript manuscript/scripts/real/run_ti_gse131847.R
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
message("Loading GSE131847_seu.rds ...")
seurat <- readRDS(file.path(repo_root, "data", "GSE131847_seu.rds"))
message(sprintf("  %d cells x %d genes", ncol(seurat), nrow(seurat)))
message("  cell_type (day) distribution:")
print(table(seurat$cell_type))

DefaultAssay(seurat) <- "SCT"

# =============================================================================
# 2. PREPARE EXPRESSION MATRICES
# =============================================================================
seurat_pca <- Embeddings(seurat, "pca")

message("Selecting top 2000 highly variable genes ...")
seurat <- FindVariableFeatures(seurat, nfeatures = 2000, verbose = FALSE)
hvg    <- VariableFeatures(seurat)

expr          <- as.matrix(GetAssayData(seurat, assay = "SCT", layer = "data")[hvg, ])
fullgene_expr <- as.matrix(GetAssayData(seurat, assay = "SCT", layer = "data"))

message(sprintf("  log-norm HVG:  %d x %d", nrow(expr),          ncol(expr)))
message(sprintf("  log-norm (all):%d x %d", nrow(fullgene_expr), ncol(fullgene_expr)))

# =============================================================================
# 3. ROOT CELL — most primitive cell among cells labelled "naive"
# =============================================================================
# Ground-truth day-of-infection is known here, so anchor the root in the
# naive population rather than picking it up genome-wide via CytoTRACE alone.
message("Computing CytoTRACE score for root cell selection (restricted to naive cells) ...")
cytotrace_pt <- run_cytotrace(expr, fullgene_expr)
naive_cells  <- colnames(seurat)[seurat$cell_type == "naive"]
naive_cells  <- intersect(naive_cells, names(cytotrace_pt))
# run_cytotrace() inverts: 0 = primitive, 1 = differentiated -- so root = min
start_cell   <- naive_cells[which.min(cytotrace_pt[naive_cells])]
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
out_dir <- file.path(repo_root, "manuscript", "results", "linear_gse131847")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "ti_pseudotime_gse131847.csv")
write.csv(ti_results, out, row.names = TRUE)
message("\nSaved: ", out)
message("Done.")
