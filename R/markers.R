# =============================================================================
# markers.R — Lightweight marker retrieval adapters
#
# All get_markers_*() functions return a marker_set: a small S3 object pairing
# early and terminal gene vectors with source metadata.  The object's $early
# and $terminal slots are plain character vectors, so they drop directly into
# the early_markers / terminal_markers arguments of any BioTrajX function.
#
# External dependencies are optional (Suggests):
#   msigdbr         — get_markers_msigdb()
#   AnnotationDbi   — get_markers_go()
#   org.Hs.eg.db    — get_markers_go() default organism
# =============================================================================


# ── Internal constructor ──────────────────────────────────────────────────────

.marker_set <- function(early, terminal, source = "custom", metadata = list()) {
  structure(
    list(
      early    = unique(as.character(early)),
      terminal = unique(as.character(terminal)),
      source   = source,
      metadata = metadata
    ),
    class = "marker_set"
  )
}


# ── print.marker_set ───────────────────────────────────────────────────────────

#' Print a marker_set
#'
#' @param x A `marker_set` object.
#' @param n Number of example genes to display per group. Default 6.
#' @param ... Unused.
#' @export
print.marker_set <- function(x, n = 6L, ...) {
  preview <- function(v) {
    if (length(v) == 0L) return("(none)")
    ex <- paste(head(v, n), collapse = ", ")
    if (length(v) > n) ex <- paste0(ex, ", ...")
    sprintf("%s  [%d genes]", ex, length(v))
  }
  cat("marker_set  <source:", x$source, ">\n")
  cat("  early   :", preview(x$early),    "\n")
  cat("  terminal:", preview(x$terminal), "\n")
  invisible(x)
}


# ── get_markers_msigdb ────────────────────────────────────────────────────────

#' Retrieve marker genes from MSigDB via msigdbr
#'
#' Queries the Molecular Signatures Database and returns a [marker_set] whose
#' `$early` and `$terminal` slots can be passed directly to any BioTrajX
#' function.
#'
#' @param early Character vector of gene set name(s) representing the early /
#'   root state.  When `exact = FALSE` each element is treated as a
#'   [grepl()] pattern matched against `gs_name`.
#' @param terminal Character vector of gene set name(s) representing the
#'   terminal / differentiated state.
#' @param collection MSigDB collection identifier, or `NULL` (default) to
#'   search all collections.  Commonly used values: `"C2"` (curated sets),
#'   `"C5"` (GO-based), `"C7"` (immunologic signatures), `"C8"`
#'   (cell-type signatures).  Passed as the `collection` argument to
#'   [msigdbr::msigdbr()] (msigdbr >= 10.0.0); for older versions the
#'   value is forwarded via the deprecated `category` argument automatically.
#' @param subcollection Optional MSigDB sub-collection (e.g. `"GO:BP"`).
#' @param species Species name passed to [msigdbr::msigdbr()].
#'   Default `"Homo sapiens"`.
#' @param exact Logical.  If `TRUE` (default) set names must match exactly;
#'   if `FALSE` each element of `early` / `terminal` is used as a pattern
#'   via [grepl()], and matching sets are unioned.
#'
#' @return A `marker_set` with `$early` and `$terminal` character vectors.
#'
#' @examples
#' \dontrun{
#' ms <- get_markers_msigdb(
#'   early      = "GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_UP",
#'   terminal   = "GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_UP",
#'   collection = "C7"
#' )
#' compute_single_doe_linear(expr, pt, ms$early, ms$terminal, ...)
#' }
#' @export
get_markers_msigdb <- function(early,
                                terminal,
                                collection    = NULL,
                                subcollection = NULL,
                                species       = "Homo sapiens",
                                exact         = TRUE) {
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("Package 'msigdbr' is required. Install with: install.packages('msigdbr')")

  # msigdbr >= 10.0.0 uses `collection`; earlier versions used `category`.
  # Detect which API is available and call accordingly.
  gs <- tryCatch(
    suppressWarnings(
      msigdbr::msigdbr(species = species, collection = collection,
                       subcollection = subcollection)
    ),
    error = function(e) {
      # Fall back to the old category/subcategory API
      msigdbr::msigdbr(species = species, category = collection,
                       subcategory = subcollection)
    }
  )

  .pull <- function(names) {
    if (exact) {
      rows <- gs[gs$gs_name %in% names, ]
    } else {
      pat  <- paste(names, collapse = "|")
      rows <- gs[grepl(pat, gs$gs_name, ignore.case = TRUE), ]
    }
    if (nrow(rows) == 0L)
      stop("No gene sets found for: ", paste(names, collapse = ", "),
           "\n  collection=", if (is.null(collection)) "NULL (all)" else paste0("'", collection, "'"),
           "  species='", species, "'")
    unique(rows$gene_symbol)
  }

  .marker_set(
    early    = .pull(early),
    terminal = .pull(terminal),
    source   = "MSigDB",
    metadata = list(
      collection    = collection,
      subcollection = subcollection,
      species       = species,
      early_sets    = early,
      terminal_sets = terminal
    )
  )
}


# ── get_markers_go ────────────────────────────────────────────────────────────

#' Retrieve marker genes from GO terms via AnnotationDbi
#'
#' Maps GO term IDs to gene symbols using an OrgDb annotation package and
#' returns a [marker_set].  Requires `AnnotationDbi` and an OrgDb package
#' (defaults to `org.Hs.eg.db` for human).
#'
#' @param early_terms Character vector of GO IDs (e.g. `"GO:0045588"`) for
#'   the early / root state.
#' @param terminal_terms Character vector of GO IDs for the terminal state.
#' @param org_db An OrgDb object (e.g. `org.Hs.eg.db::org.Hs.eg.db`).
#'   If `NULL`, defaults to human via `org.Hs.eg.db`.
#' @param ont Ontology filter: `"BP"`, `"CC"`, `"MF"`, or `"ALL"` (default).
#'   Applied to the `ONTOLOGYALL` column returned by [AnnotationDbi::select()].
#'
#' @return A `marker_set` with `$early` and `$terminal` character vectors.
#'
#' @examples
#' \dontrun{
#' library(org.Hs.eg.db)
#' ms <- get_markers_go(
#'   early_terms    = "GO:0045624",  # positive regulation of T-helper cell diff.
#'   terminal_terms = "GO:0002250",  # adaptive immune response
#'   org_db         = org.Hs.eg.db
#' )
#' compute_single_doe_linear(expr, pt, ms$early, ms$terminal, ...)
#' }
#' @export
get_markers_go <- function(early_terms,
                            terminal_terms,
                            org_db = NULL,
                            ont    = "ALL") {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE))
    stop("Package 'AnnotationDbi' is required. Install with:\n",
         "  BiocManager::install('AnnotationDbi')")

  if (is.null(org_db)) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop("Package 'org.Hs.eg.db' is required. Install with:\n",
           "  BiocManager::install('org.Hs.eg.db')")
    org_db <- org.Hs.eg.db::org.Hs.eg.db
  }

  .pull <- function(terms) {
    res <- tryCatch(
      suppressMessages(
        AnnotationDbi::select(org_db,
                              keys    = terms,
                              keytype = "GOALL",
                              columns = "SYMBOL")
      ),
      error = function(e) stop("AnnotationDbi query failed: ", e$message)
    )
    if (ont != "ALL" && "ONTOLOGYALL" %in% names(res))
      res <- res[res$ONTOLOGYALL == ont, ]
    genes <- unique(res$SYMBOL)
    genes <- genes[!is.na(genes)]
    if (length(genes) == 0L)
      stop("No genes found for GO terms: ", paste(terms, collapse = ", "))
    genes
  }

  .marker_set(
    early    = .pull(early_terms),
    terminal = .pull(terminal_terms),
    source   = "GO",
    metadata = list(
      org_db         = class(org_db)[1L],
      ont            = ont,
      early_terms    = early_terms,
      terminal_terms = terminal_terms
    )
  )
}


# ── get_markers_cellmarker ────────────────────────────────────────────────────

# Internal: download and normalise a CellMarker table
.fetch_cellmarker <- function(species) {
  urls <- c(
    Human = "http://xteam.xbio.top/CellMarker/download/Human_cell_markers.txt",
    Mouse = "http://xteam.xbio.top/CellMarker/download/Mouse_cell_markers.txt"
  )
  url <- urls[[species]]
  message("Downloading CellMarker (", species, ") from:\n  ", url)
  tmp <- tempfile(fileext = ".txt")
  tryCatch(
    utils::download.file(url, tmp, quiet = TRUE, mode = "wb"),
    error = function(e) stop(
      "Download failed: ", e$message, "\n",
      "Load the file manually instead:\n",
      "  df <- read.table('Cell_marker_", species, ".txt',\n",
      "                   sep = '\\t', header = TRUE, quote = '', fill = TRUE)\n",
      "  get_markers_cellmarker(..., df = df)"
    )
  )
  df <- read.table(tmp, sep = "\t", header = TRUE, quote = "",
                   fill = TRUE, stringsAsFactors = FALSE)

  # Normalise column names: CellMarker v2 uses snake_case, v1 uses camelCase.
  # Remap v2 names to v1 so the rest of the function works with one set of defaults.
  remap <- c(cell_name = "cellName", gene_symbol = "geneSymbol",
             tissue_type = "tissueType")
  for (old in names(remap)) {
    new <- remap[[old]]
    if (old %in% names(df) && !new %in% names(df))
      names(df)[names(df) == old] <- new
  }

  df
}

# Internal: expand rows where a column contains comma-separated values
.expand_comma_col <- function(df, col) {
  if (!col %in% names(df) || !any(grepl(",", df[[col]], fixed = TRUE)))
    return(df)
  df <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    vals  <- trimws(strsplit(df[[col]][i], ",", fixed = TRUE)[[1L]])
    vals  <- vals[nchar(vals) > 0L]
    row_i <- df[rep(i, length(vals)), , drop = FALSE]
    row_i[[col]] <- vals
    row_i
  }))
  rownames(df) <- NULL
  df
}

#' Retrieve marker genes from CellMarker
#'
#' Fetches cell-type marker genes from the CellMarker database and returns a
#' [marker_set].  When `df = NULL` (default) the relevant species table is
#' downloaded automatically from
#' <http://xteam.xbio.top/CellMarker/download.jsp>.  Pass a pre-loaded data
#' frame to `df` to work offline or with a custom table.
#'
#' @param early Character vector of cell type name(s) for the early / root
#'   state, matched against `cell_type_col`.
#' @param terminal Character vector of cell type name(s) for the terminal
#'   state.
#' @param df Optional data frame loaded from a CellMarker download.  When
#'   `NULL` (default) the file is downloaded using `species`.
#' @param species `"Human"` (default) or `"Mouse"`.  Ignored when `df` is
#'   provided.
#' @param cell_type_col Column name containing cell type labels.
#'   Default `"cellName"`.
#' @param marker_col Column name containing gene symbols.
#'   Default `"geneSymbol"`.
#' @param tissue Optional character vector of tissue types to pre-filter on.
#' @param tissue_col Column name for tissue types.  Default `"tissueType"`.
#'
#' @return A `marker_set` with `$early` and `$terminal` character vectors.
#'
#' @examples
#' \dontrun{
#' # Auto-download (requires internet)
#' ms <- get_markers_cellmarker(
#'   early    = "Naive T cell",
#'   terminal = "Exhausted T cell",
#'   species  = "Human",
#'   tissue   = "Blood"
#' )
#'
#' # Or supply a pre-downloaded file
#' df <- read.table("Human_cell_markers.txt",
#'                  sep = "\t", header = TRUE, quote = "", fill = TRUE)
#' ms <- get_markers_cellmarker("Naive T cell", "Exhausted T cell", df = df)
#'
#' compute_single_doe_linear(expr, pt, ms$early, ms$terminal, ...)
#' }
#' @export
get_markers_cellmarker <- function(early,
                                    terminal,
                                    df            = NULL,
                                    species       = c("Human", "Mouse"),
                                    cell_type_col = "cellName",
                                    marker_col    = "geneSymbol",
                                    tissue        = NULL,
                                    tissue_col    = "tissueType") {
  if (is.null(df)) {
    species <- match.arg(species)
    df <- .fetch_cellmarker(species)
  }

  if (!is.data.frame(df))
    stop("'df' must be a data.frame.")

  # Expand comma-separated gene symbols (present in some CellMarker versions)
  df <- .expand_comma_col(df, marker_col)

  for (col in c(cell_type_col, marker_col)) {
    if (!col %in% names(df))
      stop("Column '", col, "' not found in df. ",
           "Available: ", paste(names(df), collapse = ", "))
  }

  if (!is.null(tissue)) {
    if (!tissue_col %in% names(df))
      stop("Column '", tissue_col, "' not found in df.")
    df <- df[df[[tissue_col]] %in% tissue, ]
  }

  .pull <- function(types) {
    rows <- df[df[[cell_type_col]] %in% types, ]
    if (nrow(rows) == 0L)
      stop("No entries found for: ", paste(types, collapse = ", "),
           "\n  Check cell type names and tissue filter.\n",
           "  Available types: ",
           paste(head(unique(df[[cell_type_col]]), 10L), collapse = ", "), " ...")
    genes <- unique(as.character(rows[[marker_col]]))
    genes[!is.na(genes) & nchar(genes) > 0L]
  }

  .marker_set(
    early    = .pull(early),
    terminal = .pull(terminal),
    source   = "CellMarker",
    metadata = list(
      early_types    = early,
      terminal_types = terminal,
      species        = if (is.null(df)) species else NULL,
      tissue         = tissue
    )
  )
}


# ── markers_to_list ───────────────────────────────────────────────────────────

#' Convert a marker_set or data frame to a named list
#'
#' Produces the named-list format required by [compute_single_doe_branched()]
#' and [compute_multi_doe_branched()].
#'
#' - For a `marker_set`: returns `list(early = ..., terminal = ...)`.
#' - For a `data.frame`: splits gene symbols by cell-type label, one element
#'   per unique cell type.  Useful when a single table covers multiple branches.
#'
#' @param x A `marker_set` or a `data.frame`.
#' @param cell_type_col `data.frame` method only — column of cell type labels.
#'   Default `"cell_type"`.
#' @param marker_col `data.frame` method only — column of gene symbols.
#'   Default `"gene_symbol"`.
#' @param ... Unused.
#'
#' @return A named list of character vectors.
#'
#' @examples
#' ms <- get_markers_cellmarker("Naive T cell", "Exhausted T cell", df = cm)
#' lst <- markers_to_list(ms)
#' # lst$early and lst$terminal ready for compute_single_doe_branched()
#'
#' @export
markers_to_list <- function(x, ...) UseMethod("markers_to_list")

#' @rdname markers_to_list
#' @export
markers_to_list.marker_set <- function(x, ...) {
  list(early = x$early, terminal = x$terminal)
}

#' @rdname markers_to_list
#' @export
markers_to_list.data.frame <- function(x,
                                        cell_type_col = "cell_type",
                                        marker_col    = "gene_symbol",
                                        ...) {
  for (col in c(cell_type_col, marker_col)) {
    if (!col %in% names(x))
      stop("Column '", col, "' not found. Available: ",
           paste(names(x), collapse = ", "))
  }
  lapply(
    split(x[[marker_col]], x[[cell_type_col]]),
    function(g) unique(as.character(g[!is.na(g) & nchar(as.character(g)) > 0L]))
  )
}


# ── filter_markers ────────────────────────────────────────────────────────────

#' Filter a marker_set to high-quality genes using expression data
#'
#' Applies a two-step expression-based filter to a [marker_set]:
#'
#' 1. **Dropout filter** — removes genes detected (non-zero) in fewer than
#'    `min_detection` fraction of cells.  Genes that are zero in most cells
#'    contribute noise rather than signal to module scores.
#' 2. **Top-N selection** — among genes that pass the dropout filter, ranks by
#'    mean non-zero expression and keeps the top `top_n`.  Higher mean
#'    expression means a stronger, more detectable signal per cell.
#'
#' Both steps are purely expression-based with no pseudotime or cluster
#' assumptions, so there is no circularity with downstream DOE scoring.
#' Accepts dense matrices, sparse [Matrix::dgCMatrix] objects, or [Seurat]
#' objects (normalised data layer is extracted automatically).
#'
#' @param ms A [marker_set] returned by `get_markers_msigdb()`,
#'   `get_markers_cellmarker()`, `get_markers_go()`, or `.marker_set()`.
#' @param expr A genes-by-cells expression matrix (dense or sparse), or a
#'   [Seurat] object.  Should contain log-normalised counts.
#' @param top_n Integer.  Maximum number of genes to retain per group after
#'   dropout filtering.  Genes are ranked by mean non-zero expression and the
#'   top `top_n` are kept.  Pass `NULL` to skip this step and only apply the
#'   dropout filter.  Default `30L`.
#' @param min_detection Numeric in \[0, 1\].  Minimum fraction of cells in
#'   which a gene must be detected (non-zero) to be retained.  Default `0.10`.
#'
#' @return A [marker_set] with filtered `$early` and `$terminal` vectors.
#'   The `$source` field is appended with `" [filtered]"` and the original
#'   gene counts and filter parameters are stored in `$metadata`.
#'
#' @examples
#' \dontrun{
#' ms <- get_markers_cellmarker(
#'   early    = "Hematopoietic stem cell",
#'   terminal = c("Erythroid cell", "Red blood cell (erythrocyte)"),
#'   species  = "Mouse"
#' )
#' ms <- filter_markers(ms, seurat_obj, top_n = 30, min_detection = 0.10)
#' compute_single_doe_linear(seurat_obj, pseudotime, ms$early, ms$terminal)
#' }
#' @export
filter_markers <- function(ms, expr, top_n = 30L, min_detection = 0.10) {

  if (!inherits(ms, "marker_set"))
    stop("'ms' must be a marker_set object returned by get_markers_*().")

  # Accept Seurat objects — extract the normalised data layer
  if (inherits(expr, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE))
      stop("Package 'Seurat' is required when 'expr' is a Seurat object.")
    expr <- Seurat::GetAssayData(expr, layer = "data")
  }

  if (!is.matrix(expr) && !inherits(expr, "Matrix"))
    stop("'expr' must be a matrix, sparse Matrix (dgCMatrix), or Seurat object.")

  if (!is.null(top_n)) {
    top_n <- as.integer(top_n)
    if (top_n < 1L) stop("'top_n' must be a positive integer or NULL.")
  }

  if (min_detection < 0 || min_detection > 1)
    stop("'min_detection' must be between 0 and 1.")

  # Row-wise helpers — work for both dense and sparse matrices
  .detection_rate <- function(mat) {
    if (inherits(mat, "Matrix")) Matrix::rowMeans(mat != 0)
    else rowMeans(mat != 0)
  }

  .mean_nonzero <- function(mat) {
    if (inherits(mat, "Matrix")) {
      rs  <- Matrix::rowSums(mat)
      nnz <- Matrix::rowSums(mat != 0)
      ifelse(nnz > 0L, rs / nnz, 0)
    } else {
      apply(mat, 1L, function(x) { nz <- x[x > 0L]; if (length(nz)) mean(nz) else 0 })
    }
  }

  .filter_one <- function(genes, label) {
    genes_in  <- intersect(genes, rownames(expr))
    n_missing <- length(genes) - length(genes_in)
    if (n_missing > 0L)
      message(sprintf("  [%s] %d gene(s) not in expr — dropped", label, n_missing))
    if (length(genes_in) == 0L) {
      warning(sprintf("[%s] No genes overlap with the expression matrix.", label))
      return(character(0L))
    }

    # Step 1: dropout filter
    det     <- .detection_rate(expr[genes_in, , drop = FALSE])
    pass    <- genes_in[det >= min_detection]
    n_drop  <- length(genes_in) - length(pass)
    if (n_drop > 0L)
      message(sprintf("  [%s] %d gene(s) below %.0f%% detection — dropped",
                      label, n_drop, min_detection * 100))

    if (length(pass) == 0L) {
      warning(sprintf(
        "[%s] No genes passed min_detection = %.2f. Returning all intersected genes.",
        label, min_detection))
      pass <- genes_in
    }

    # Step 2: rank by mean non-zero expression, keep top_n
    if (!is.null(top_n) && length(pass) > top_n) {
      mean_nz <- .mean_nonzero(expr[pass, , drop = FALSE])
      pass    <- pass[order(mean_nz, decreasing = TRUE)[seq_len(top_n)]]
      message(sprintf("  [%s] kept top %d by mean non-zero expression", label, top_n))
    }

    pass
  }

  message(sprintf(
    "filter_markers: early=%d, terminal=%d genes before filtering",
    length(ms$early), length(ms$terminal)))

  early_out <- .filter_one(ms$early,    "early")
  term_out  <- .filter_one(ms$terminal, "terminal")

  message(sprintf(
    "filter_markers: early=%d, terminal=%d genes after filtering",
    length(early_out), length(term_out)))

  .marker_set(
    early    = early_out,
    terminal = term_out,
    source   = paste0(ms$source, " [filtered]"),
    metadata = c(ms$metadata, list(
      filter_top_n         = top_n,
      filter_min_detection = min_detection,
      n_early_before       = length(ms$early),
      n_term_before        = length(ms$terminal)
    ))
  )
}
