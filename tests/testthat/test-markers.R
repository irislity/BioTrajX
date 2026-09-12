# =============================================================================
# test-markers.R — Unit tests for marker_set and get_markers_*()
# =============================================================================
library(testthat)

# ── helpers ───────────────────────────────────────────────────────────────────

make_cellmarker_df <- function() {
  data.frame(
    cellName   = c("Naive T cell", "Naive T cell", "Exhausted T cell",
                   "Exhausted T cell", "B cell"),
    geneSymbol = c("TCF7", "LEF1", "GZMB", "PRF1", "CD19"),
    tissueType = c("Blood", "Blood", "Blood", "Tumor", "Blood"),
    stringsAsFactors = FALSE
  )
}

# ── print.marker_set ───────────────────────────────────────────────────────────

test_that("print.marker_set: outputs source and gene counts", {
  ms <- BioTrajX:::.marker_set(c("TCF7", "LEF1"), c("GZMB", "PRF1"),
                                source = "test")
  expect_output(print(ms), "marker_set")
  expect_output(print(ms), "early")
  expect_output(print(ms), "terminal")
})

test_that("print.marker_set: returns invisibly", {
  ms <- BioTrajX:::.marker_set(c("A"), c("B"), source = "test")
  expect_invisible(print(ms))
})

# ── get_markers_cellmarker ────────────────────────────────────────────────────

test_that("get_markers_cellmarker: returns marker_set", {
  cm <- make_cellmarker_df()
  ms <- get_markers_cellmarker("Naive T cell", "Exhausted T cell", df = cm)
  expect_s3_class(ms, "marker_set")
})

test_that("get_markers_cellmarker: correct early and terminal genes", {
  cm <- make_cellmarker_df()
  ms <- get_markers_cellmarker("Naive T cell", "Exhausted T cell", df = cm)
  expect_setequal(ms$early,    c("TCF7", "LEF1"))
  expect_setequal(ms$terminal, c("GZMB", "PRF1"))
})

test_that("get_markers_cellmarker: tissue filter works", {
  cm <- make_cellmarker_df()
  ms <- get_markers_cellmarker("Naive T cell", "Exhausted T cell",
                                df = cm, tissue = "Blood")
  # PRF1 is only in Tumor tissue — should be excluded
  expect_false("PRF1" %in% ms$terminal)
  expect_true("GZMB" %in% ms$terminal)
})

test_that("get_markers_cellmarker: source field is 'CellMarker'", {
  cm <- make_cellmarker_df()
  ms <- get_markers_cellmarker("Naive T cell", "Exhausted T cell", df = cm)
  expect_equal(ms$source, "CellMarker")
})

test_that("get_markers_cellmarker: non-data.frame df errors", {
  expect_error(
    get_markers_cellmarker("Naive T cell", "Exhausted T cell", df = list()),
    "data.frame"
  )
})

test_that("get_markers_cellmarker: comma-separated geneSymbol is expanded", {
  cm <- data.frame(
    cellName   = c("Naive T cell", "Exhausted T cell"),
    geneSymbol = c("TCF7, LEF1, CCR7", "GZMB, PRF1"),
    tissueType = c("Blood", "Blood"),
    stringsAsFactors = FALSE
  )
  ms <- get_markers_cellmarker("Naive T cell", "Exhausted T cell", df = cm)
  expect_true("TCF7" %in% ms$early)
  expect_true("LEF1" %in% ms$early)
  expect_true("GZMB" %in% ms$terminal)
})

test_that("get_markers_cellmarker: missing cell_type_col errors", {
  cm <- make_cellmarker_df()
  expect_error(
    get_markers_cellmarker("Naive T cell", "Exhausted T cell",
                            df = cm, cell_type_col = "no_such_col"),
    "no_such_col"
  )
})

test_that("get_markers_cellmarker: unknown cell type errors", {
  cm <- make_cellmarker_df()
  expect_error(
    get_markers_cellmarker("Naive T cell", "PlaceholderCell", df = cm),
    "PlaceholderCell"
  )
})

test_that("get_markers_cellmarker: custom column names work", {
  cm <- data.frame(
    type = c("stem", "stem", "effector"),
    gene = c("MYC", "SOX2", "GZMB"),
    stringsAsFactors = FALSE
  )
  ms <- get_markers_cellmarker("stem", "effector", df = cm,
                                cell_type_col = "type", marker_col = "gene")
  expect_setequal(ms$early,    c("MYC", "SOX2"))
  expect_setequal(ms$terminal, "GZMB")
})

# ── markers_to_list ───────────────────────────────────────────────────────────

test_that("markers_to_list.marker_set: returns named list with early and terminal", {
  ms  <- BioTrajX:::.marker_set(c("TCF7", "LEF1"), c("GZMB"), source = "test")
  lst <- markers_to_list(ms)
  expect_type(lst, "list")
  expect_named(lst, c("early", "terminal"), ignore.order = TRUE)
  expect_setequal(lst$early,    c("TCF7", "LEF1"))
  expect_setequal(lst$terminal, "GZMB")
})

test_that("markers_to_list.data.frame: splits by cell type", {
  df <- data.frame(
    cell_type   = c("A", "A", "B", "B"),
    gene_symbol = c("g1", "g2", "g3", "g4"),
    stringsAsFactors = FALSE
  )
  lst <- markers_to_list(df)
  expect_named(lst, c("A", "B"), ignore.order = TRUE)
  expect_setequal(lst$A, c("g1", "g2"))
  expect_setequal(lst$B, c("g3", "g4"))
})

test_that("markers_to_list.data.frame: custom column names work", {
  df <- data.frame(ct = c("X", "Y"), gn = c("ACTB", "GAPDH"),
                   stringsAsFactors = FALSE)
  lst <- markers_to_list(df, cell_type_col = "ct", marker_col = "gn")
  expect_named(lst, c("X", "Y"), ignore.order = TRUE)
})

test_that("markers_to_list.data.frame: drops NA genes", {
  df <- data.frame(
    cell_type   = c("A", "A", "A"),
    gene_symbol = c("g1", NA, ""),
    stringsAsFactors = FALSE
  )
  lst <- markers_to_list(df)
  expect_equal(lst$A, "g1")
})

test_that("markers_to_list.data.frame: missing column errors", {
  df <- data.frame(x = 1:3, y = letters[1:3], stringsAsFactors = FALSE)
  expect_error(markers_to_list(df), "cell_type")
})

test_that("markers_to_list: deduplicates genes", {
  ms <- BioTrajX:::.marker_set(c("A", "A", "B"), c("C"), source = "test")
  expect_equal(length(ms$early), 2L)  # constructor deduplicates
})

# ── get_markers_cellmarker auto-download ──────────────────────────────────────

test_that("get_markers_cellmarker: auto-downloads Human table when df = NULL", {
  skip_if_offline()
  ms <- tryCatch(
    get_markers_cellmarker("Naive T cell", "Exhausted T cell",
                            species = "Human", tissue = "Blood"),
    error = function(e) NULL
  )
  # If the download succeeded, we should get a valid marker_set
  skip_if(is.null(ms), "CellMarker server unreachable — skipping")
  expect_s3_class(ms, "marker_set")
  expect_gt(length(ms$early),    0L)
  expect_gt(length(ms$terminal), 0L)
})

# ── get_markers_msigdb ────────────────────────────────────────────────────────

test_that("get_markers_msigdb: errors without msigdbr", {
  skip_if(requireNamespace("msigdbr", quietly = TRUE),
          "msigdbr is installed — skipping absence test")
  expect_error(
    get_markers_msigdb("SETNAME_A", "SETNAME_B"),
    "msigdbr"
  )
})

test_that("get_markers_msigdb: returns marker_set with correct source", {
  skip_if_not_installed("msigdbr")
  # Use stable C7 immunologic signature sets (present in msigdbr >= 10.0.0)
  ms <- get_markers_msigdb(
    early      = "GOLDRATH_NAIVE_VS_EFF_CD8_TCELL_UP",
    terminal   = "GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_UP",
    collection = "C7"
  )
  expect_s3_class(ms, "marker_set")
  expect_equal(ms$source, "MSigDB")
  expect_gt(length(ms$early),    0L)
  expect_gt(length(ms$terminal), 0L)
  expect_type(ms$early,    "character")
  expect_type(ms$terminal, "character")
})

test_that("get_markers_msigdb: unknown set name errors", {
  skip_if_not_installed("msigdbr")
  expect_error(
    get_markers_msigdb("TOTALLY_NONEXISTENT_SET_XYZ",
                        "ALSO_NONEXISTENT_SET_ABC",
                        collection = "C7"),
    "No gene sets found"
  )
})

test_that("get_markers_msigdb: exact=FALSE pattern match works", {
  skip_if_not_installed("msigdbr")
  ms <- get_markers_msigdb("GOLDRATH_NAIVE", "GOLDRATH_NAIVE_VS_EFF",
                            collection = "C7", exact = FALSE)
  expect_gt(length(ms$early), 0L)
})

# ── get_markers_go ────────────────────────────────────────────────────────────

test_that("get_markers_go: errors without AnnotationDbi", {
  skip_if(requireNamespace("AnnotationDbi", quietly = TRUE),
          "AnnotationDbi is installed — skipping absence test")
  expect_error(get_markers_go("GO:0002376", "GO:0002250"), "AnnotationDbi")
})

test_that("get_markers_go: returns marker_set with correct source", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  ms <- get_markers_go(
    early_terms    = "GO:0045624",  # pos. reg. of T-helper diff.
    terminal_terms = "GO:0002250"   # adaptive immune response
  )
  expect_s3_class(ms, "marker_set")
  expect_equal(ms$source, "GO")
  expect_gt(length(ms$early),    0L)
  expect_gt(length(ms$terminal), 0L)
})

test_that("get_markers_go: ont filter reduces gene set", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  ms_all <- get_markers_go("GO:0002376", "GO:0002376", ont = "ALL")
  ms_bp  <- get_markers_go("GO:0002376", "GO:0002376", ont = "BP")
  # BP-only should be a subset of ALL
  expect_true(all(ms_bp$early %in% ms_all$early))
})

test_that("get_markers_go: invalid GO term errors", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  expect_error(
    get_markers_go("GO:9999999", "GO:9999999"),
    regexp = "failed|No genes found|select"
  )
})
