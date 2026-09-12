#!/usr/bin/env Rscript

library(grid)
library(hexSticker)
library(showtext)

.libPaths(c(
  .libPaths(),
  "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library"
))

# ── fonts ─────────────────────────────────────────────────────────────────────
font_add_google("Montserrat", "mont")
showtext_auto()

set.seed(21)

# ── helpers ───────────────────────────────────────────────────────────────────

ring_polygon <- function(
  cx, cy,
  r_out, r_in,
  from_deg, to_deg,
  n = 160
) {
  a_fwd <- seq(
    from_deg * pi / 180,
    to_deg * pi / 180,
    length.out = n
  )

  a_rev <- rev(a_fwd)

  list(
    x = c(
      cx + r_out * cos(a_fwd),
      cx + r_in  * cos(a_rev)
    ),
    y = c(
      cy + r_out * sin(a_fwd),
      cy + r_in  * sin(a_rev)
    )
  )
}

draw_cell <- function(
  x, y, r,
  fill,
  halo_alpha = 0.17,
  nucleus_alpha = 0.35
) {
  # soft halo
  grid.circle(
    x, y,
    r * 1.65,
    gp = gpar(
      fill = adjustcolor(fill, halo_alpha),
      col = NA
    )
  )

  # outer cell
  grid.circle(
    x, y,
    r,
    gp = gpar(
      fill = adjustcolor(fill, 0.86),
      col = adjustcolor(fill, 0.95),
      lwd = 0.8
    )
  )

  # pale inner body
  grid.circle(
    x, y,
    r * 0.68,
    gp = gpar(
      fill = adjustcolor("#FFFFFF", 0.17),
      col = adjustcolor("#FFFFFF", 0.23),
      lwd = 0.5
    )
  )

  # nucleus
  grid.circle(
    x - r * 0.10,
    y + r * 0.05,
    r * 0.28,
    gp = gpar(
      fill = adjustcolor(fill, nucleus_alpha),
      col = NA
    )
  )
}

# ── palette ───────────────────────────────────────────────────────────────────

deep_teal <- "#006C7E"
teal      <- "#20B9B4"
aqua      <- "#64D8D0"
blue      <- "#4D8FD8"
lavender  <- "#8C6BD6"
plum      <- "#A947B5"
coral     <- "#FF6F7D"
offwhite  <- "#F7FBFA"
white     <- "#FFFFFF"

# Viridis for trajectory and cells
viridis_pal <- colorRampPalette(
  c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")
)(200)

# Teal-coral for gauge only
gauge_palette <- colorRampPalette(
  c(teal, aqua, blue, lavender, plum, coral)
)(200)

# ── trajectory geometry ───────────────────────────────────────────────────────

traj_x <- function(t) {
  0.15 + 0.70 * t
}

traj_y <- function(t) {
  0.70 -
    0.030 * t +
    0.050 * sin(2 * pi * t)
}

# Smooth trajectory spine
ts <- seq(0, 1, length.out = 400)
sx <- traj_x(ts)
sy <- traj_y(ts)

# Main cells
n_main <- 52

pt_main <- sort(
  c(
    runif(round(n_main * 0.72), 0.02, 0.78),
    runif(round(n_main * 0.28), 0.40, 0.90)
  )
)

cx_main <- traj_x(pt_main) + rnorm(length(pt_main), 0, 0.025)
cy_main <- traj_y(pt_main) + rnorm(length(pt_main), 0, 0.033)

cx_main <- pmax(0.10, pmin(0.83, cx_main))
cy_main <- pmax(0.57, pmin(0.84, cy_main))

cr_main <- 0.009 +
  0.007 * runif(length(pt_main)) +
  0.004 * exp(-((pt_main - 0.45) / 0.25)^2)

cr_main <- pmax(0.006, pmin(0.021, cr_main))

col_main <- viridis_pal[
  pmax(1, pmin(200, round(pt_main * 199) + 1))
]

# Dense terminal cluster
n_terminal <- 30

terminal_center_x <- traj_x(0.89)
terminal_center_y <- traj_y(0.89) + 0.015

# Generate extra candidates so filtering still leaves enough cells
n_candidates <- n_terminal * 5

pt_candidates <- runif(n_candidates, 0.86, 1)

cx_candidates <- terminal_center_x +
  rnorm(n_candidates, 0, 0.034)

cy_candidates <- terminal_center_y +
  rnorm(n_candidates, 0, 0.030)

# Keep cells within a compact elliptical cluster
ellipse_distance <- sqrt(
  ((cx_candidates - terminal_center_x) / 0.070)^2 +
  ((cy_candidates - terminal_center_y) / 0.060)^2
)

keep <- which(ellipse_distance <= 1)

# Keep exactly n_terminal cells
keep <- keep[seq_len(min(n_terminal, length(keep)))]

pt_terminal <- pt_candidates[keep]
cx_terminal <- cx_candidates[keep]
cy_terminal <- cy_candidates[keep]

cr_terminal <- pmax(
  0.006,
  pmin(
    0.018,
    rnorm(length(keep), 0.011, 0.0025)
  )
)

col_terminal <- viridis_pal[
  pmax(
    1,
    pmin(
      200,
      round(pt_terminal * 199) + 1
    )
  )
]

# ── gauge geometry ────────────────────────────────────────────────────────────

gx <- 0.50
gy <- 0.330

gauge_outer <- 0.195
gauge_inner <- 0.128

# ── create transparent icon PNG ──────────────────────────────────────────────

icon_file <- tempfile(fileext = ".png")

png(
  filename = icon_file,
  width = 900,
  height = 760,
  res = 180,
  bg = "transparent"
)

grid.newpage()

# ── viridis trajectory spine ──────────────────────────────────────────────────
for (i in seq_len(length(ts) - 1)) {
  col_i <- viridis_pal[
    round(
      1 +
        (i - 1) /
        (length(ts) - 2) *
        (length(viridis_pal) - 1)
    )
  ]

  grid.lines(
    x = sx[i:(i + 1)],
    y = sy[i:(i + 1)],
    gp = gpar(
      col = adjustcolor(col_i, 0.42),
      lwd = 1.6,
      lineend = "round"
    )
  )
}

# ── main cells ────────────────────────────────────────────────────────────────

draw_order <- order(cr_main)

for (i in draw_order) {
  draw_cell(
    x = cx_main[i],
    y = cy_main[i],
    r = cr_main[i],
    fill = col_main[i]
  )
}

# ── terminal cluster glow ────────────────────────────────────────────────────

draw_order_terminal <- order(cr_terminal)

for (i in draw_order_terminal) {
  draw_cell(
    x = cx_terminal[i],
    y = cy_terminal[i],
    r = cr_terminal[i],
    fill = col_terminal[i],
    halo_alpha = 0.20,
    nucleus_alpha = 0.34
  )
}

# ── white arrowhead ───────────────────────────────────────────────────────────

t_arrow <- seq(0.82, 1, length.out = 80)

grid.lines(
  x = traj_x(t_arrow),
  y = traj_y(t_arrow),
  arrow = arrow(
    type = "closed",
    length = unit(0.12, "inches")
  ),
  gp = gpar(
    col = adjustcolor(white, 0.75),
    lwd = 2.8,
    lineend = "round",
    linejoin = "round"
  )
)

# ── gauge track ───────────────────────────────────────────────────────────────

gauge_track <- ring_polygon(
  gx,
  gy,
  gauge_outer,
  gauge_inner,
  0,
  180
)

grid.polygon(
  gauge_track$x,
  gauge_track$y,
  gp = gpar(
    fill = adjustcolor(deep_teal, 0.30),
    col = adjustcolor(white, 0.45),
    lwd = 1.2
  )
)

# Gauge gradient segments
gauge_angles <- seq(0, 180, length.out = 41)

for (i in seq_len(length(gauge_angles) - 1)) {
  segment <- ring_polygon(
    gx,
    gy,
    gauge_outer,
    gauge_inner,
    gauge_angles[i],
    gauge_angles[i + 1],
    n = 8
  )

  segment_col <- gauge_palette[
    round(
      1 +
        (i - 1) /
        (length(gauge_angles) - 2) *
        (length(gauge_palette) - 1)
    )
  ]

  grid.polygon(
    segment$x,
    segment$y,
    gp = gpar(
      fill = segment_col,
      col = NA
    )
  )
}

# Gauge tick marks
tick_angles <- seq(10, 170, by = 20)

for (deg in tick_angles) {
  rad <- deg * pi / 180

  is_major <- deg %in% c(10, 90, 170)

  tick_outer <- gauge_outer - 0.012
  tick_inner <- tick_outer - if (is_major) 0.026 else 0.016

  grid.lines(
    x = c(
      gx + tick_outer * cos(rad),
      gx + tick_inner * cos(rad)
    ),
    y = c(
      gy + tick_outer * sin(rad),
      gy + tick_inner * sin(rad)
    ),
    gp = gpar(
      col = white,
      lwd = if (is_major) 2.1 else 1.35,
      lineend = "round"
    )
  )
}

# ── DOE label inside gauge ────────────────────────────────────────────────────

grid.text(
  "DOE",
  x = gx,
  y = gy + 0.072,
  gp = gpar(
    fontfamily = "mont",
    fontface   = "bold",
    fontsize   = 16,
    col        = deep_teal
  )
)

# ── gauge needle ──────────────────────────────────────────────────────────────

needle_deg    <- 42
needle_rad    <- needle_deg * pi / 180
needle_length <- gauge_inner * 0.95

grid.lines(
  x = c(gx, gx + needle_length * cos(needle_rad)),
  y = c(gy, gy + needle_length * sin(needle_rad)),
  gp = gpar(col = white, lwd = 5, lineend = "round")
)

grid.lines(
  x = c(gx, gx + needle_length * cos(needle_rad)),
  y = c(gy, gy + needle_length * sin(needle_rad)),
  gp = gpar(col = deep_teal, lwd = 3.2, lineend = "round")
)

grid.circle(
  gx, gy, 0.018,
  gp = gpar(fill = deep_teal, col = white, lwd = 2)
)

grid.circle(
  gx, gy, 0.007,
  gp = gpar(fill = white, col = NA)
)


# ── wordmark inside icon (purple X) ──────────────────────────────────────────

grid.text(
  "BioTraj",
  x    = 0.695,
  y    = 0.175,
  just = "right",
  gp   = gpar(
    fontfamily = "mont",
    fontface   = "bold",
    fontsize   = 48,
    col        = deep_teal
  )
)

grid.text(
  "X",
  x    = 0.700,
  y    = 0.175,
  just = "left",
  gp   = gpar(
    fontfamily = "mont",
    fontface   = "bold",
    fontsize   = 48,
    col        = "#8C3DCC"
  )
)

dev.off()

# ── final hex sticker ─────────────────────────────────────────────────────────

sticker(
  subplot = icon_file,

  package = "",
  p_size  = 1,
  p_color = adjustcolor(white, 0),

  s_x     = 1.00,
  s_y     = 1.08,
  s_width = 0.95,
  s_height = 0.84,

  h_fill  = "#EAF8F6",
  h_color = teal,
  h_size  = 1.8,

  filename = "inst/BioTrajX_hex_teal_coral_DOE.png",
  dpi = 900
)

sticker(
  subplot = icon_file,

  package = "",
  p_size  = 1,
  p_color = adjustcolor(white, 0),

  s_x     = 1.00,
  s_y     = 1.08,
  s_width = 0.95,
  s_height = 0.84,

  h_fill  = "#EAF8F6",
  h_color = teal,
  h_size  = 1.8,

  filename = "inst/BioTrajX_hex_teal_coral_DOE.svg",
  dpi = 900
)
