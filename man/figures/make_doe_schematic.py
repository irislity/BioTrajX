#!/usr/bin/env python3
"""
make_doe_schematic.py

Recreates the D / O / E metric schematic with three card panels.
Embeds example plots from manuscript/figures/synthetic/ (gaussian_linear
scenario) and uses the current DOE-score formula design:

  D  -- ReLU Spearman:  max(0, ±ρ_S)                 [was (1+ρ)/2]
  O  -- null-calibrated isotonic R²                   [was raw isotonic R²]
  E  -- precision@k for "early" cells  (not "naive")  [naming change]
         + E_comp = harmonic mean(E_early, E_term)

Font style matches the original schematic:
  - Text: Helvetica (Arial equivalent), via usetex + helvet package
  - Math: Computer Modern (LaTeX default serif math font)
  - Titles: all-black bold; equations left-aligned; inline \textbf{} in desc

Output:
  manuscript/figures/schematics/doe_schematic.pdf
  manuscript/figures/schematics/doe_schematic.png
"""

from pathlib import Path

import numpy as np
import pymupdf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

matplotlib.rcParams.update({
    "font.family":      "Arial",
    "font.sans-serif":  ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "cm",   # Computer Modern — identical to LaTeX's math font
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO    = Path(__file__).resolve().parents[2]
SYN     = REPO / "figures" / "synthetic"
OUT_DIR = REPO / "figures" / "schematics"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _pdf_to_array(path: Path, dpi: int = 220) -> np.ndarray:
    doc = pymupdf.open(str(path))
    pix = doc[0].get_pixmap(matrix=pymupdf.Matrix(dpi / 72, dpi / 72), alpha=False)
    return np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3)


# ---------------------------------------------------------------------------
# Card definitions
# ---------------------------------------------------------------------------

GRAY_CARD   = "none"          # transparent background
GRAY_BORDER = "#bbbbbb"
GRAY_TEXT   = "#222222"

# Border colors matched to panel_A.svg palette
CARD_COLORS = {
    "D": "#a2c4c9",   # teal
    "O": "#f4c542",   # amber
    "E": "#93c47d",   # green
}

CARDS = [
    # ── D : Directionality ──────────────────────────────────────────────
    dict(
        title  = r"\textbf{D}: Directionality",
        border = CARD_COLORS["D"],
        pdf    = SYN / "01_gaussian_linear_c_D_metric.pdf",
        formulas = [
            r"$D_{\mathrm{early}} := \max(0,\; -\rho_S(\mathbf{t},\, s_{\mathrm{early}})) \in [0,1],$",
            r"$D_{\mathrm{term}}\;\;\, := \max(0,\;\; \rho_S(\mathbf{t},\, s_{\mathrm{term}}))\;\, \in [0,1].$",
            r"$D_{\mathrm{comp}} := \frac{2}{1/D_{\mathrm{early}} + 1/D_{\mathrm{term}}}$",
        ],
        desc = [
            "Based on Spearman correlation between",
            "pseudotime and early or terminal marker scores.",
        ],
    ),
    # ── O : Order consistency ───────────────────────────────────────────
    dict(
        title  = r"\textbf{O}: Order consistency",
        border = CARD_COLORS["O"],
        pdf    = SYN / "01_gaussian_linear_d_O_metric.pdf",
        formulas = [
            r"$O = \max\!\left(0,\; \frac{R^2 - R^2_{\mathrm{null}}}{1 - R^2_{\mathrm{null}}}\right) \in [0,1],$",
            r"$R^2:\ \text{isotonic regression on pseudobulk bin means}$",
        ],
        desc = [
            "Tests monotonic trend of marker expression",
            r"across pseudotime via null-adjusted $R^2$.",
        ],
    ),
    # ── E : Endpoint validity ───────────────────────────────────────────
    dict(
        title  = r"\textbf{E}: Endpoint validity",
        border = CARD_COLORS["E"],
        pdf    = SYN / "01_gaussian_linear_e_E_metric.pdf",
        formulas = [
            r"$E_{\mathrm{early}} := \frac{\#\{\text{early cells in top-}k_N\}}{k_N} \in [0,1]$",
            r"$E_{\mathrm{term}}\;\;\, := \frac{\#\{\text{terminal cells in top-}k_T\}}{k_T} \in [0,1]$",
            r"$E_{\mathrm{comp}} := \frac{2}{1/E_{\mathrm{early}} + 1/E_{\mathrm{term}}}$",
        ],
        desc = [
            "Precision@k for first/last k cells;",
            "composite endpoint score is harmonic mean.",
        ],
    ),
]


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

def _draw_card(
    fig: plt.Figure,
    card: dict,
    x0: float, y0: float,
    cw: float, ch: float,
) -> None:
    ax = fig.add_axes([x0, y0, cw, ch])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Transparent card with colored border
    border_color = card.get("border", GRAY_BORDER)
    ax.add_patch(mpatches.FancyBboxPatch(
        (0.01, 0.01), 0.98, 0.98,
        boxstyle="round,pad=0.025",
        facecolor=GRAY_CARD, edgecolor=border_color,
        linewidth=2.0, transform=ax.transAxes, clip_on=False, zorder=0,
    ))

    # ── Title ────────────────────────────────────────────────────────
    letter = card["title"].split(":")[0].replace(r"\textbf{", "").replace("}", "")
    rest   = ":" + card["title"].split(":", 1)[1]
    ax.text(0.5, 0.978, f"{letter}{rest}",
            ha="center", va="top", fontsize=16.5, fontweight="bold",
            color=GRAY_TEXT, transform=ax.transAxes)

    # ── Embedded plot image ───────────────────────────────────────────
    img = _pdf_to_array(card["pdf"])
    if "Directionality" in card["title"]:
        inset_rect = [0.04, 0.43, 0.92, 0.48]
    else:
        inset_rect = [0.04, 0.40, 0.92, 0.51]
    plot_ax = ax.inset_axes(inset_rect)
    plot_ax.imshow(img, aspect="auto", interpolation="bicubic")
    plot_ax.axis("off")

    # ── Formulas ─────────────────────────────────────────────────────
    n_form = len(card["formulas"])
    if n_form <= 2:
        form_top  = 0.385
        form_step = 0.078
        formula_fs = 14.8
        desc_gap   = 0.080
    else:
        form_top  = 0.390
        form_step = 0.072
        formula_fs = 14.2
        desc_gap   = 0.065

    for i, formula in enumerate(card["formulas"]):
        y = form_top - i * form_step
        if "Order" in card["title"] and i == 1:
            y -= 0.040
        ax.text(0.5, y, formula,
                ha="center", va="top", fontsize=formula_fs,
                color=GRAY_TEXT, transform=ax.transAxes)

    # ── Description ──────────────────────────────────────────────────
    last_form_y = form_top - (n_form - 1) * form_step
    if "Order" in card["title"]:
        last_form_y -= 0.040
        desc_gap    += 0.035
    if "Endpoint" in card["title"]:
        desc_gap    += 0.030
    desc_top  = last_form_y - desc_gap
    desc_step = 0.048

    for i, text in enumerate(card["desc"][:2]):
        ax.text(0.08, desc_top - i * desc_step, text,
                ha="left", va="top", fontsize=14.2,
                fontweight="normal", color=GRAY_TEXT, transform=ax.transAxes)


def make_figure() -> plt.Figure:
    fig = plt.figure(figsize=(18, 6.3))
    fig.patch.set_facecolor("white")

    card_w   = 0.295
    card_h   = 0.80
    card_y   = 0.13
    gap      = 0.025
    x_starts = [0.02, 0.02 + card_w + gap, 0.02 + 2 * (card_w + gap)]

    for card, x0 in zip(CARDS, x_starts):
        _draw_card(fig, card, x0, card_y, card_w, card_h)

    # Composite DOE formula centred below the cards
    fig.text(
        0.5, 0.065,
        r"$\mathrm{DOE\ score} = \dfrac{1}{3}"
        r"\!\left(D_{\mathrm{comp}} + O + E_{\mathrm{comp}}\right)\;\in [0,1]$",
        ha="center", va="center", fontsize=15.5, color=GRAY_TEXT,
    )

    return fig


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    fig = make_figure()
    pdf_out = OUT_DIR / "doe_schematic.pdf"
    png_out = OUT_DIR / "doe_schematic.png"
    fig.savefig(pdf_out, bbox_inches="tight", dpi=150)
    fig.savefig(png_out, bbox_inches="tight", dpi=150)
    print(f"PDF → {pdf_out}")
    print(f"PNG → {png_out}")
