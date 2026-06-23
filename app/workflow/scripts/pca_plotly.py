#!/usr/bin/env python3
"""Interactive (Plotly) PCA plots from a VST-transformed count matrix.

Produces:
  * <outdir>/pca.html                  - all samples
  * <outdir>/pca_<control>.html        - one per unique control condition,
                                         containing that control plus its
                                         treatment condition(s) from contrasts

PCA matches DESeq2::plotPCA: it is computed on the top-N most variable genes
(default 500), recomputed independently for each plot's sample subset.
"""

import argparse
import os

import numpy as np
import pandas as pd
import plotly.express as px
from PIL import Image


def compute_pca(vst, ntop=500):
    """Return (PC coords indexed by sample, percent variance for PC1/PC2).

    `vst` is a genes x samples DataFrame.
    """
    variances = vst.var(axis=1)
    ntop = min(ntop, len(variances))
    top = variances.sort_values(ascending=False).head(ntop).index
    mat = vst.loc[top].T  # samples x genes

    # Center each gene, then PCA via SVD (equivalent to prcomp, as plotPCA uses).
    x = mat.values - mat.values.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    coords = u * s
    pct = (s ** 2) / np.sum(s ** 2) * 100.0

    n = min(2, coords.shape[1])
    df = pd.DataFrame(coords[:, :n], index=mat.index,
                      columns=[f"PC{i + 1}" for i in range(n)])
    if n < 2:  # degenerate (single sample): pad so plotting still works
        df["PC2"] = 0.0
        pct = np.append(pct, 0.0)
    return df, pct[:2]


# PNG sizing: render on a natural-sized canvas (so fonts/markers look the same
# as the HTML) and upscale via `scale` for resolution, rather than rendering
# onto a huge canvas where everything ends up tiny. 7x6 in at 300 dpi.
PNG_DPI = 300
PNG_BASE_DPI = 96  # plotly's logical px-per-inch
PNG_WIDTH_PX = 7 * PNG_BASE_DPI
PNG_HEIGHT_PX = 6 * PNG_BASE_DPI
PNG_SCALE = PNG_DPI / PNG_BASE_DPI  # -> 2100x1800 px output


def make_plot(vst, design, samples, title, html_path):
    """Compute and write an interactive PCA HTML plus a 300-dpi PNG.

    The PNG is written alongside `html_path` with a `.png` extension.
    """
    sub = vst.loc[:, samples]
    coords, pct = compute_pca(sub)

    coords = coords.copy()
    coords["sample"] = coords.index
    coords["condition"] = design.reindex(coords.index)["condition"].astype(str)

    fig = px.scatter(
        coords, x="PC1", y="PC2",
        color="condition", hover_name="sample",
        title=title,
        labels={
            "PC1": f"PC1 ({pct[0]:.1f}%)",
            "PC2": f"PC2 ({pct[1]:.1f}%)",
        },
    )
    fig.update_traces(
        marker=dict(size=12, line=dict(width=1, color="DarkSlateGrey"))
    )
    fig.update_layout(template="plotly_white", legend_title_text="condition")

    os.makedirs(os.path.dirname(os.path.abspath(html_path)), exist_ok=True)
    # Self-contained HTML so it renders offline, plus a static 300-dpi PNG.
    fig.write_html(html_path, include_plotlyjs=True)

    png_path = os.path.splitext(html_path)[0] + ".png"
    fig.write_image(png_path, width=PNG_WIDTH_PX, height=PNG_HEIGHT_PX,
                    scale=PNG_SCALE)
    # write_image only sets pixel dimensions; stamp the DPI metadata so the
    # file reports a true 7x6 in at 300 dpi.
    with Image.open(png_path) as img:
        img.save(png_path, dpi=(PNG_DPI, PNG_DPI))
    print(f"[pca] wrote {html_path} and {png_path} ({len(samples)} samples)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vst", required=True,
                    help="VST matrix CSV (genes x samples)")
    ap.add_argument("--design", required=True,
                    help="Design CSV with columns: sample,condition")
    ap.add_argument("--contrasts", required=True,
                    help="Contrasts CSV with columns: control,treatment")
    ap.add_argument("--outdir", required=True,
                    help="Output directory for the HTML plots")
    args = ap.parse_args()

    vst = pd.read_csv(args.vst, index_col=0)
    design = pd.read_csv(args.design).set_index("sample")
    contrasts = pd.read_csv(args.contrasts)

    os.makedirs(args.outdir, exist_ok=True)

    # Overall PCA across every sample.
    make_plot(vst, design, list(vst.columns),
              "PCA (VST-transformed counts)",
              os.path.join(args.outdir, "pca.html"))

    # One PCA per unique control condition: the control plus its treatment(s).
    for ctrl in sorted(contrasts["control"].unique()):
        treatments = (
            contrasts.loc[contrasts["control"] == ctrl, "treatment"]
            .unique()
            .tolist()
        )
        conditions = [ctrl] + [t for t in treatments if t != ctrl]
        in_group = design.index[design["condition"].isin(conditions)]
        samples = [s for s in vst.columns if s in set(in_group)]

        title = f"PCA: {ctrl} vs {', '.join(str(t) for t in treatments)}"
        make_plot(vst, design, samples, title,
                  os.path.join(args.outdir, f"pca_{ctrl}.html"))


if __name__ == "__main__":
    main()
