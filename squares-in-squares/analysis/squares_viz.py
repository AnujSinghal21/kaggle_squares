"""Helpers for inspecting a Squares in Squares submission: per-n scores, charts, SVG packings."""
import glob
import math
import os
import random

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SUBMISSIONS = os.path.join(os.path.dirname(HERE), "submissions")
OUT = os.path.join(HERE, "out")

# chart colours from the reference palette (light / dark text on a transparent background)
BLUE, ORANGE, GRID, INK = "#2a78d6", "#eb6834", "#d9d8d3", "#52514e"

# ~24 solid colours for the SVG squares
SVG_COLORS = [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948",
    "#0f9ed5", "#a44a3f", "#7a9a01", "#c1559f", "#3b5bdb", "#f08c00", "#099268", "#d6336c",
    "#5f3dc4", "#e8590c", "#1971c2", "#2f9e44", "#9c36b5", "#c92a2a", "#0c8599", "#f59f00",
]


def latest_submission(path=None):
    """The given path, else the newest submissions/*.csv."""
    if path:
        return path
    files = glob.glob(os.path.join(SUBMISSIONS, "*.csv"))
    if not files:
        raise FileNotFoundError(f"no submissions in {SUBMISSIONS}")
    return max(files, key=os.path.getmtime)


def load_submission(path=None):
    """-> (path, {n: float array of shape (n, 3) with columns x, y, angle})."""
    path = latest_submission(path)
    df = pd.read_csv(path)  # float64, like Kaggle's scorer
    packings = {}
    for n, g in df.groupby("n"):
        g = g.sort_values("square_id")
        if len(g) == n:
            packings[int(n)] = g[["x", "y", "angle"]].to_numpy(float)
    return path, packings


def corners(p):
    """(n, 4, 2) corners of every unit square in packing p."""
    c, s = np.cos(p[:, 2]), np.sin(p[:, 2])
    off = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]])
    x = p[:, 0, None] + off[None, :, 0] * c[:, None] - off[None, :, 1] * s[:, None]
    y = p[:, 1, None] + off[None, :, 0] * s[:, None] + off[None, :, 1] * c[:, None]
    return np.stack([x, y], axis=-1)


def side(p):
    """Side of the tightest axis-aligned container (what the metric measures)."""
    pts = corners(p).reshape(-1, 2)
    return float(max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1])))


def per_n_table(packings, total_puzzles=1000):
    """One row per puzzle n, sorted by score, worst (largest) first. Missing n are charged the grid."""
    rows = []
    for n in range(1, total_puzzles + 1):
        grid_L = math.ceil(math.sqrt(n))
        if n in packings:
            p = packings[n]
            L = side(p)
            tilted = float(np.mean(np.abs(np.sin(2 * p[:, 2])) > 1e-6))
        else:
            L, tilted = float(grid_L), 0.0
        rows.append(dict(n=n, L=L, score=L / math.sqrt(n) - 1, grid_score=grid_L / math.sqrt(n) - 1,
                         tilted_frac=tilted, present=n in packings))
    df = pd.DataFrame(rows)
    df["gain_vs_grid"] = df["grid_score"] - df["score"]
    df["share_of_total"] = df["score"] / df["score"].sum()
    return df.sort_values("score", ascending=False).reset_index(drop=True)


def _style(ax, title, xlabel, ylabel):
    ax.set_title(title, loc="left", fontsize=11, color=INK)
    ax.set_xlabel(xlabel, color=INK)
    ax.set_ylabel(ylabel, color=INK)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=INK)


def plot_per_n(df):
    """Score of every puzzle n (the sawtooth), with the grid score for reference."""
    import matplotlib.pyplot as plt
    d = df.sort_values("n")
    fig, ax = plt.subplots(figsize=(11, 3.8))
    ax.plot(d["n"], d["grid_score"], color=GRID, lw=1.2, label="grid baseline")
    ax.plot(d["n"], d["score"], color=BLUE, lw=1.2, label="this submission")
    _style(ax, f"Score per puzzle, L/√n − 1 (total {d['score'].sum():.4f}, lower is better)", "n", "score")
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    return fig


def plot_buckets(df, width=50):
    """Total score, remaining headroom vs the grid, and mean score per bucket of `width` values of n."""
    import matplotlib.pyplot as plt
    d = df.copy()
    d["bucket"] = (d["n"] - 1) // width
    g = d.groupby("bucket").agg(total=("score", "sum"), grid_total=("grid_score", "sum"),
                                mean=("score", "mean"), gain=("gain_vs_grid", "sum"))
    labels = [f"{b * width + 1}–{(b + 1) * width}" for b in g.index]
    x = np.arange(len(g))
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    axes[0].bar(x, g["total"], color=BLUE, width=0.7)
    _style(axes[0], "Total score per bucket of n", f"n bucket (width {width})", "sum of score")
    axes[1].bar(x, g["mean"], color=BLUE, width=0.7)
    _style(axes[1], "Mean score per puzzle", f"n bucket (width {width})", "mean score")
    axes[2].bar(x, g["gain"], color=ORANGE, width=0.7)
    _style(axes[2], "Gain over the grid so far", f"n bucket (width {width})", "grid score − score")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    fig.tight_layout()
    return fig, g.assign(label=labels)


def svg_packing(n, packings, size=640, seed=None, save=True):
    """SVG of packing n: random solid colours per square, container box with its side printed."""
    if n not in packings:
        raise KeyError(f"n={n} is not in this submission")
    p = packings[n]
    pts = corners(p).reshape(-1, 2)
    x0, y0 = pts[:, 0].min(), pts[:, 1].min()
    L = side(p)
    pad, head = 12, 26
    scale = (size - 2 * pad) / L
    rng = random.Random(n if seed is None else seed)

    def px(x, y):  # y up -> y down
        return pad + (x - x0) * scale, head + pad + (L - (y - y0)) * scale

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size + head}" '
             f'viewBox="0 0 {size} {size + head}" font-family="sans-serif">',
             f'<rect width="100%" height="100%" fill="#fcfcfb"/>']
    for quad in corners(p):
        poly = " ".join("%.3f,%.3f" % px(x, y) for x, y in quad)
        parts.append(f'<polygon points="{poly}" fill="{rng.choice(SVG_COLORS)}" stroke="#fcfcfb" stroke-width="0.8"/>')
    bx, by = px(x0, y0 + L)
    parts.append(f'<rect x="{bx:.2f}" y="{by:.2f}" width="{L * scale:.2f}" height="{L * scale:.2f}" '
                 f'fill="none" stroke="#0b0b0b" stroke-width="1.5"/>')
    score = L / math.sqrt(n) - 1
    parts.append(f'<text x="{pad}" y="17" font-size="14" fill="#0b0b0b">n = {n}   side L = {L:.9f}   '
                 f'score = {score:.6f}   (grid {math.ceil(math.sqrt(n))})</text>')
    parts.append("</svg>")
    svg = "\n".join(parts)
    if save:
        os.makedirs(OUT, exist_ok=True)
        with open(os.path.join(OUT, f"n_{n:04d}.svg"), "w") as f:
            f.write(svg)
    return svg


# ---------------------------------------------------------------------------------------------
# Headroom: how far each puzzle is from the best packing any team has submitted.
# The host's "Per-Puzzle Leaderboard" notebook (Kaggle, inversion) is the source. Its data sits in a
# private dataset, so refresh it by downloading that notebook's .ipynb (with outputs) and calling
# import_leaderboard_ipynb(). It refreshes roughly once a day on Kaggle's side.
# ---------------------------------------------------------------------------------------------
LEADERBOARD_CSV = os.path.join(HERE, "data", "per_puzzle_leaderboard.csv")
DEFAULT_LB_IPYNB = os.path.expanduser("~/Downloads/squares-in-squares-per-puzzle-leaderboard.ipynb")


def import_leaderboard_ipynb(ipynb_path=DEFAULT_LB_IPYNB, out_csv=LEADERBOARD_CSV):
    """Extract the per-puzzle rows (n, L, team, date) embedded in the downloaded notebook and save them as CSV."""
    import json
    import re
    with open(ipynb_path) as f:
        nb = json.load(f)
    for cell in nb["cells"]:
        for out in cell.get("outputs", []):
            html = out.get("data", {}).get("text/html")
            if not html:
                continue
            html = "".join(html) if isinstance(html, list) else html
            m = re.search(r'(\[\{"n":\s*\d+.*?\}\])', html, re.S)
            if m:
                rows = json.loads(m.group(1))
                os.makedirs(os.path.dirname(out_csv), exist_ok=True)
                pd.DataFrame(rows).rename(columns={"L": "best_L_str"}).to_csv(out_csv, index=False)
                return len(rows)
    raise ValueError(f"no per-puzzle rows found in {ipynb_path}")


def headroom_table(df, leaderboard_csv=LEADERBOARD_CSV):
    """Join our per-n table with the best score any team holds for each n.

    headroom = our score - best score (>= 0 means a team already has a smaller box than ours).
    Sorted by headroom, largest first. Rows with team '-' still hold the grid baseline (nobody beat it).
    """
    lb = pd.read_csv(leaderboard_csv, dtype={"best_L_str": str})
    lb["best_L"] = lb["best_L_str"].astype(float)
    lb["best_score"] = lb["best_L"] / np.sqrt(lb["n"]) - 1
    d = df.merge(lb[["n", "best_L", "best_score", "team", "date"]], on="n")
    d["headroom"] = d["score"] - d["best_score"]
    d["headroom_L"] = d["L"] - d["best_L"]
    d["headroom_share"] = d["headroom"].clip(lower=0) / d["headroom"].clip(lower=0).sum()
    return d.sort_values("headroom", ascending=False).reset_index(drop=True)


def plot_headroom(d, width=50):
    """Headroom per n, per bucket of n, and how concentrated it is."""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    s = d.sort_values("n")
    axes[0].bar(s["n"], s["headroom"].clip(lower=0), color=BLUE, width=1.0)
    _style(axes[0], "Headroom per puzzle", "n", "our score − best team's score")
    b = ((s["n"] - 1) // width)
    g = s.groupby(b)["headroom"].apply(lambda x: x.clip(lower=0).sum())
    axes[1].bar(np.arange(len(g)), g.values, color=BLUE, width=0.7)
    axes[1].set_xticks(np.arange(len(g)))
    axes[1].set_xticklabels([f"{i * width + 1}–{(i + 1) * width}" for i in g.index], rotation=60, ha="right", fontsize=7)
    _style(axes[1], f"Headroom per bucket of {width} n", "n bucket", "sum of headroom")
    h = d["headroom"].clip(lower=0).to_numpy()  # already sorted descending
    axes[2].plot(np.arange(1, len(h) + 1), np.cumsum(h) / h.sum(), color=BLUE, lw=1.6)
    _style(axes[2], "How concentrated: share of total headroom", "worst puzzles first (count)", "cumulative share")
    axes[2].set_ylim(0, 1.02)
    fig.tight_layout()
    return fig
