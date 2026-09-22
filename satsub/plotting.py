"""Plotly figure builders for the SatSub GUI.

Colors follow a fixed categorical order (never auto-cycled) and a single-hue
sequential ramp for magnitude (identity %), per the project's chart style:
- Series 1 (transitions / primary magnitude): blue  #2a78d6
- Series 2 (transversions / secondary):       orange #eb6834
- Series 3/4 (physicochemical extra groups):  aqua #1baf7a, yellow #eda100
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

# Fixed categorical order (do not reorder / do not auto-cycle).
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
CATEGORICAL_ORDER = [BLUE, ORANGE, AQUA, YELLOW, "#e87ba4", "#008300", "#4a3aa7", "#e34948"]

SURFACE = "#fcfcfb"
PRIMARY_INK = "#0b0b0b"
SECONDARY_INK = "#52514e"
MUTED_INK = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"

SEQUENTIAL_BLUE = [
    [0.0, "#cde2fb"],
    [0.15, "#b7d3f6"],
    [0.3, "#9ec5f4"],
    [0.45, "#6da7ec"],
    [0.6, "#3987e5"],
    [0.75, "#256abf"],
    [0.9, "#184f95"],
    [1.0, "#0d366b"],
]

_BASE_LAYOUT = dict(
    plot_bgcolor=SURFACE,
    paper_bgcolor=SURFACE,
    font=dict(color=PRIMARY_INK, family="system-ui, -apple-system, 'Segoe UI', sans-serif"),
    margin=dict(l=60, r=30, t=60, b=60),
)


def _style_axes(fig: go.Figure) -> None:
    fig.update_xaxes(showgrid=True, gridcolor=GRIDLINE, zeroline=False, linecolor=BASELINE, tickfont=dict(color=MUTED_INK))
    fig.update_yaxes(showgrid=True, gridcolor=GRIDLINE, zeroline=False, linecolor=BASELINE, tickfont=dict(color=MUTED_INK))


def _lowess(x, y, frac: float = 0.6667) -> tuple[np.ndarray, np.ndarray] | None:
    """Locally weighted scatterplot smoothing (LOWESS): a degree-1 local
    regression at each observed x, weighted by a tricube kernel over its
    nearest neighbors.

    Unlike a single straight-line fit, this can bend and flatten, which is
    what a saturation plot needs to show: points tracking distance linearly
    versus points plateauing as multiple/back substitutions erase signal.

    Returns (x_sorted, y_smoothed) evaluated at every valid observation, or
    None when there are too few points (<4) or all x-values are identical.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    n = x.size
    if n < 4 or np.ptp(x) == 0:
        return None

    order = np.argsort(x)
    xs, ys = x[order], y[order]
    k = int(np.clip(np.ceil(frac * n), 3, n))
    design = np.column_stack([np.ones(n), xs])

    y_smooth = np.empty(n)
    for i in range(n):
        dist = np.abs(xs - xs[i])
        h = np.partition(dist, k - 1)[k - 1]
        if h <= 0:
            h = np.finfo(float).eps
        w = np.clip(1 - np.clip(dist / h, 0, 1) ** 3, 0, None) ** 3  # tricube weights
        wx = design * w[:, None]
        try:
            beta, *_ = np.linalg.lstsq(wx.T @ design, wx.T @ ys, rcond=None)
        except np.linalg.LinAlgError:
            y_smooth[i] = ys[i]
            continue
        y_smooth[i] = beta[0] + beta[1] * xs[i]

    return xs, y_smooth


def saturation_scatter(
    df: pd.DataFrame,
    distance_col: str,
    distance_label: str,
    title: str,
) -> tuple[go.Figure, int]:
    """Scatter of transition (Ps) and transversion (Pv) proportions against a
    genetic-distance estimate, one point per sequence pair. Pairs with an
    undefined distance (fully saturated under that model) are dropped from
    the plot; the count dropped is returned alongside the figure."""
    valid = df.dropna(subset=[distance_col])
    n_dropped = len(df) - len(valid)

    hover = (
        "%{customdata[0]} vs %{customdata[1]}"
        "<br>" + distance_label + ": %{x:.4f}"
        "<br>%{customdata[2]}: %{y:.4f}<extra></extra>"
    )

    fig = go.Figure()

    def _add_series(y_col: str, color: str, series_name: str) -> None:
        trend = _lowess(valid[distance_col], valid[y_col])
        if trend is not None:
            x_line, y_line = trend
            fig.add_trace(
                go.Scatter(
                    x=x_line,
                    y=y_line,
                    mode="lines",
                    line=dict(color=color, width=2, dash="dash", shape="spline", smoothing=0.6),
                    name=f"{series_name} trend",
                    showlegend=False,
                    hoverinfo="skip",
                )
            )
        fig.add_trace(
            go.Scatter(
                x=valid[distance_col],
                y=valid[y_col],
                mode="markers",
                name=series_name,
                marker=dict(color=color, size=9, line=dict(color=SURFACE, width=1)),
                customdata=np.stack(
                    [valid["seq1"], valid["seq2"], np.full(len(valid), series_name)], axis=-1
                ),
                hovertemplate=hover,
            )
        )

    _add_series("Ps", BLUE, "Transitions (s)")
    _add_series("Pv", ORANGE, "Transversions (v)")

    fig.update_layout(
        title=title,
        xaxis_title=f"Genetic distance ({distance_label})",
        yaxis_title="Observed proportion of substitutions",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        **_BASE_LAYOUT,
    )
    _style_axes(fig)
    return fig, n_dropped


def generic_heatmap(
    matrix: pd.DataFrame,
    title: str,
    colorbar_title: str = "Value",
    value_label: str = "Value",
    value_format: str = ".2f",
    value_suffix: str = "",
    row_label: str = "%{y}",
    col_label: str = "%{x}",
) -> go.Figure:
    """Sequential-blue heatmap for any square/rectangular magnitude matrix
    (pairwise identity, directional base-pair counts, etc.)."""
    fig = go.Figure(
        data=go.Heatmap(
            z=matrix.values,
            x=list(matrix.columns),
            y=list(matrix.index),
            colorscale=SEQUENTIAL_BLUE,
            colorbar=dict(title=colorbar_title),
            hovertemplate=(
                f"{row_label} → {col_label}<br>{value_label}: %{{z:{value_format}}}{value_suffix}<extra></extra>"
            ),
        )
    )
    fig.update_layout(title=title, **_BASE_LAYOUT)
    fig.update_yaxes(autorange="reversed")
    return fig


def identity_heatmap(matrix: pd.DataFrame, title: str = "Pairwise identity (%)") -> go.Figure:
    return generic_heatmap(
        matrix, title, colorbar_title="% identity", value_label="Identity", value_format=".2f", value_suffix="%"
    )


def directional_pair_heatmap(matrix: pd.DataFrame, title: str) -> go.Figure:
    """Heatmap of averaged directional base-pair counts: row = base in the
    earlier-listed sequence of each pair, column = base in the later one."""
    return generic_heatmap(
        matrix,
        title,
        colorbar_title="Mean count / pair",
        value_label="Mean count",
        value_format=".1f",
        row_label="row (earlier seq) %{y}",
        col_label="col (later seq) %{x}",
    )


def categorical_bar(categories: list[str], values: list[float], title: str, y_title: str) -> go.Figure:
    """Bar chart for a small set of named categories (e.g. base composition,
    physicochemical groups), colored in the fixed categorical order."""
    colors = [CATEGORICAL_ORDER[i % len(CATEGORICAL_ORDER)] for i in range(len(categories))]
    fig = go.Figure(
        data=go.Bar(
            x=categories,
            y=values,
            marker=dict(color=colors),
            hovertemplate="%{x}: %{y:.2f}<extra></extra>",
        )
    )
    fig.update_layout(title=title, yaxis_title=y_title, showlegend=False, **_BASE_LAYOUT)
    _style_axes(fig)
    return fig


def grouped_bar_by_base(
    x_categories: list[str], series: dict, title: str, y_title: str
) -> go.Figure:
    """Grouped bar chart for base-composition data split by an x-axis
    category (e.g. codon position). `series` maps a base letter (A/C/G/T)
    to a list of values aligned with `x_categories`; base colors follow the
    same fixed A/C/G/T -> categorical-slot mapping used elsewhere."""
    fig = go.Figure()
    for i, base in enumerate("ACGT"):
        if base not in series:
            continue
        fig.add_trace(
            go.Bar(
                x=x_categories,
                y=series[base],
                name=base,
                marker=dict(color=CATEGORICAL_ORDER[i]),
                hovertemplate=f"{base} - " + "%{x}: %{y:.2f}<extra></extra>",
            )
        )
    fig.update_layout(
        title=title,
        yaxis_title=y_title,
        barmode="group",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        **_BASE_LAYOUT,
    )
    _style_axes(fig)
    return fig


def rscu_diverging_bar(df: pd.DataFrame, title: str = "Relative synonymous codon usage (RSCU)") -> go.Figure:
    """RSCU bar chart, one bar per codon (grouped by amino acid, in the
    order given), colored on a blue<->red diverging scale centered at 1.0
    (the value expected under equal usage within a synonymous family)."""
    labels = [f"{row.codon} ({row.amino_acid})" for row in df.itertuples()]
    diverging_scale = [
        [0.0, "#0d366b"],
        [0.25, "#3987e5"],
        [0.5, "#f0efec"],
        [0.75, "#e88f8f"],
        [1.0, "#7a1f1f"],
    ]
    fig = go.Figure(
        data=go.Bar(
            x=labels,
            y=df["rscu"],
            marker=dict(
                color=df["rscu"],
                colorscale=diverging_scale,
                cmid=1.0,
                colorbar=dict(title="RSCU"),
                line=dict(color=SURFACE, width=0.5),
            ),
            customdata=df["count_avg"],
            hovertemplate="%{x}<br>RSCU: %{y:.2f}<br>Mean count: %{customdata:.1f}<extra></extra>",
        )
    )
    fig.add_hline(y=1.0, line=dict(color=MUTED_INK, width=1, dash="dot"))
    fig.update_layout(title=title, yaxis_title="RSCU", height=480, **_BASE_LAYOUT)
    fig.update_xaxes(tickangle=-90, tickfont=dict(size=9))
    _style_axes(fig)
    return fig


def monochrome_bar(categories: list[str], values: list[float], title: str, y_title: str) -> go.Figure:
    """Single-hue bar chart for larger category sets (e.g. 20 amino acids)
    where distinct per-bar colors would not carry meaning."""
    fig = go.Figure(
        data=go.Bar(
            x=categories,
            y=values,
            marker=dict(color=BLUE),
            hovertemplate="%{x}: %{y:.2f}<extra></extra>",
        )
    )
    fig.update_layout(title=title, yaxis_title=y_title, showlegend=False, **_BASE_LAYOUT)
    _style_axes(fig)
    return fig
