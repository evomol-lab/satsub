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


def _linear_trend(x, y) -> tuple[np.ndarray, np.ndarray] | None:
    """Ordinary-least-squares trend line through (x, y).

    Returns two endpoints spanning the observed x-range, or None when there
    are too few points (<2) or all x-values are identical (no defined slope).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    if x.size < 2 or np.ptp(x) == 0:
        return None
    slope, intercept = np.polyfit(x, y, 1)
    x_line = np.array([x.min(), x.max()])
    return x_line, slope * x_line + intercept


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
        trend = _linear_trend(valid[distance_col], valid[y_col])
        if trend is not None:
            x_line, y_line = trend
            fig.add_trace(
                go.Scatter(
                    x=x_line,
                    y=y_line,
                    mode="lines",
                    line=dict(color=color, width=2, dash="dash"),
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


def identity_heatmap(matrix: pd.DataFrame, title: str = "Pairwise identity (%)") -> go.Figure:
    fig = go.Figure(
        data=go.Heatmap(
            z=matrix.values,
            x=list(matrix.columns),
            y=list(matrix.index),
            colorscale=SEQUENTIAL_BLUE,
            colorbar=dict(title="% identity"),
            hovertemplate="%{y} vs %{x}<br>Identity: %{z:.2f}%<extra></extra>",
        )
    )
    fig.update_layout(title=title, **_BASE_LAYOUT)
    fig.update_yaxes(autorange="reversed")
    return fig


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
