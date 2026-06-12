# molecular_calculator/ui/components/gmm_charts.py
"""Plotly figures for the GMM Analysis page.

Pure figure builders (no Streamlit), so they are unit-testable. Colors use the
ColorBrewer Set2 qualitative palette, matching the Impulator GMM charts.
"""

import numpy as np
import plotly.express as px
import plotly.graph_objects as go

_PALETTE = px.colors.qualitative.Set2


def create_density_overlay(values_1d, analysis, column_name: str) -> go.Figure:
    """Histogram of one property with each group's fitted Gaussian overlaid.

    Caller must ensure ``values_1d`` has non-zero variance (the upstream
    ``gmm_sentinel_check`` guarantees this in the normal page flow).
    """
    v = np.asarray(values_1d, dtype=float).ravel()
    grid = np.linspace(float(v.min()), float(v.max()), 200)
    means, _weights, _sigmas, pdfs = analysis.component_curves(grid)

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=v, histnorm="probability density", name="Your data",
        marker=dict(color="lightgray"), opacity=0.6,
    ))
    for k in range(len(means)):
        color = _PALETTE[k % len(_PALETTE)]
        fig.add_trace(go.Scatter(
            x=grid, y=pdfs[k], mode="lines", line=dict(width=2.5, color=color),
            fill="tozeroy", name=f"Group {k + 1} (avg {means[k]:.2f})",
        ))
    fig.update_layout(
        title=f"Distribution of {column_name}, split into groups",
        xaxis_title=column_name, yaxis_title="Density",
        template="plotly_white", barmode="overlay",
    )
    return fig


def create_cluster_scatter(analysis, x_idx: int, y_idx: int,
                           x_name: str, y_name: str) -> go.Figure:
    """Scatter of two properties, points colored by group, with group centers."""
    v = analysis.raw_values
    labels = analysis.labels
    centers = analysis.means_real_units

    fig = go.Figure()
    for g in range(analysis.n_components):
        m = labels == g
        color = _PALETTE[g % len(_PALETTE)]
        fig.add_trace(go.Scatter(
            x=v[m, x_idx], y=v[m, y_idx], mode="markers",
            marker=dict(color=color, size=7, opacity=0.7), name=f"Group {g + 1}",
        ))
    fig.add_trace(go.Scatter(
        x=centers[:, x_idx], y=centers[:, y_idx], mode="markers",
        marker=dict(symbol="star", size=16, color="black",
                    line=dict(width=1, color="white")),
        name="Group centers",
    ))
    fig.update_layout(
        title=f"{y_name} vs {x_name}, colored by group",
        xaxis_title=x_name, yaxis_title=y_name, template="plotly_white",
    )
    return fig


def create_bic_aic_plot(sweep_df, recommended_k: int) -> go.Figure:
    """Line plot of BIC and AIC vs number of groups, with the recommendation marked."""
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=sweep_df["n_groups"], y=sweep_df["bic"], mode="lines+markers", name="BIC",
    ))
    fig.add_trace(go.Scatter(
        x=sweep_df["n_groups"], y=sweep_df["aic"], mode="lines+markers", name="AIC",
    ))
    fig.add_vline(
        x=recommended_k, line_dash="dash", line_color="green",
        annotation_text=f"Recommended: {recommended_k} groups",
    )
    fig.update_layout(
        title="Model quality vs number of groups (lower is better)",
        xaxis_title="Number of groups", yaxis_title="Score (BIC / AIC)",
        template="plotly_white",
    )
    return fig
