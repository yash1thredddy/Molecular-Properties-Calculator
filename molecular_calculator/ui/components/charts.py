"""Chart components for molecular data visualization.

This module provides reusable chart components built on Plotly,
with support for the structure viewer integration.

This is the CANONICAL source for chart creation functions.
All chart creation should use functions from this module.
"""

import logging
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from typing import Optional, List, Dict, Any, Tuple
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

from molecular_calculator.ui.components.molecule_viewer import (
    embed_structure_viewer,
    render_structure_viewer_hint,
    prepare_chart_customdata,
)

logger = logging.getLogger(__name__)


# =============================================================================
# Chart Creation Functions
# =============================================================================

def create_scatter_plot(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: Optional[str] = None,
    size_col: Optional[str] = None,
    shape_col: Optional[str] = None,
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None,
    title: Optional[str] = None,
    trendline: bool = False,
    color_scale: str = "Viridis",
    marker_size: int = 8,
    size_max: int = 20,
    opacity: float = 0.7,
    key_prefix: str = "scatter"
) -> Tuple[go.Figure, Optional[Dict[str, Any]]]:
    """Create a scatter plot with structure viewer support and optional trendline.

    Args:
        df: DataFrame with data
        x_col: X-axis column
        y_col: Y-axis column
        color_col: Optional column for color encoding
        size_col: Optional column for size encoding
        shape_col: Optional column for shape/symbol encoding
        smiles_col: SMILES column for structure viewer
        name_col: Name/ID column for display
        title: Chart title
        trendline: Whether to add a trendline
        color_scale: Color scale for continuous color
        marker_size: Base marker size (used when size_col is None)
        size_max: Maximum marker size (used when size_col is set)
        opacity: Marker opacity (0.0 to 1.0)
        key_prefix: Key prefix for widget uniqueness (reserved for future use)

    Returns:
        Tuple of (Plotly Figure, regression_stats dict or None)
    """
    # Note: key_prefix is reserved for future widget key generation
    _ = key_prefix

    # Prepare data - drop NaN for size column if specified
    if size_col:
        plot_df = df.dropna(subset=[size_col]).copy()
    else:
        plot_df = df.copy()

    # Build scatter parameters
    scatter_params = {
        'x': x_col,
        'y': y_col,
        'title': title or f"{y_col.replace('_', ' ')} vs {x_col.replace('_', ' ')}",
        'opacity': opacity,
    }

    if color_col:
        scatter_params['color'] = color_col
        # Check if color column is numeric for continuous scale
        numeric_cols = plot_df.select_dtypes(include=[np.number]).columns.tolist()
        if color_col in numeric_cols:
            scatter_params['color_continuous_scale'] = color_scale

    if size_col:
        scatter_params['size'] = size_col
        scatter_params['size_max'] = size_max

    if shape_col:
        scatter_params['symbol'] = shape_col

    if trendline:
        scatter_params['trendline'] = 'ols'

    # Add customdata for structure viewer if SMILES column is available
    if smiles_col and smiles_col in plot_df.columns:
        plot_df['_row_index'] = range(len(plot_df))
        if name_col and name_col in plot_df.columns:
            scatter_params['custom_data'] = [smiles_col, name_col, '_row_index']
        else:
            scatter_params['custom_data'] = [smiles_col, '_row_index']

    # Create figure
    fig = px.scatter(plot_df, **scatter_params)

    # Set marker size if not using size column
    if not size_col:
        fig.update_traces(marker=dict(size=marker_size))

    # Calculate regression statistics if trendline is shown
    regression_stats = None
    if trendline:
        regression_stats = _calculate_regression_stats(plot_df, x_col, y_col)
        if regression_stats:
            _add_regression_annotation(fig, regression_stats, x_col, y_col)

    # Update layout
    fig.update_layout(
        xaxis_title=x_col.replace('_', ' '),
        yaxis_title=y_col.replace('_', ' '),
        template="plotly_white",
        hovermode="closest"
    )

    return fig, regression_stats


def create_histogram(
    df: pd.DataFrame,
    column: str,
    bins: int = 30,
    color_col: Optional[str] = None,
    title: Optional[str] = None,
    show_mean: bool = True,
    show_median: bool = True
) -> go.Figure:
    """Create a histogram with optional mean/median lines.

    Args:
        df: DataFrame with data
        column: Column to plot
        bins: Number of bins
        color_col: Optional column for color grouping
        title: Chart title
        show_mean: Whether to show mean line
        show_median: Whether to show median line

    Returns:
        Plotly Figure
    """
    fig = px.histogram(
        df,
        x=column,
        color=color_col,
        nbins=bins,
        title=title or f"Distribution of {column.replace('_', ' ')}",
        template="plotly_white"
    )

    # Add mean line
    if show_mean:
        mean_val = df[column].mean()
        if not pd.isna(mean_val):
            fig.add_vline(
                x=mean_val,
                line_dash="dash",
                line_color="red",
                annotation_text=f"Mean: {mean_val:.2f}"
            )

    # Add median line
    if show_median:
        median_val = df[column].median()
        if not pd.isna(median_val):
            fig.add_vline(
                x=median_val,
                line_dash="dot",
                line_color="blue",
                annotation_text=f"Median: {median_val:.2f}"
            )

    fig.update_layout(
        xaxis_title=column.replace('_', ' '),
        yaxis_title="Count"
    )

    return fig


def create_box_plot(
    df: pd.DataFrame,
    y_col: str,
    x_col: Optional[str] = None,
    color_col: Optional[str] = None,
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None,
    title: Optional[str] = None,
    show_points: bool = False
) -> go.Figure:
    """Create a box plot with structure viewer support.

    Args:
        df: DataFrame with data
        y_col: Y-axis column (values)
        x_col: Optional X-axis column (categories/groups)
        color_col: Optional column for color grouping
        smiles_col: SMILES column for structure viewer
        name_col: Name/ID column for display
        title: Chart title
        show_points: Whether to show all data points

    Returns:
        Plotly Figure
    """
    plot_df = df.copy()

    # Build title
    box_title = title or f"Box Plot of {y_col.replace('_', ' ')}"
    if x_col:
        box_title += f" by {x_col.replace('_', ' ')}"

    # Build customdata for structure viewer
    custom_data_cols = None
    if smiles_col and smiles_col in plot_df.columns:
        plot_df['_row_index'] = range(len(plot_df))
        if name_col and name_col in plot_df.columns:
            custom_data_cols = [smiles_col, name_col, '_row_index']
        else:
            custom_data_cols = [smiles_col, '_row_index']

    fig = px.box(
        plot_df,
        x=x_col,
        y=y_col,
        color=color_col if color_col else x_col,
        title=box_title,
        points="all" if show_points else False,
        custom_data=custom_data_cols,
        template="plotly_white"
    )

    return fig


def create_violin_plot(
    df: pd.DataFrame,
    y_col: str,
    x_col: Optional[str] = None,
    color_col: Optional[str] = None,
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None,
    title: Optional[str] = None,
    show_points: bool = False,
    show_box: bool = True
) -> go.Figure:
    """Create a violin plot with structure viewer support.

    Args:
        df: DataFrame with data
        y_col: Y-axis column (values)
        x_col: Optional X-axis column (categories)
        color_col: Optional column for color grouping
        smiles_col: SMILES column for structure viewer
        name_col: Name/ID column for display
        title: Chart title
        show_points: Whether to show all data points
        show_box: Whether to show box plot inside violin

    Returns:
        Plotly Figure
    """
    plot_df = df.copy()

    # Build title
    violin_title = title or f"Violin Plot of {y_col.replace('_', ' ')}"
    if x_col:
        violin_title += f" by {x_col.replace('_', ' ')}"

    # Build customdata for structure viewer
    custom_data_cols = None
    if smiles_col and smiles_col in plot_df.columns:
        plot_df['_row_index'] = range(len(plot_df))
        if name_col and name_col in plot_df.columns:
            custom_data_cols = [smiles_col, name_col, '_row_index']
        else:
            custom_data_cols = [smiles_col, '_row_index']

    fig = px.violin(
        plot_df,
        x=x_col,
        y=y_col,
        color=color_col if color_col else x_col,
        box=show_box,
        title=violin_title,
        points="all" if show_points else False,
        custom_data=custom_data_cols,
        template="plotly_white"
    )

    return fig


def create_bar_chart(
    df: pd.DataFrame,
    x_col: str,
    y_col: Optional[str] = None,
    color_col: Optional[str] = None,
    title: Optional[str] = None,
    orientation: str = "v",
    use_value_counts: bool = True
) -> go.Figure:
    """Create a bar chart.

    Args:
        df: DataFrame with data
        x_col: X-axis column (category)
        y_col: Y-axis column (values) - if None and use_value_counts=True, shows counts
        color_col: Optional column for color encoding
        title: Chart title
        orientation: 'v' for vertical, 'h' for horizontal
        use_value_counts: If True and y_col is None, show value counts

    Returns:
        Plotly Figure
    """
    if y_col is None and use_value_counts:
        # Create bar chart from value counts
        value_counts = df[x_col].value_counts().reset_index()
        value_counts.columns = [x_col, 'Count']

        fig = px.bar(
            value_counts,
            x=x_col,
            y='Count',
            title=title or f"Count by {x_col.replace('_', ' ')}",
            color=x_col if len(value_counts) <= 20 else None,
            template="plotly_white",
            orientation=orientation
        )
    else:
        fig = px.bar(
            df,
            x=x_col,
            y=y_col,
            color=color_col,
            title=title or f"{y_col} by {x_col}".replace('_', ' '),
            template="plotly_white",
            orientation=orientation
        )

    return fig


def create_line_plot(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: Optional[str] = None,
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None,
    title: Optional[str] = None,
    show_markers: bool = True
) -> go.Figure:
    """Create a line plot with structure viewer support.

    Args:
        df: DataFrame with data
        x_col: X-axis column
        y_col: Y-axis column
        color_col: Optional column for color/grouping
        smiles_col: SMILES column for structure viewer
        name_col: Name/ID column for display
        title: Chart title
        show_markers: Whether to show markers on the line

    Returns:
        Plotly Figure
    """
    plot_df = df.sort_values(x_col).copy()

    # Build customdata for structure viewer
    custom_data_cols = None
    if smiles_col and smiles_col in plot_df.columns:
        plot_df['_row_index'] = range(len(plot_df))
        if name_col and name_col in plot_df.columns:
            custom_data_cols = [smiles_col, name_col, '_row_index']
        else:
            custom_data_cols = [smiles_col, '_row_index']

    fig = px.line(
        plot_df,
        x=x_col,
        y=y_col,
        color=color_col,
        title=title or f"{y_col.replace('_', ' ')} vs {x_col.replace('_', ' ')}",
        markers=show_markers,
        custom_data=custom_data_cols,
        template="plotly_white"
    )

    return fig


def create_3d_scatter(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    z_col: str,
    color_col: Optional[str] = None,
    size_col: Optional[str] = None,
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None,
    title: Optional[str] = None,
    color_scale: str = "Viridis",
    marker_size: int = 5,
    opacity: float = 0.8
) -> go.Figure:
    """Create a 3D scatter plot with structure viewer support.

    Args:
        df: DataFrame with data
        x_col: X-axis column
        y_col: Y-axis column
        z_col: Z-axis column
        color_col: Optional column for color encoding
        size_col: Optional column for size encoding
        smiles_col: SMILES column for structure viewer
        name_col: Name/ID column for display
        title: Chart title
        color_scale: Color scale for continuous color
        marker_size: Base marker size
        opacity: Marker opacity

    Returns:
        Plotly Figure
    """
    # Prepare customdata for structure viewer
    if smiles_col:
        df, customdata_cols = prepare_chart_customdata(df, smiles_col, name_col)
    else:
        customdata_cols = None

    fig = px.scatter_3d(
        df,
        x=x_col,
        y=y_col,
        z=z_col,
        color=color_col,
        size=size_col,
        title=title or f"3D: {x_col} vs {y_col} vs {z_col}",
        color_continuous_scale=color_scale,
        opacity=opacity,
        custom_data=customdata_cols if customdata_cols else None
    )

    # Update marker size if not using size column
    if not size_col:
        fig.update_traces(marker=dict(size=marker_size))

    fig.update_layout(
        scene=dict(
            xaxis_title=x_col.replace('_', ' '),
            yaxis_title=y_col.replace('_', ' '),
            zaxis_title=z_col.replace('_', ' '),
        )
    )

    return fig


def create_correlation_heatmap(
    df: pd.DataFrame,
    columns: Optional[List[str]] = None,
    title: str = "Correlation Heatmap",
    color_scale: str = "RdBu_r"
) -> go.Figure:
    """Create a correlation heatmap.

    Args:
        df: DataFrame with data
        columns: Columns to include (None for all numeric)
        title: Chart title
        color_scale: Color scale

    Returns:
        Plotly Figure
    """
    if columns:
        corr_df = df[columns].select_dtypes(include=[np.number])
    else:
        corr_df = df.select_dtypes(include=[np.number])

    corr_matrix = corr_df.corr()

    fig = px.imshow(
        corr_matrix,
        title=title,
        color_continuous_scale=color_scale,
        aspect="auto",
        text_auto=".2f"
    )

    fig.update_layout(
        xaxis_title="",
        yaxis_title=""
    )

    return fig


def create_pair_plot(
    df: pd.DataFrame,
    columns: List[str],
    color_col: Optional[str] = None,
    title: str = "Pair Plot"
) -> go.Figure:
    """Create a scatter matrix (pair plot).

    Args:
        df: DataFrame with data
        columns: Columns to include
        color_col: Optional column for color encoding
        title: Chart title

    Returns:
        Plotly Figure
    """
    fig = px.scatter_matrix(
        df,
        dimensions=columns,
        color=color_col,
        title=title,
        template="plotly_white"
    )

    fig.update_traces(diagonal_visible=False)

    return fig


# =============================================================================
# Helper Functions
# =============================================================================

def _calculate_regression_stats(
    df: pd.DataFrame,
    x_col: str,
    y_col: str
) -> Optional[Dict[str, Any]]:
    """Calculate regression statistics for a scatter plot.

    Args:
        df: DataFrame with data
        x_col: X-axis column
        y_col: Y-axis column

    Returns:
        Dictionary with regression statistics or None if calculation fails
    """
    try:
        X = df[[x_col]].dropna()
        y_vals = df[y_col].dropna()
        common_idx = X.index.intersection(y_vals.index)
        X = X.loc[common_idx]
        y_vals = y_vals.loc[common_idx]

        if len(X) < 2:
            return None

        reg = LinearRegression()
        reg.fit(X, y_vals)
        y_pred = reg.predict(X)
        r2 = r2_score(y_vals, y_pred)
        slope = reg.coef_[0]
        intercept = reg.intercept_
        correlation = df[x_col].corr(df[y_col])

        return {
            'slope': slope,
            'intercept': intercept,
            'r2': r2,
            'correlation': correlation if not pd.isna(correlation) else None,
            'n': len(X)
        }

    except Exception as e:
        logger.warning(f"Failed to calculate regression stats: {e}")
        return None


def _add_regression_annotation(
    fig: go.Figure,
    stats: Dict[str, Any],
    x_col: str,
    y_col: str
) -> None:
    """Add regression statistics annotation to a figure.

    Args:
        fig: Plotly figure
        stats: Regression statistics dictionary
        x_col: X-axis column name
        y_col: Y-axis column name
    """
    slope = stats['slope']
    intercept = stats['intercept']
    r2 = stats['r2']
    correlation = stats['correlation']

    sign = "+" if intercept >= 0 else "-"
    y_name = y_col.replace('_', ' ')
    x_name = x_col.replace('_', ' ')
    equation = f"{y_name} = {slope:.3f} × {x_name} {sign} {abs(intercept):.3f}"

    # Handle NaN correlation (per CLAUDE.md: don't set to 0.0, use N/A instead)
    corr_text = f"{correlation:.3f}" if correlation is not None else "N/A"

    fig.add_annotation(
        x=0.02, y=0.98,
        xref="paper", yref="paper",
        text=f"<b>Regression Statistics:</b><br>Equation: {equation}<br>R² = {r2:.3f}<br>Correlation: {corr_text}",
        showarrow=False,
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="black",
        borderwidth=1,
        align="left",
        font=dict(size=11)
    )


# =============================================================================
# Utility Functions
# =============================================================================

def render_chart_with_viewer(
    fig: go.Figure,
    chart_id: str = "chart",
    x_col: Optional[str] = None,
    y_col: Optional[str] = None,
    z_col: Optional[str] = None,
    name_col: Optional[str] = None,
    show_hint: bool = True,
    width: str = 'stretch'
) -> None:
    """Render a Plotly chart with structure viewer integration.

    Args:
        fig: Plotly Figure
        chart_id: Unique identifier for the chart
        x_col: X-axis column name
        y_col: Y-axis column name
        z_col: Z-axis column name (for 3D)
        name_col: Name/ID column
        show_hint: Whether to show the viewer hint
        width: Chart width ('stretch' for full width, 'content' for auto)
    """
    if show_hint:
        render_structure_viewer_hint()

    st.plotly_chart(fig, width=width, key=f"{chart_id}_plotly")

    embed_structure_viewer(
        chart_id=chart_id,
        x_col=x_col,
        y_col=y_col,
        z_col=z_col,
        name_col=name_col
    )


def get_available_chart_types() -> Dict[str, str]:
    """Get dictionary of available chart types.

    Returns:
        Dictionary mapping chart type names to descriptions
    """
    return {
        "Scatter Plot": "Show relationship between two variables with optional trendline",
        "Histogram": "Show distribution of a single variable",
        "Box Plot": "Show distribution with quartiles",
        "Violin Plot": "Show distribution with density",
        "Bar Chart": "Compare values across categories",
        "Line Plot": "Show trends over continuous data",
        "3D Scatter": "Show relationship between three variables",
        "Correlation Heatmap": "Show correlations between all numeric columns",
        "Pair Plot": "Show all pairwise relationships",
    }


def render_chart_controls(
    df: pd.DataFrame,
    key_prefix: str = "chart"
) -> Dict[str, Any]:
    """Render chart configuration controls.

    Args:
        df: DataFrame with data
        key_prefix: Prefix for widget keys

    Returns:
        Dictionary of selected options
    """
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    all_cols = df.columns.tolist()

    config = {}

    col1, col2 = st.columns(2)

    with col1:
        config['chart_type'] = st.selectbox(
            "Chart Type",
            options=list(get_available_chart_types().keys()),
            key=f"{key_prefix}_type"
        )

        config['x_col'] = st.selectbox(
            "X Axis",
            options=numeric_cols,
            key=f"{key_prefix}_x"
        )

    with col2:
        config['y_col'] = st.selectbox(
            "Y Axis",
            options=numeric_cols,
            index=min(1, len(numeric_cols) - 1),
            key=f"{key_prefix}_y"
        )

        config['color_col'] = st.selectbox(
            "Color By",
            options=["None"] + all_cols,
            key=f"{key_prefix}_color"
        )

    if config['color_col'] == "None":
        config['color_col'] = None

    return config
