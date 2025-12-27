# ==============================================================================
# Centralized Theme Configuration
# ==============================================================================
"""
Centralized theme settings for consistent styling across all visualizations.
Provides color palettes, chart themes, and styling utilities.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any

import plotly.graph_objects as go
import plotly.io as pio


# ==============================================================================
# Color Palettes
# ==============================================================================

@dataclass
class ColorPalette:
    """Color palette configuration."""
    primary: str = "#1f77b4"      # Blue
    secondary: str = "#ff7f0e"    # Orange
    success: str = "#2ca02c"      # Green
    warning: str = "#d62728"      # Red
    info: str = "#9467bd"         # Purple
    muted: str = "#7f7f7f"        # Gray

    # Categorical colors for data series
    categorical: List[str] = field(default_factory=lambda: [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    ])

    # Sequential color scales for continuous data
    sequential_blue: str = "Blues"
    sequential_green: str = "Greens"
    sequential_red: str = "Reds"
    sequential_purple: str = "Purples"

    # Diverging color scales
    diverging: str = "RdBu"
    diverging_alt: str = "RdYlGn"

    # Default continuous scale
    continuous: str = "Viridis"


# Available color scales for user selection
COLOR_SCALES: Dict[str, str] = {
    "Viridis": "Viridis",
    "Plasma": "Plasma",
    "Inferno": "Inferno",
    "Magma": "Magma",
    "Cividis": "Cividis",
    "Blues": "Blues",
    "Greens": "Greens",
    "Reds": "Reds",
    "Purples": "Purples",
    "Oranges": "Oranges",
    "YlOrRd": "YlOrRd",
    "YlGnBu": "YlGnBu",
    "RdBu": "RdBu",
    "RdYlGn": "RdYlGn",
    "Spectral": "Spectral",
    "Turbo": "Turbo",
}


# ==============================================================================
# Chart Theme Configuration
# ==============================================================================

@dataclass
class ChartTheme:
    """Theme configuration for Plotly charts."""

    # Font settings
    font_family: str = "Arial, sans-serif"
    font_size: int = 12
    title_font_size: int = 16
    axis_font_size: int = 11

    # Colors
    background_color: str = "white"
    plot_background_color: str = "white"
    grid_color: str = "#e5e5e5"
    text_color: str = "#333333"

    # Layout
    margin: Dict[str, int] = field(default_factory=lambda: {
        "l": 60, "r": 40, "t": 60, "b": 60
    })

    # Legend
    legend_orientation: str = "h"
    legend_x: float = 0.5
    legend_y: float = -0.15
    legend_xanchor: str = "center"

    # Hover
    hover_mode: str = "closest"

    # Color palette
    colors: ColorPalette = field(default_factory=ColorPalette)


# Default theme instance
DEFAULT_THEME = ChartTheme()


# ==============================================================================
# Theme Application Functions
# ==============================================================================

def apply_theme(
    fig: go.Figure,
    theme: Optional[ChartTheme] = None,
    title: Optional[str] = None,
) -> go.Figure:
    """
    Apply consistent theme to a Plotly figure.

    Args:
        fig: Plotly figure to style
        theme: Theme configuration (uses default if None)
        title: Optional title to set

    Returns:
        Styled figure
    """
    if theme is None:
        theme = DEFAULT_THEME

    update_dict = {
        "font": {
            "family": theme.font_family,
            "size": theme.font_size,
            "color": theme.text_color,
        },
        "paper_bgcolor": theme.background_color,
        "plot_bgcolor": theme.plot_background_color,
        "margin": theme.margin,
        "hovermode": theme.hover_mode,
        "legend": {
            "orientation": theme.legend_orientation,
            "x": theme.legend_x,
            "y": theme.legend_y,
            "xanchor": theme.legend_xanchor,
        },
        "xaxis": {
            "gridcolor": theme.grid_color,
            "tickfont": {"size": theme.axis_font_size},
        },
        "yaxis": {
            "gridcolor": theme.grid_color,
            "tickfont": {"size": theme.axis_font_size},
        },
    }

    if title:
        update_dict["title"] = {
            "text": title,
            "font": {"size": theme.title_font_size},
            "x": 0.5,
            "xanchor": "center",
        }

    fig.update_layout(**update_dict)
    return fig


def apply_3d_theme(
    fig: go.Figure,
    theme: Optional[ChartTheme] = None,
    title: Optional[str] = None,
) -> go.Figure:
    """
    Apply theme to 3D Plotly figures.

    Args:
        fig: Plotly 3D figure to style
        theme: Theme configuration
        title: Optional title

    Returns:
        Styled 3D figure
    """
    if theme is None:
        theme = DEFAULT_THEME

    # Base theme
    apply_theme(fig, theme, title)

    # 3D-specific settings
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                backgroundcolor=theme.plot_background_color,
                gridcolor=theme.grid_color,
                tickfont=dict(size=theme.axis_font_size),
            ),
            yaxis=dict(
                backgroundcolor=theme.plot_background_color,
                gridcolor=theme.grid_color,
                tickfont=dict(size=theme.axis_font_size),
            ),
            zaxis=dict(
                backgroundcolor=theme.plot_background_color,
                gridcolor=theme.grid_color,
                tickfont=dict(size=theme.axis_font_size),
            ),
        )
    )
    return fig


def get_color_scale(name: str) -> str:
    """
    Get color scale by name with fallback.

    Args:
        name: Color scale name

    Returns:
        Valid color scale name
    """
    return COLOR_SCALES.get(name, "Viridis")


def get_categorical_colors(n: int, theme: Optional[ChartTheme] = None) -> List[str]:
    """
    Get n categorical colors from the theme.

    Args:
        n: Number of colors needed
        theme: Theme configuration

    Returns:
        List of color hex codes
    """
    if theme is None:
        theme = DEFAULT_THEME

    colors = theme.colors.categorical
    # Cycle through colors if more needed
    return [colors[i % len(colors)] for i in range(n)]


# ==============================================================================
# Streamlit Custom CSS
# ==============================================================================

STREAMLIT_CUSTOM_CSS = """
<style>
    /* Metric cards styling */
    .stMetric {
        background-color: #f8f9fa;
        padding: 10px;
        border-radius: 5px;
    }

    /* Chart container styling */
    .chart-container {
        background-color: white;
        padding: 15px;
        border-radius: 8px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    }

    /* Export button styling */
    .export-button {
        margin-top: 10px;
    }

    /* Info box styling */
    .info-box {
        background-color: #e7f3ff;
        padding: 12px;
        border-radius: 6px;
        border-left: 4px solid #1f77b4;
        margin: 10px 0;
    }

    /* Warning box styling */
    .warning-box {
        background-color: #fff3e0;
        padding: 12px;
        border-radius: 6px;
        border-left: 4px solid #ff7f0e;
        margin: 10px 0;
    }

    /* Hide Streamlit footer */
    footer {visibility: hidden;}
</style>
"""


def inject_custom_css():
    """Inject custom CSS into Streamlit app."""
    import streamlit as st
    st.markdown(STREAMLIT_CUSTOM_CSS, unsafe_allow_html=True)


# ==============================================================================
# Chart Size Presets
# ==============================================================================

CHART_SIZES: Dict[str, Dict[str, int]] = {
    "small": {"width": 400, "height": 300},
    "medium": {"width": 600, "height": 450},
    "large": {"width": 800, "height": 600},
    "wide": {"width": 1000, "height": 400},
    "tall": {"width": 500, "height": 700},
    "square": {"width": 500, "height": 500},
}


def get_chart_size(preset: str) -> Dict[str, int]:
    """
    Get chart dimensions by preset name.

    Args:
        preset: Size preset name

    Returns:
        Dictionary with width and height
    """
    return CHART_SIZES.get(preset, CHART_SIZES["medium"])
