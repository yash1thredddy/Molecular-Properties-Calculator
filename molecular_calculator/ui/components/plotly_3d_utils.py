"""3D Plotly visualization utilities.

Shared helper functions for 3D Plotly visualizations,
extracted from regression_3d_page.py for reuse.
"""

import numpy as np
import plotly.graph_objects as go
from typing import Tuple, Optional


def add_axis_arrows(
    fig: go.Figure,
    x_range: Tuple[float, float],
    y_range: Tuple[float, float],
    z_range: Tuple[float, float],
    show_labels: bool = True
) -> None:
    """Add axis arrows to a Plotly 3D figure for better orientation.

    Uses thick Scatter3d lines for shafts and small Cone heads.
    Colors follow common convention: X=red, Y=green, Z=blue.

    Args:
        fig: Plotly 3D figure
        x_range: (min, max) tuple for x-axis
        y_range: (min, max) tuple for y-axis
        z_range: (min, max) tuple for z-axis
        show_labels: Whether to show X, Y, Z labels
    """
    x0, x1 = float(x_range[0]), float(x_range[1])
    y0, y1 = float(y_range[0]), float(y_range[1])
    z0, z1 = float(z_range[0]), float(z_range[1])

    # Compute lengths with robust fallback when range is zero
    xr = (x1 - x0) if x1 > x0 else 1.0
    yr = (y1 - y0) if y1 > y0 else 1.0
    zr = (z1 - z0) if z1 > z0 else 1.0

    # Arrow sizes (shaft + head fractions of range)
    frac = 0.18
    head_frac = 0.22

    # Start at the min corner to avoid occluding the cloud
    xs, ys, zs = x0, y0, z0

    # Endpoint for shafts (slightly before full length to place heads)
    x_end = xs + xr * (frac * (1 - head_frac))
    y_end = ys + yr * (frac * (1 - head_frac))
    z_end = zs + zr * (frac * (1 - head_frac))

    # Colors (Plotly defaults palette for contrast on dark bg)
    col_x = '#EF553B'  # red
    col_y = '#00CC96'  # green
    col_z = '#636EFA'  # blue

    # Shafts: X, Y, Z
    fig.add_trace(go.Scatter3d(
        x=[xs, x_end], y=[ys, ys], z=[zs, zs], mode='lines',
        line=dict(color=col_x, width=7), showlegend=False, hoverinfo='skip'
    ))
    fig.add_trace(go.Scatter3d(
        x=[xs, xs], y=[ys, y_end], z=[zs, zs], mode='lines',
        line=dict(color=col_y, width=7), showlegend=False, hoverinfo='skip'
    ))
    fig.add_trace(go.Scatter3d(
        x=[xs, xs], y=[ys, ys], z=[zs, z_end], mode='lines',
        line=dict(color=col_z, width=7), showlegend=False, hoverinfo='skip'
    ))

    # Heads: three small cones with constant color
    fig.add_trace(go.Cone(
        x=[x_end], y=[ys], z=[zs], u=[xr * frac * head_frac], v=[0], w=[0],
        anchor='tail', sizemode='absolute', sizeref=max(xr, yr, zr) * 0.06,
        showscale=False, colorscale=[[0, col_x], [1, col_x]],
        name='X axis', hoverinfo='skip'
    ))
    fig.add_trace(go.Cone(
        x=[xs], y=[y_end], z=[zs], u=[0], v=[yr * frac * head_frac], w=[0],
        anchor='tail', sizemode='absolute', sizeref=max(xr, yr, zr) * 0.06,
        showscale=False, colorscale=[[0, col_y], [1, col_y]],
        name='Y axis', hoverinfo='skip'
    ))
    fig.add_trace(go.Cone(
        x=[xs], y=[ys], z=[z_end], u=[0], v=[0], w=[zr * frac * head_frac],
        anchor='tail', sizemode='absolute', sizeref=max(xr, yr, zr) * 0.06,
        showscale=False, colorscale=[[0, col_z], [1, col_z]],
        name='Z axis', hoverinfo='skip'
    ))

    if show_labels:
        # Place simple text labels slightly beyond heads
        fig.add_trace(go.Scatter3d(
            x=[x_end + xr * 0.03], y=[ys], z=[zs], mode='text', text=['X'],
            textfont=dict(color=col_x, size=12), showlegend=False, hoverinfo='skip'
        ))
        fig.add_trace(go.Scatter3d(
            x=[xs], y=[y_end + yr * 0.03], z=[zs], mode='text', text=['Y'],
            textfont=dict(color=col_y, size=12), showlegend=False, hoverinfo='skip'
        ))
        fig.add_trace(go.Scatter3d(
            x=[xs], y=[ys], z=[z_end + zr * 0.03], mode='text', text=['Z'],
            textfont=dict(color=col_z, size=12), showlegend=False, hoverinfo='skip'
        ))


def add_origin_crosshair(
    fig: go.Figure,
    x_range: Tuple[float, float],
    y_range: Tuple[float, float],
    z_range: Tuple[float, float],
    opacity: float = 0.3
) -> None:
    """Add faint origin crosshair lines for orientation.

    Args:
        fig: Plotly 3D figure
        x_range: (min, max) tuple for x-axis
        y_range: (min, max) tuple for y-axis
        z_range: (min, max) tuple for z-axis
        opacity: Line opacity
    """
    x0, x1 = float(x_range[0]), float(x_range[1])
    y0, y1 = float(y_range[0]), float(y_range[1])
    z0, z1 = float(z_range[0]), float(z_range[1])

    # Choose anchor values: use 0 if within range, else mid-point
    cx = 0.0 if (x0 <= 0.0 <= x1) else (x0 + x1) / 2.0
    cy = 0.0 if (y0 <= 0.0 <= y1) else (y0 + y1) / 2.0
    cz = 0.0 if (z0 <= 0.0 <= z1) else (z0 + z1) / 2.0

    col = 'rgba(127,140,141,1.0)'  # gray
    width = 2

    # X-axis line at y=cy, z=cz
    fig.add_trace(go.Scatter3d(
        x=[x0, x1], y=[cy, cy], z=[cz, cz], mode='lines',
        line=dict(color=col, width=width), opacity=opacity, showlegend=False
    ))
    # Y-axis line at x=cx, z=cz
    fig.add_trace(go.Scatter3d(
        x=[cx, cx], y=[y0, y1], z=[cz, cz], mode='lines',
        line=dict(color=col, width=width), opacity=opacity, showlegend=False
    ))
    # Z-axis line at x=cx, y=cy
    fig.add_trace(go.Scatter3d(
        x=[cx, cx], y=[cy, cy], z=[z0, z1], mode='lines',
        line=dict(color=col, width=width), opacity=opacity, showlegend=False
    ))


def add_residual_vectors(
    fig: go.Figure,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    z_predicted: np.ndarray,
    color: str = '#F1C40F',
    width: int = 2,
    opacity: float = 0.7,
    max_vectors: int = 300
) -> None:
    """Draw residual vectors from actual (x,y,z) to predicted on plane (x,y,z_predicted).

    Args:
        fig: Plotly 3D figure
        x: X coordinates
        y: Y coordinates
        z: Actual Z values
        z_predicted: Predicted Z values on the plane
        color: Vector color
        width: Line width
        opacity: Line opacity
        max_vectors: Maximum number of vectors to draw (for performance)
    """
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()
    z = np.asarray(z).ravel()
    z_predicted = np.asarray(z_predicted).ravel()
    n = x.shape[0]

    if n == 0:
        return

    # Subsample for performance
    if n > max_vectors:
        idx = np.linspace(0, n - 1, max_vectors, dtype=int)
        x, y, z, z_predicted = x[idx], y[idx], z[idx], z_predicted[idx]

    # Build segment arrays (NaN-separated)
    xs = np.empty(x.size * 3)
    xs[0::3] = x
    xs[1::3] = x
    xs[2::3] = np.nan

    ys = np.empty(y.size * 3)
    ys[0::3] = y
    ys[1::3] = y
    ys[2::3] = np.nan

    zs = np.empty(z.size * 3)
    zs[0::3] = z
    zs[1::3] = z_predicted
    zs[2::3] = np.nan

    fig.add_trace(go.Scatter3d(
        x=xs, y=ys, z=zs, mode='lines',
        line=dict(color=color, width=width), opacity=opacity,
        name='Residual vectors', showlegend=True
    ))


def add_predicted_markers(
    fig: go.Figure,
    x: np.ndarray,
    y: np.ndarray,
    z_predicted: np.ndarray,
    color: str = '#19D3F3',
    size: int = 3,
    opacity: float = 0.8,
    max_points: int = 1000
) -> None:
    """Add predicted markers on the fitted plane at (x, y, z_predicted).

    Args:
        fig: Plotly 3D figure
        x: X coordinates
        y: Y coordinates
        z_predicted: Predicted Z values on the plane
        color: Marker color
        size: Marker size
        opacity: Marker opacity
        max_points: Maximum number of points to draw (for performance)
    """
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()
    z_predicted = np.asarray(z_predicted).ravel()
    n = x.shape[0]

    if n == 0:
        return

    if n > max_points:
        idx = np.linspace(0, n - 1, max_points, dtype=int)
        x, y, z_predicted = x[idx], y[idx], z_predicted[idx]

    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z_predicted, mode='markers',
        marker=dict(size=size, color=color, opacity=opacity),
        name='Predicted (on plane)', showlegend=True
    ))


def create_fitted_plane_mesh(
    x: np.ndarray,
    y: np.ndarray,
    predict_func,
    padding_fraction: float = 0.1,
    color: str = 'royalblue',
    opacity: float = 0.6,
    name: str = 'Fitted Plane'
) -> go.Mesh3d:
    """Create a mesh for the fitted plane.

    Args:
        x: X coordinates of data
        y: Y coordinates of data
        predict_func: Function that takes (x, y) and returns z
        padding_fraction: How much to extend the plane beyond data
        color: Plane color
        opacity: Plane opacity
        name: Trace name

    Returns:
        Plotly Mesh3d trace
    """
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()

    x_min, x_max = float(np.min(x)), float(np.max(x))
    y_min, y_max = float(np.min(y)), float(np.max(y))

    x_pad = padding_fraction * (x_max - x_min) if x_max > x_min else 1.0
    y_pad = padding_fraction * (y_max - y_min) if y_max > y_min else 1.0

    corners_x = np.array([x_min - x_pad, x_max + x_pad, x_max + x_pad, x_min - x_pad])
    corners_y = np.array([y_min - y_pad, y_min - y_pad, y_max + y_pad, y_max + y_pad])
    corners_z = predict_func(corners_x, corners_y)

    return go.Mesh3d(
        x=corners_x,
        y=corners_y,
        z=corners_z,
        i=[0, 0],
        j=[1, 2],
        k=[2, 3],
        opacity=opacity,
        color=color,
        name=name,
        hoverinfo='skip'
    )
