"""
Enhanced 3D Visualization for OLS Regression

This module provides advanced 3D plotting features with improved:
- Camera controls and rotation
- Visual effects (lighting, shadows, materials)
- Interactive features
- Animation capabilities
- Export options

Developed by: Yashwanth Reddy for ITR-UIC
Part of: Chemo-Informatics Toolkit
"""

import plotly.graph_objects as go
import plotly.express as px
import numpy as np
from typing import Optional, Dict, Any, List, Tuple
from molecular_calculator import ThreeDOLSRegression


class Enhanced3DPlot:
    """
    Enhanced 3D plotting for regression analysis with advanced visual features
    """

    # Professional color schemes
    COLOR_SCHEMES = {
        'Default': {
            'points': 'Reds',
            'plane': 'Blues',
            'background': 'white'
        },
        'Professional': {
            'points': 'Viridis',
            'plane': 'YlGnBu',
            'background': 'white'
        },
        'Dark Mode': {
            'points': 'Plasma',
            'plane': 'Twilight',
            'background': '#1e1e1e'
        },
        'High Contrast': {
            'points': 'Hot',
            'plane': 'Ice',
            'background': 'white'
        },
        'Publication': {
            'points': 'RdYlBu_r',
            'plane': 'Blues',
            'background': 'white'
        },
        'Colorblind Safe': {
            'points': 'Cividis',
            'plane': 'Purples',
            'background': 'white'
        }
    }

    # Camera preset positions
    CAMERA_PRESETS = {
        'Default': {'eye': {'x': 1.5, 'y': 1.5, 'z': 1.3}},
        'Top View': {'eye': {'x': 0, 'y': 0, 'z': 2.5}},
        'Side View X': {'eye': {'x': 2.5, 'y': 0, 'z': 0}},
        'Side View Y': {'eye': {'x': 0, 'y': 2.5, 'z': 0}},
        'Isometric': {'eye': {'x': 1.25, 'y': 1.25, 'z': 1.25}},
        'Bird\'s Eye': {'eye': {'x': 0.5, 'y': 0.5, 'z': 2.5}},
        'Close-up': {'eye': {'x': 1.0, 'y': 1.0, 'z': 1.0}},
        'Far View': {'eye': {'x': 2.5, 'y': 2.5, 'z': 2.5}}
    }

    def __init__(self, model: ThreeDOLSRegression, x_name: str, y_name: str, z_name: str):
        """
        Initialize enhanced 3D plot

        Args:
            model: Fitted ThreeDOLSRegression model
            x_name: Name of X variable
            y_name: Name of Y variable
            z_name: Name of Z variable
        """
        self.model = model
        self.x_name = x_name
        self.y_name = y_name
        self.z_name = z_name

    def create_enhanced_plot(
        self,
        color_scheme: str = 'Default',
        camera_preset: str = 'Default',
        show_residuals: bool = True,
        residual_lines: bool = False,
        plane_opacity: float = 0.7,
        point_size: int = 6,
        mesh_resolution: int = 30,
        show_grid: bool = True,
        show_axes_labels: bool = True,
        title: Optional[str] = None,
        height: int = 700,
        width: Optional[int] = None
    ) -> go.Figure:
        """
        Create enhanced 3D visualization with advanced features

        Args:
            color_scheme: Color scheme from COLOR_SCHEMES
            camera_preset: Camera position from CAMERA_PRESETS
            show_residuals: Color points by residuals
            residual_lines: Draw lines from points to plane
            plane_opacity: Opacity of fitted plane (0-1)
            point_size: Size of scatter points
            mesh_resolution: Resolution of plane mesh
            show_grid: Show grid lines
            show_axes_labels: Show axis labels
            title: Custom title
            height: Plot height in pixels
            width: Plot width in pixels

        Returns:
            Plotly Figure object
        """
        # Get color scheme
        colors = self.COLOR_SCHEMES.get(color_scheme, self.COLOR_SCHEMES['Default'])

        # Create figure
        fig = go.Figure()

        # Get plane mesh with higher resolution
        X_mesh, Y_mesh, Z_mesh = self.model.get_plane_mesh(num_points=mesh_resolution)

        # Add fitted plane with enhanced styling
        fig.add_trace(go.Surface(
            x=X_mesh,
            y=Y_mesh,
            z=Z_mesh,
            opacity=plane_opacity,
            colorscale=colors['plane'],
            showscale=False,
            name='Fitted Plane',
            hovertemplate='<b>Fitted Plane</b><br>%{x:.3f}, %{y:.3f}<br>Predicted: %{z:.3f}<extra></extra>',
            lighting=dict(
                ambient=0.6,
                diffuse=0.8,
                fresnel=0.2,
                specular=0.5,
                roughness=0.5
            ),
            lightposition=dict(
                x=1000,
                y=1000,
                z=2000
            ),
            contours={
                "z": {
                    "show": True,
                    "usecolormap": True,
                    "highlightcolor": "limegreen",
                    "project": {"z": True}
                }
            }
        ))

        # Prepare point colors
        if show_residuals:
            residuals_abs = np.abs(self.model.residuals)
            point_colors = residuals_abs
            colorbar_title = '|Residual|'
        else:
            point_colors = self.model.z
            colorbar_title = self.z_name

        # Add data points with enhanced styling
        hover_text = [
            f"<b>{self.x_name}:</b> {x:.3f}<br>" +
            f"<b>{self.y_name}:</b> {y:.3f}<br>" +
            f"<b>{self.z_name}:</b> {z:.3f}<br>" +
            f"<b>Predicted:</b> {pred:.3f}<br>" +
            f"<b>Residual:</b> {res:+.3f}"
            for x, y, z, pred, res in zip(
                self.model.x,
                self.model.y,
                self.model.z,
                self.model.predict(self.model.x, self.model.y),
                self.model.residuals
            )
        ]

        fig.add_trace(go.Scatter3d(
            x=self.model.x,
            y=self.model.y,
            z=self.model.z,
            mode='markers',
            marker=dict(
                size=point_size,
                color=point_colors,
                colorscale=colors['points'],
                showscale=True,
                colorbar=dict(
                    title=colorbar_title,
                    x=1.15,
                    len=0.7,
                    thickness=20,
                    tickfont=dict(size=10)
                ),
                line=dict(
                    color='rgba(0, 0, 0, 0.8)',
                    width=1.5
                ),
                opacity=1.0
            ),
            name='Data Points',
            text=hover_text,
            hovertemplate='%{text}<extra></extra>'
        ))

        # Add residual lines if requested
        if residual_lines:
            self._add_residual_lines(fig, colors)

        # Get camera settings
        camera = self.CAMERA_PRESETS.get(camera_preset, self.CAMERA_PRESETS['Default'])

        # Update layout with enhanced styling
        layout_config = dict(
            title=dict(
                text=title or f"3D OLS Regression: {self.z_name} vs {self.x_name} and {self.y_name}",
                font=dict(size=16, family='Arial, sans-serif', color='#333'),
                x=0.5,
                xanchor='center'
            ),
            scene=dict(
                xaxis=self._create_axis_config(self.x_name, show_grid, show_axes_labels),
                yaxis=self._create_axis_config(self.y_name, show_grid, show_axes_labels),
                zaxis=self._create_axis_config(self.z_name, show_grid, show_axes_labels),
                camera=camera,
                aspectmode='cube',
                bgcolor=colors['background']
            ),
            height=height,
            width=width,
            showlegend=True,
            legend=dict(
                x=0.02,
                y=0.98,
                bgcolor='rgba(255, 255, 255, 0.8)',
                bordercolor='rgba(0, 0, 0, 0.2)',
                borderwidth=1,
                font=dict(size=10)
            ),
            paper_bgcolor=colors['background'],
            plot_bgcolor=colors['background'],
            hovermode='closest',
            margin=dict(l=0, r=0, t=40, b=0)
        )

        fig.update_layout(**layout_config)

        # Add modebar buttons
        fig.update_layout(
            modebar=dict(
                bgcolor='rgba(255, 255, 255, 0.7)',
                color='#333',
                activecolor='#2196F3'
            )
        )

        return fig

    def _add_residual_lines(self, fig: go.Figure, colors: Dict):
        """Add lines from data points to fitted plane"""
        predicted = self.model.predict(self.model.x, self.model.y)

        # Create line segments
        for i in range(len(self.model.x)):
            # Determine color based on residual sign
            color = 'rgba(255, 0, 0, 0.3)' if self.model.residuals[i] > 0 else 'rgba(0, 0, 255, 0.3)'

            fig.add_trace(go.Scatter3d(
                x=[self.model.x[i], self.model.x[i]],
                y=[self.model.y[i], self.model.y[i]],
                z=[self.model.z[i], predicted[i]],
                mode='lines',
                line=dict(color=color, width=1),
                showlegend=False,
                hoverinfo='skip'
            ))

    def _create_axis_config(self, label: str, show_grid: bool, show_labels: bool) -> Dict:
        """Create axis configuration"""
        config = {
            'title': dict(
                text=label.replace('_', ' ') if show_labels else '',
                font=dict(size=12, family='Arial, sans-serif')
            ),
            'showgrid': show_grid,
            'gridcolor': 'rgba(128, 128, 128, 0.2)',
            'showline': True,
            'linecolor': 'rgba(128, 128, 128, 0.5)',
            'linewidth': 1,
            'zeroline': True,
            'zerolinecolor': 'rgba(128, 128, 128, 0.3)',
            'showbackground': True,
            'backgroundcolor': 'rgba(240, 240, 240, 0.1)',
            'showspikes': False
        }
        return config

    def create_animation_frames(
        self,
        num_frames: int = 36,
        color_scheme: str = 'Default',
        **kwargs
    ) -> List[go.Frame]:
        """
        Create animation frames for auto-rotation

        Args:
            num_frames: Number of frames for full 360° rotation
            color_scheme: Color scheme to use
            **kwargs: Additional arguments for create_enhanced_plot

        Returns:
            List of plotly frames
        """
        frames = []
        angles = np.linspace(0, 360, num_frames)

        for i, angle in enumerate(angles):
            # Calculate camera position for rotation around Z-axis
            rad = np.radians(angle)
            camera = {
                'eye': {
                    'x': 2.5 * np.cos(rad),
                    'y': 2.5 * np.sin(rad),
                    'z': 1.3
                }
            }

            # Create frame data (simplified - just update camera)
            frame = go.Frame(
                name=f'frame_{i}',
                layout=dict(
                    scene=dict(camera=camera)
                )
            )
            frames.append(frame)

        return frames

    def create_multiple_views(self, color_scheme: str = 'Default') -> Dict[str, go.Figure]:
        """
        Create multiple viewpoint figures

        Args:
            color_scheme: Color scheme to use

        Returns:
            Dictionary of figures with different camera angles
        """
        views = {}

        for preset_name, camera_config in self.CAMERA_PRESETS.items():
            fig = self.create_enhanced_plot(
                color_scheme=color_scheme,
                camera_preset=preset_name,
                title=f"{self.z_name} vs {self.x_name} & {self.y_name} ({preset_name})"
            )
            views[preset_name] = fig

        return views

    def export_to_html(
        self,
        filename: str,
        include_plotlyjs: str = 'cdn',
        auto_open: bool = False,
        **plot_kwargs
    ):
        """
        Export plot to standalone HTML file

        Args:
            filename: Output filename
            include_plotlyjs: How to include plotly.js ('cdn', True, False)
            auto_open: Open in browser after export
            **plot_kwargs: Arguments for create_enhanced_plot
        """
        fig = self.create_enhanced_plot(**plot_kwargs)
        fig.write_html(
            filename,
            include_plotlyjs=include_plotlyjs,
            auto_open=auto_open,
            config={
                'displayModeBar': True,
                'displaylogo': False,
                'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
                'toImageButtonOptions': {
                    'format': 'png',
                    'filename': 'regression_3d',
                    'height': 1200,
                    'width': 1600,
                    'scale': 2
                }
            }
        )

    def get_plot_config(self) -> Dict[str, Any]:
        """
        Get recommended plot configuration for Streamlit

        Returns:
            Configuration dictionary for st.plotly_chart
        """
        return {
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
            'scrollZoom': True,
            'toImageButtonOptions': {
                'format': 'png',
                'filename': f'3d_regression_{self.z_name}',
                'height': 1200,
                'width': 1600,
                'scale': 2
            }
        }


def create_comparison_plot(
    models: List[ThreeDOLSRegression],
    model_names: List[str],
    x_name: str,
    y_name: str,
    z_name: str,
    color_scheme: str = 'Default'
) -> go.Figure:
    """
    Create comparison plot with multiple fitted planes

    Args:
        models: List of fitted models
        model_names: Names for each model
        x_name: X variable name
        y_name: Y variable name
        z_name: Z variable name
        color_scheme: Color scheme to use

    Returns:
        Plotly figure with multiple planes
    """
    fig = go.Figure()

    # Color schemes for multiple models
    plane_colors = ['Blues', 'Greens', 'Oranges', 'Purples', 'Reds']

    for idx, (model, name) in enumerate(zip(models, model_names)):
        # Add plane
        X_mesh, Y_mesh, Z_mesh = model.get_plane_mesh(num_points=20)

        fig.add_trace(go.Surface(
            x=X_mesh,
            y=Y_mesh,
            z=Z_mesh,
            opacity=0.5,
            colorscale=plane_colors[idx % len(plane_colors)],
            showscale=False,
            name=f'{name} (R²={model.r_squared:.3f})',
            hovertemplate=f'<b>{name}</b><br>Predicted: %{{z:.3f}}<extra></extra>'
        ))

        # Add points for first model only
        if idx == 0:
            fig.add_trace(go.Scatter3d(
                x=model.x,
                y=model.y,
                z=model.z,
                mode='markers',
                marker=dict(size=4, color='black'),
                name='Data Points',
                showlegend=True
            ))

    # Update layout
    fig.update_layout(
        title=f"Model Comparison: {z_name} vs {x_name} and {y_name}",
        scene=dict(
            xaxis_title=x_name,
            yaxis_title=y_name,
            zaxis_title=z_name,
            camera=dict(eye=dict(x=1.5, y=1.5, z=1.3))
        ),
        height=700,
        showlegend=True
    )

    return fig
