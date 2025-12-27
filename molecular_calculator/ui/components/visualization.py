"""Shared visualization component.

This module provides a comprehensive visualization component that can be used
by both batch processing and data visualization pages, avoiding code duplication.
"""

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from typing import Optional, Set, Dict, List, Any
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# Import theme utilities
try:
    from ..theme import apply_theme, apply_3d_theme, get_color_scale, COLOR_SCALES
    THEME_AVAILABLE = True
except ImportError:
    THEME_AVAILABLE = False

# Import ThreeDOLSRegression from the legacy molecular_calculator.py file
# Using importlib to avoid naming conflict with the molecular_calculator package
import importlib.util
import os as _os

_module_path = _os.path.join(
    _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))),
    'molecular_calculator.py'
)
_spec = importlib.util.spec_from_file_location("molecular_calculator_legacy", _module_path)
_legacy_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_legacy_module)
ThreeDOLSRegression = _legacy_module.ThreeDOLSRegression

# Import structure viewer component
try:
    from structure_viewer_component import (
        get_structure_viewer_component,
        get_structure_viewer_hint
    )
    STRUCTURE_VIEWER_AVAILABLE = True
except ImportError:
    STRUCTURE_VIEWER_AVAILABLE = False

# Import export utilities
try:
    from ..utils.export import get_download_data, EXPORT_FORMATS
    EXPORT_AVAILABLE = True
except ImportError:
    EXPORT_AVAILABLE = False


def _render_export_buttons(
    fig: go.Figure,
    chart_type: str,
    key_prefix: str,
) -> None:
    """
    Render export buttons for a Plotly chart.

    Args:
        fig: Plotly figure to export
        chart_type: Type of chart for filename
        key_prefix: Unique key prefix
    """
    if not EXPORT_AVAILABLE:
        return

    with st.expander("📥 Export Chart", expanded=False):
        export_cols = st.columns(4)

        # Clean chart type for filename
        base_filename = chart_type.lower().replace(' ', '_').replace('/', '_')

        formats = ['png', 'svg', 'pdf', 'html']

        for i, fmt in enumerate(formats):
            with export_cols[i]:
                try:
                    data, mime, ext = get_download_data(fig, fmt)
                    st.download_button(
                        label=f"{fmt.upper()}",
                        data=data,
                        file_name=f"{base_filename}{ext}",
                        mime=mime,
                        key=f"{key_prefix}_export_{fmt}",
                        use_container_width=True,
                    )
                except Exception as e:
                    st.button(
                        f"{fmt.upper()}",
                        disabled=True,
                        key=f"{key_prefix}_export_{fmt}_disabled",
                        help=f"Export failed: {str(e)}",
                        use_container_width=True,
                    )


def render_distribution_plots(
    df: pd.DataFrame,
    calculated_columns: Optional[Set[str]] = None,
    key_prefix: str = "dist"
) -> None:
    """Render quick distribution analysis plots.

    Args:
        df: Results DataFrame
        calculated_columns: Optional set of calculated column names to prioritize
        key_prefix: Unique key prefix for widgets
    """
    st.subheader("📈 Quick Distribution Analysis")

    # Priority properties for visualization
    priority_properties = ['QED', 'LogP', 'Molecular_Weight', 'TPSA']
    lei_properties = ['BEI', 'SEI', 'NSEI', 'NBEI', 'nBEI', 'mBEI', 'LEH', 'LEP']

    # Get available properties
    available_props = []

    if calculated_columns:
        # Use calculated columns
        for prop in priority_properties:
            if prop in df.columns and prop in calculated_columns:
                available_props.append(prop)
        # Add LEIs that were calculated
        for lei in lei_properties:
            if lei in df.columns and len(available_props) < 4:
                available_props.append(lei)
    else:
        # Use any available numeric columns
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        for prop in priority_properties:
            if prop in numeric_cols:
                available_props.append(prop)
        for lei in lei_properties:
            if lei in numeric_cols and len(available_props) < 4:
                available_props.append(lei)
        # Fill with other numeric columns if needed
        for col in numeric_cols:
            if col not in available_props and len(available_props) < 4:
                available_props.append(col)

    if not available_props:
        st.info("No numeric properties available for distribution plots.")
        return

    cols = st.columns(2)
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4']

    for i, prop in enumerate(available_props[:4]):
        with cols[i % 2]:
            fig = px.histogram(
                df,
                x=prop,
                title=f"{prop.replace('_', ' ')} Distribution",
                nbins=20,
                color_discrete_sequence=[colors[i % len(colors)]]
            )

            # Add mean line
            mean_val = df[prop].mean()
            if not pd.isna(mean_val):
                fig.add_vline(
                    x=mean_val,
                    line_dash="dash",
                    line_color="red",
                    annotation_text=f"Mean: {mean_val:.2f}",
                    annotation_position="top"
                )

            # Add reference lines for drug-like properties
            if prop == 'QED':
                fig.add_vline(x=0.5, line_dash="dot", line_color="orange",
                            annotation_text="QED > 0.5 (Drug-like)", annotation_position="bottom")
            elif prop == 'LogP':
                fig.add_vline(x=5, line_dash="dot", line_color="orange",
                            annotation_text="LogP < 5 (Lipinski)", annotation_position="bottom")
            elif prop == 'Molecular_Weight':
                fig.add_vline(x=500, line_dash="dot", line_color="orange",
                            annotation_text="MW < 500 (Lipinski)", annotation_position="bottom")
            elif prop == 'TPSA':
                fig.add_vline(x=140, line_dash="dot", line_color="orange",
                            annotation_text="TPSA < 140 (Drug-like)", annotation_position="bottom")

            fig.update_layout(
                xaxis_title=prop.replace('_', ' '),
                yaxis_title="Count",
                height=400,
                showlegend=False,
                plot_bgcolor='rgba(0,0,0,0)',
                paper_bgcolor='rgba(0,0,0,0)'
            )
            st.plotly_chart(fig, use_container_width=True, key=f"{key_prefix}_{prop}")


def render_interactive_visualization(
    df: pd.DataFrame,
    key_prefix: str = "viz",
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None
) -> None:
    """Render interactive data visualization section.

    This is a shared component used by both batch processing and data visualization pages.

    Args:
        df: DataFrame with data to visualize
        key_prefix: Unique key prefix for all widgets (must be unique per page)
        smiles_col: Optional SMILES column for structure viewer integration
        name_col: Optional name/ID column for display
    """
    # Get column types
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    all_cols = df.columns.tolist()
    categorical_cols = [col for col in all_cols
                       if df[col].dtype == 'object' or df[col].nunique() <= 20]

    if len(numeric_cols) < 1:
        st.warning("No numeric columns available for visualization.")
        return

    st.subheader("🎨 Interactive Data Visualization")

    # Chart type definitions
    chart_types = {
        'Scatter Plot': {'requires': ['x', 'y'], 'desc': 'Compare two numeric variables with optional trendline'},
        '3D OLS Regression': {'requires': ['x', 'y', 'z'], 'desc': 'Fit 3D plane: Z = b₀ + b₁·X + b₂·Y'},
        'Histogram': {'requires': ['x'], 'desc': 'Distribution of a single variable'},
        'Box Plot': {'requires': ['y'], 'desc': 'Statistical summary with quartiles'},
        'Violin Plot': {'requires': ['y'], 'desc': 'Distribution shape with density'},
        'Bar Chart': {'requires': ['cat'], 'desc': 'Count of categories'},
        'Line Plot': {'requires': ['x', 'y'], 'desc': 'Trends over continuous data'},
        'Heatmap': {'requires': [], 'desc': 'Correlation matrix of all numeric variables'}
    }

    # Chart type selection
    col1, col2 = st.columns([1, 2])
    with col1:
        selected_chart = st.selectbox(
            "📊 Chart Type:",
            options=list(chart_types.keys()),
            key=f"{key_prefix}_chart_type"
        )
    with col2:
        st.caption(f"ℹ️ {chart_types[selected_chart]['desc']}")

    requires = chart_types[selected_chart]['requires']

    # Initialize axis variables
    x_axis = None
    y_axis = None
    z_axis = None
    color_col = None
    group_by_col = None

    # Variable selection based on chart type
    # Add extra column for trendline toggle on scatter plots
    if selected_chart == 'Scatter Plot':
        num_cols = 4  # X, Y, Color, Trendline
    elif 'z' in requires:
        num_cols = 4
    else:
        num_cols = 3
    input_cols = st.columns(num_cols)

    with input_cols[0]:
        if 'x' in requires:
            x_axis = st.selectbox(
                "🔢 X-axis:",
                options=numeric_cols,
                key=f"{key_prefix}_x_axis"
            )
        elif 'cat' in requires:
            x_axis = st.selectbox(
                "📋 Category:",
                options=all_cols,
                key=f"{key_prefix}_cat_axis"
            )

    with input_cols[1]:
        if 'y' in requires:
            y_axis = st.selectbox(
                "🔢 Y-axis:",
                options=numeric_cols,
                index=min(1, len(numeric_cols) - 1) if len(numeric_cols) > 1 else 0,
                key=f"{key_prefix}_y_axis"
            )
        elif selected_chart == 'Histogram':
            # For histogram, x is the variable
            pass
        elif selected_chart in ['Box Plot', 'Violin Plot']:
            y_axis = st.selectbox(
                "🔢 Variable:",
                options=numeric_cols,
                key=f"{key_prefix}_box_var"
            )

    if 'z' in requires:
        with input_cols[2]:
            z_axis = st.selectbox(
                "🔢 Z-axis (Dependent):",
                options=numeric_cols,
                index=min(2, len(numeric_cols) - 1) if len(numeric_cols) > 2 else 0,
                key=f"{key_prefix}_z_axis"
            )
        col_idx = 3
    else:
        col_idx = 2

    # Color column selection
    if selected_chart not in ['Heatmap', '3D OLS Regression']:
        with input_cols[col_idx]:
            color_options = ['None'] + all_cols
            color_selection = st.selectbox(
                "🎨 Color by:",
                options=color_options,
                key=f"{key_prefix}_color"
            )
            color_col = None if color_selection == 'None' else color_selection

    # Trendline toggle for Scatter Plot (next to Color by)
    show_trendline = False
    if selected_chart == 'Scatter Plot':
        with input_cols[3]:
            st.markdown("<br>", unsafe_allow_html=True)  # Add spacing to align with dropdowns
            show_trendline = st.toggle(
                "📈 Trendline",
                value=False,
                key=f"{key_prefix}_trendline",
                help="Show OLS regression trendline with R² and equation"
            )

    # Advanced options based on chart type
    size_col = None
    shape_col = None
    size_max = 20
    fixed_size = 8
    nbins = 30
    show_points = False
    color_scale = 'Viridis'

    if selected_chart == 'Scatter Plot':
        with st.expander("⚙️ Advanced Scatter Options", expanded=False):
            adv_cols = st.columns(3)

            with adv_cols[0]:
                size_options = ['None (Fixed size)'] + numeric_cols
                size_selection = st.selectbox(
                    "📏 Size by:",
                    options=size_options,
                    key=f"{key_prefix}_size_col",
                    help="Map point size to a numeric column"
                )
                size_col = None if size_selection == 'None (Fixed size)' else size_selection

            with adv_cols[1]:
                if size_col:
                    size_max = st.slider(
                        "Max point size:",
                        min_value=5,
                        max_value=50,
                        value=20,
                        key=f"{key_prefix}_size_max"
                    )
                else:
                    fixed_size = st.slider(
                        "Point size:",
                        min_value=3,
                        max_value=20,
                        value=8,
                        key=f"{key_prefix}_fixed_size"
                    )

            with adv_cols[2]:
                shape_options = ['None'] + categorical_cols
                shape_selection = st.selectbox(
                    "🔷 Shape by:",
                    options=shape_options,
                    key=f"{key_prefix}_shape_col",
                    help="Map marker shape to a categorical column"
                )
                shape_col = None if shape_selection == 'None' else shape_selection

    elif selected_chart == 'Histogram':
        with st.expander("⚙️ Histogram Options", expanded=False):
            nbins = st.slider(
                "Number of bins:",
                min_value=5,
                max_value=100,
                value=30,
                key=f"{key_prefix}_nbins"
            )

    elif selected_chart in ['Box Plot', 'Violin Plot']:
        with st.expander("⚙️ Box/Violin Options", expanded=False):
            opt_cols = st.columns(2)
            with opt_cols[0]:
                group_options = ['None'] + categorical_cols
                group_selection = st.selectbox(
                    "📊 Group by:",
                    options=group_options,
                    key=f"{key_prefix}_group_by",
                    help="Group data by a categorical column"
                )
                group_by_col = None if group_selection == 'None' else group_selection

            with opt_cols[1]:
                show_points = st.checkbox(
                    "Show all data points",
                    value=False,
                    key=f"{key_prefix}_show_points"
                )

    # Color scale options - show for any chart that can use colors
    show_color_options = (
        selected_chart in ['Heatmap', '3D OLS Regression'] or
        (color_col and color_col in numeric_cols) or
        selected_chart == 'Scatter Plot'  # Always show for scatter plots
    )

    if show_color_options:
        with st.expander("🎨 Color Scale Options", expanded=False):
            # Use centralized color scales if available
            if THEME_AVAILABLE:
                color_scales = list(COLOR_SCALES.keys())
            else:
                color_scales = ['Viridis', 'Plasma', 'Inferno', 'Magma', 'Cividis',
                              'Blues', 'Reds', 'Greens', 'Purples', 'Oranges',
                              'RdBu', 'RdYlGn', 'Spectral', 'Turbo']
            color_scale = st.selectbox(
                "Color scale:",
                options=color_scales,
                key=f"{key_prefix}_color_scale"
            )
            reverse_scale = st.checkbox("Reverse scale", key=f"{key_prefix}_reverse_scale")
            if reverse_scale:
                color_scale = color_scale + '_r'

    # Chart types that support clicking on points
    point_based_charts = ['Scatter Plot', '3D OLS Regression', 'Line Plot', 'Box Plot', 'Violin Plot']

    # Show structure viewer hint if SMILES column is available
    can_show_structures = smiles_col and smiles_col in df.columns and STRUCTURE_VIEWER_AVAILABLE
    if can_show_structures and selected_chart in point_based_charts:
        st.markdown(get_structure_viewer_hint(), unsafe_allow_html=True)

    # Create and display chart
    fig = None
    chart_rendered = False

    try:
        if selected_chart == 'Scatter Plot' and x_axis and y_axis:
            fig = _create_scatter_plot(
                df, x_axis, y_axis, color_col, size_col, shape_col,
                size_max, fixed_size, show_trendline, color_scale, numeric_cols, key_prefix,
                smiles_col=smiles_col, name_col=name_col
            )

        elif selected_chart == '3D OLS Regression' and x_axis and y_axis and z_axis:
            _render_3d_ols_regression(
                df, x_axis, y_axis, z_axis, color_scale, key_prefix,
                smiles_col=smiles_col, name_col=name_col
            )
            chart_rendered = True
            # Embed structure viewer for 3D chart
            if can_show_structures:
                viewer_html = get_structure_viewer_component(
                    chart_id=f"{key_prefix}_3d",
                    x_col=x_axis,
                    y_col=y_axis,
                    z_col=z_axis,
                    name_col=name_col
                )
                components.html(viewer_html, height=0)
            return  # 3D OLS handles its own rendering

        elif selected_chart == 'Histogram' and x_axis:
            fig = _create_histogram(df, x_axis, color_col, nbins)

        elif selected_chart == 'Box Plot' and y_axis:
            fig = _create_box_plot(df, y_axis, group_by_col, color_col, show_points,
                                   smiles_col=smiles_col, name_col=name_col)

        elif selected_chart == 'Violin Plot' and y_axis:
            fig = _create_violin_plot(df, y_axis, group_by_col, color_col, show_points,
                                      smiles_col=smiles_col, name_col=name_col)

        elif selected_chart == 'Bar Chart' and x_axis:
            fig = _create_bar_chart(df, x_axis)

        elif selected_chart == 'Line Plot' and x_axis and y_axis:
            fig = _create_line_plot(df, x_axis, y_axis, color_col,
                                    smiles_col=smiles_col, name_col=name_col)

        elif selected_chart == 'Heatmap':
            fig = _create_heatmap(df, numeric_cols, color_scale)

        if fig:
            # Apply centralized theme if available
            if THEME_AVAILABLE:
                fig = apply_theme(fig)
            fig.update_layout(
                height=500,
                plot_bgcolor='rgba(0,0,0,0)',
                paper_bgcolor='rgba(0,0,0,0)'
            )
            st.plotly_chart(fig, use_container_width=True, key=f"{key_prefix}_main_chart")
            chart_rendered = True

            # Export buttons
            _render_export_buttons(fig, selected_chart, key_prefix)

            # Embed structure viewer for point-based charts
            if can_show_structures and selected_chart in point_based_charts:
                viewer_html = get_structure_viewer_component(
                    chart_id=f"{key_prefix}_chart",
                    x_col=x_axis,
                    y_col=y_axis,
                    name_col=name_col
                )
                components.html(viewer_html, height=0)

    except Exception as e:
        st.error(f"Error creating chart: {str(e)}")


def _create_scatter_plot(
    df: pd.DataFrame,
    x_axis: str,
    y_axis: str,
    color_col: Optional[str],
    size_col: Optional[str],
    shape_col: Optional[str],
    size_max: int,
    fixed_size: int,
    show_trendline: bool,
    color_scale: str,
    numeric_cols: List[str],
    key_prefix: str,
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None
) -> go.Figure:
    """Create a scatter plot with all options."""
    scatter_params = {
        'x': x_axis,
        'y': y_axis,
        'title': f"{y_axis.replace('_', ' ')} vs {x_axis.replace('_', ' ')}"
    }

    if color_col:
        scatter_params['color'] = color_col
        if color_col in numeric_cols:
            scatter_params['color_continuous_scale'] = color_scale

    if size_col:
        plot_df = df.dropna(subset=[size_col]).copy()
        scatter_params['size'] = size_col
        scatter_params['size_max'] = size_max
    else:
        plot_df = df.copy()

    if shape_col:
        scatter_params['symbol'] = shape_col

    if show_trendline:
        scatter_params['trendline'] = 'ols'

    # Add customdata for structure viewer if SMILES column is available
    if smiles_col and smiles_col in plot_df.columns:
        # Add row index for tracking
        plot_df['_row_index'] = range(len(plot_df))
        if name_col and name_col in plot_df.columns:
            scatter_params['custom_data'] = [smiles_col, name_col, '_row_index']
        else:
            scatter_params['custom_data'] = [smiles_col, '_row_index']

    fig = px.scatter(plot_df, **scatter_params)

    if not size_col:
        fig.update_traces(marker=dict(size=fixed_size))

    # Add regression statistics annotation
    if show_trendline:
        X = plot_df[[x_axis]].dropna()
        y_vals = plot_df[y_axis].dropna()
        common_idx = X.index.intersection(y_vals.index)
        X = X.loc[common_idx]
        y_vals = y_vals.loc[common_idx]

        if len(X) > 1:
            reg = LinearRegression()
            reg.fit(X, y_vals)
            y_pred = reg.predict(X)
            r2 = r2_score(y_vals, y_pred)
            slope = reg.coef_[0]
            intercept = reg.intercept_
            correlation = plot_df[x_axis].corr(plot_df[y_axis])

            sign = "+" if intercept >= 0 else "-"
            y_name = y_axis.replace('_', ' ')
            x_name = x_axis.replace('_', ' ')
            equation = f"{y_name} = {slope:.3f} × {x_name} {sign} {abs(intercept):.3f}"

            fig.add_annotation(
                x=0.02, y=0.98,
                xref="paper", yref="paper",
                text=f"<b>Regression Statistics:</b><br>" +
                     f"Equation: {equation}<br>" +
                     f"R² = {r2:.3f}<br>" +
                     f"Correlation: {correlation:.3f}",
                showarrow=False,
                bgcolor="rgba(255,255,255,0.9)",
                bordercolor="black",
                borderwidth=1,
                align="left",
                font=dict(size=11)
            )

    return fig


def _render_3d_ols_regression(
    df: pd.DataFrame,
    x_axis: str,
    y_axis: str,
    z_axis: str,
    color_scale: str,
    key_prefix: str,
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None
) -> None:
    """Render 3D OLS regression with fitted plane."""
    # Prepare columns to keep for plot data
    cols_needed = [x_axis, y_axis, z_axis]
    if smiles_col and smiles_col in df.columns:
        cols_needed.append(smiles_col)
    if name_col and name_col in df.columns:
        cols_needed.append(name_col)

    plot_data = df[cols_needed].dropna(subset=[x_axis, y_axis, z_axis])

    if len(plot_data) < 3:
        st.warning("Need at least 3 data points for 3D regression")
        return

    ols_model = ThreeDOLSRegression(
        x=plot_data[x_axis].values,
        y=plot_data[y_axis].values,
        z=plot_data[z_axis].values
    )

    stats = ols_model.get_statistics()
    equation = ols_model.get_equation_string(decimals=3)

    # Create 3D scatter plot
    fig = go.Figure()

    # Prepare customdata for structure viewer
    customdata = None
    if smiles_col and smiles_col in plot_data.columns:
        if name_col and name_col in plot_data.columns:
            customdata = list(zip(
                plot_data[smiles_col].values,
                plot_data[name_col].values,
                range(len(plot_data))
            ))
        else:
            customdata = list(zip(
                plot_data[smiles_col].values,
                range(len(plot_data))
            ))

    # Color by residuals
    residuals_abs = np.abs(ols_model.residuals)
    scatter_trace = go.Scatter3d(
        x=ols_model.x,
        y=ols_model.y,
        z=ols_model.z,
        mode='markers',
        marker=dict(
            size=5,
            color=residuals_abs,
            colorscale='Reds',
            showscale=True,
            colorbar=dict(title='|Residual|')
        ),
        name='Data Points',
        text=[f"{x_axis}: {x:.2f}<br>{y_axis}: {y:.2f}<br>{z_axis}: {z:.2f}<br>Residual: {r:.3f}"
              for x, y, z, r in zip(ols_model.x, ols_model.y, ols_model.z, ols_model.residuals)],
        hovertemplate='%{text}<extra></extra>'
    )

    if customdata:
        scatter_trace.customdata = customdata

    fig.add_trace(scatter_trace)

    # Add fitted plane
    x_min, x_max = float(np.min(ols_model.x)), float(np.max(ols_model.x))
    y_min, y_max = float(np.min(ols_model.y)), float(np.max(ols_model.y))
    x_pad = 0.1 * (x_max - x_min) if x_max > x_min else 1.0
    y_pad = 0.1 * (y_max - y_min) if y_max > y_min else 1.0

    corners_x = np.array([x_min - x_pad, x_max + x_pad, x_max + x_pad, x_min - x_pad])
    corners_y = np.array([y_min - y_pad, y_min - y_pad, y_max + y_pad, y_max + y_pad])
    corners_z = ols_model.predict(corners_x, corners_y)

    fig.add_trace(go.Mesh3d(
        x=corners_x,
        y=corners_y,
        z=corners_z,
        i=[0, 0],
        j=[1, 2],
        k=[2, 3],
        opacity=0.6,
        color='royalblue',
        name='Fitted Plane',
        hoverinfo='skip'
    ))

    fig.update_layout(
        title=f"3D OLS Regression: {z_axis} vs {x_axis} and {y_axis}",
        scene=dict(
            xaxis_title=x_axis.replace('_', ' '),
            yaxis_title=y_axis.replace('_', ' '),
            zaxis_title=z_axis.replace('_', ' '),
            aspectmode='data'
        ),
        height=700,
        showlegend=True
    )

    # Apply 3D theme if available
    if THEME_AVAILABLE:
        fig = apply_3d_theme(fig)

    st.plotly_chart(fig, use_container_width=True, config={'scrollZoom': True}, key=f"{key_prefix}_3d_chart")

    # Display regression statistics
    st.markdown("### 📊 3D OLS Regression Statistics")
    stat_cols = st.columns(5)
    with stat_cols[0]:
        st.metric("R²", f"{stats['R²']:.4f}")
    with stat_cols[1]:
        st.metric("RMSE", f"{stats['RMSE']:.4f}")
    with stat_cols[2]:
        st.metric("Intercept (b₀)", f"{stats['b0']:.4f}")
    with stat_cols[3]:
        st.metric("Coefficient b₁", f"{stats['b1']:.4f}")
    with stat_cols[4]:
        st.metric("Coefficient b₂", f"{stats['b2']:.4f}")

    st.markdown(f"**Fitted Plane Equation:** {equation}")
    st.caption("💡 Points are colored by their absolute residual. Drag to rotate the 3D plot.")


def _create_histogram(
    df: pd.DataFrame,
    x_axis: str,
    color_col: Optional[str],
    nbins: int
) -> go.Figure:
    """Create histogram with mean and median lines."""
    fig = px.histogram(
        df,
        x=x_axis,
        color=color_col,
        title=f"Distribution of {x_axis.replace('_', ' ')}",
        nbins=nbins
    )

    mean_val = df[x_axis].mean()
    median_val = df[x_axis].median()

    if not pd.isna(mean_val):
        fig.add_vline(x=mean_val, line_dash="dash", line_color="red",
                    annotation_text=f"Mean: {mean_val:.2f}")
    if not pd.isna(median_val):
        fig.add_vline(x=median_val, line_dash="dot", line_color="blue",
                    annotation_text=f"Median: {median_val:.2f}")

    return fig


def _create_box_plot(
    df: pd.DataFrame,
    y_axis: str,
    group_by_col: Optional[str],
    color_col: Optional[str],
    show_points: bool,
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None
) -> go.Figure:
    """Create box plot with options."""
    box_title = f"Box Plot of {y_axis.replace('_', ' ')}"
    if group_by_col:
        box_title += f" by {group_by_col.replace('_', ' ')}"

    plot_df = df.copy()

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
        x=group_by_col,
        y=y_axis,
        color=color_col if color_col else group_by_col,
        title=box_title,
        points="all" if show_points else False,
        custom_data=custom_data_cols
    )

    return fig


def _create_violin_plot(
    df: pd.DataFrame,
    y_axis: str,
    group_by_col: Optional[str],
    color_col: Optional[str],
    show_points: bool,
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None
) -> go.Figure:
    """Create violin plot with options."""
    violin_title = f"Violin Plot of {y_axis.replace('_', ' ')}"
    if group_by_col:
        violin_title += f" by {group_by_col.replace('_', ' ')}"

    plot_df = df.copy()

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
        x=group_by_col,
        y=y_axis,
        color=color_col if color_col else group_by_col,
        box=True,
        title=violin_title,
        points="all" if show_points else False,
        custom_data=custom_data_cols
    )

    return fig


def _create_bar_chart(df: pd.DataFrame, x_axis: str) -> go.Figure:
    """Create bar chart from value counts."""
    value_counts = df[x_axis].value_counts().reset_index()
    value_counts.columns = [x_axis, 'Count']

    fig = px.bar(
        value_counts,
        x=x_axis,
        y='Count',
        title=f"Count by {x_axis.replace('_', ' ')}",
        color=x_axis if len(value_counts) <= 20 else None
    )

    return fig


def _create_line_plot(
    df: pd.DataFrame,
    x_axis: str,
    y_axis: str,
    color_col: Optional[str],
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None
) -> go.Figure:
    """Create line plot."""
    plot_df = df.sort_values(x_axis).copy()

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
        x=x_axis,
        y=y_axis,
        color=color_col,
        title=f"{y_axis.replace('_', ' ')} vs {x_axis.replace('_', ' ')}",
        markers=True,
        custom_data=custom_data_cols
    )

    return fig


def _create_heatmap(
    df: pd.DataFrame,
    numeric_cols: List[str],
    color_scale: str
) -> go.Figure:
    """Create correlation heatmap."""
    corr_matrix = df[numeric_cols].corr()

    fig = px.imshow(
        corr_matrix,
        text_auto='.2f',
        title="Correlation Matrix",
        color_continuous_scale=color_scale,
        aspect='auto'
    )

    return fig
