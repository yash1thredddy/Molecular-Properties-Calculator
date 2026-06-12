"""Shared visualization component.

This module provides UI orchestration for data visualization,
using the chart creation functions from charts.py.

This module handles:
- Streamlit widget setup and user interaction
- Chart type selection and configuration
- Structure viewer integration
- Theme application and export buttons

Chart creation logic is delegated to charts.py to avoid duplication.
"""

import hashlib
import logging
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from typing import Optional, Set, List

# Import theme utilities
from molecular_calculator.ui.theme import (
    apply_theme,
    apply_3d_theme,
    COLOR_SCALES,
)

# Import ThreeDOLSRegression from the package models
from molecular_calculator.models import ThreeDOLSRegression

# Import chart creation functions from the canonical source
from molecular_calculator.ui.components.charts import (
    create_scatter_plot,
    create_histogram,
    create_box_plot,
    create_violin_plot,
    create_bar_chart,
    create_line_plot,
    create_correlation_heatmap,
)

# Import structure viewer component
try:
    from molecular_calculator.ui.components.structure_viewer import (
        get_structure_viewer_component,
        get_structure_viewer_hint
    )
    STRUCTURE_VIEWER_AVAILABLE = True
except ImportError:
    STRUCTURE_VIEWER_AVAILABLE = False

# Import export utilities
try:
    from molecular_calculator.utils.export import get_download_data, EXPORT_FORMATS
    EXPORT_AVAILABLE = True
except ImportError:
    EXPORT_AVAILABLE = False

logger = logging.getLogger(__name__)


def _render_export_buttons(
    fig: go.Figure,
    chart_type: str,
    key_prefix: str,
) -> None:
    """Render export buttons for a Plotly chart.

    Uses session state to cache generated exports, avoiding repeated
    kaleido/choreo calls on every page rerun.

    Args:
        fig: Plotly figure to export
        chart_type: Type of chart for filename
        key_prefix: Unique key prefix
    """
    if not EXPORT_AVAILABLE:
        return

    with st.expander("Export Chart", expanded=False):
        base_filename = chart_type.lower().replace(' ', '_').replace('/', '_')
        formats = ['png', 'svg', 'pdf', 'html']

        # Create figure fingerprint to detect when chart changes
        fig_json = fig.to_json()
        fig_hash = hashlib.md5(fig_json.encode()).hexdigest()[:8]

        # Cache key includes figure hash to invalidate when chart changes
        cache_key = f"{key_prefix}_export_cache"
        hash_key = f"{key_prefix}_export_hash"

        # Check if figure changed - invalidate cache if so
        if hash_key in st.session_state and st.session_state[hash_key] != fig_hash:
            st.session_state[cache_key] = {}

        st.session_state[hash_key] = fig_hash

        # Check if exports are already generated and cached
        if cache_key not in st.session_state:
            st.session_state[cache_key] = {}

        cached_exports = st.session_state[cache_key]

        # Generate exports button - only generate when user requests
        if st.button("🔄 Generate Export Files", key=f"{key_prefix}_generate_exports",
                     help="Click to generate downloadable files (only needed once per chart)"):
            with st.spinner("Generating export files..."):
                for fmt in formats:
                    try:
                        data, mime, ext = get_download_data(fig, fmt)
                        cached_exports[fmt] = {
                            'data': data,
                            'mime': mime,
                            'ext': ext,
                            'success': True
                        }
                    except Exception as e:
                        cached_exports[fmt] = {
                            'success': False,
                            'error': str(e)
                        }
                st.session_state[cache_key] = cached_exports
                st.rerun()

        # Show download buttons if exports are cached
        if cached_exports:
            export_cols = st.columns(4)
            for i, fmt in enumerate(formats):
                with export_cols[i]:
                    if fmt in cached_exports:
                        export_data = cached_exports[fmt]
                        if export_data.get('success'):
                            st.download_button(
                                label=f"📥 {fmt.upper()}",
                                data=export_data['data'],
                                file_name=f"{base_filename}{export_data['ext']}",
                                mime=export_data['mime'],
                                key=f"{key_prefix}_export_{fmt}",
                            )
                        else:
                            st.button(
                                f"{fmt.upper()}",
                                disabled=True,
                                key=f"{key_prefix}_export_{fmt}_disabled",
                                help=f"Export failed: {export_data.get('error', 'Unknown error')}",
                            )
        else:
            st.info("Click 'Generate Export Files' to create downloadable versions.")


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
    st.subheader("Quick Distribution Analysis")

    # Priority properties for visualization
    priority_properties = ['QED', 'LogP', 'Molecular_Weight', 'TPSA']
    lei_properties = ['BEI', 'SEI', 'NSEI', 'NBEI', 'nBEI', 'mBEI', 'LEH', 'LEP']

    # Get available properties
    available_props = []

    if calculated_columns:
        for prop in priority_properties:
            if prop in df.columns and prop in calculated_columns:
                available_props.append(prop)
        for lei in lei_properties:
            if lei in df.columns and len(available_props) < 4:
                available_props.append(lei)
    else:
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        for prop in priority_properties:
            if prop in numeric_cols:
                available_props.append(prop)
        for lei in lei_properties:
            if lei in numeric_cols and len(available_props) < 4:
                available_props.append(lei)
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
            # Use the consolidated histogram function
            fig = create_histogram(
                df,
                column=prop,
                bins=20,
                title=f"{prop.replace('_', ' ')} Distribution",
                show_mean=True,
                show_median=False  # We'll add custom reference lines instead
            )

            # Update with custom color
            fig.update_traces(marker_color=colors[i % len(colors)])

            # Add reference lines for drug-like properties
            # These thresholds are based on established drug-likeness rules:
            # - QED (Quantitative Estimate of Drug-likeness): Values > 0.5 indicate
            #   favorable drug-like properties. Reference: Bickerton et al. (2012)
            # - LogP < 5: Lipinski's Rule of Five - compounds with LogP > 5 have
            #   poor absorption. Reference: Lipinski et al. (1997)
            # - MW < 500 Da: Lipinski's Rule of Five - oral drugs typically have
            #   molecular weight < 500 Da for good absorption
            # - TPSA < 140 A^2: Veber's rules - compounds with TPSA > 140 have
            #   poor oral bioavailability. Reference: Veber et al. (2002)
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
                height=400,
                showlegend=False
            )
            st.plotly_chart(fig, width='stretch', key=f"{key_prefix}_{prop}")


def render_interactive_visualization(
    df: pd.DataFrame,
    key_prefix: str = "viz",
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None
) -> None:
    """Render interactive data visualization section.

    This is a shared component used by both batch processing and data visualization pages.
    It orchestrates the UI while delegating chart creation to charts.py.

    Args:
        df: DataFrame with data to visualize
        key_prefix: Unique key prefix for all widgets (must be unique per page)
        smiles_col: Optional SMILES column for structure viewer integration
        name_col: Optional name/ID column for display
    """
    # Formula builder: let users add computed columns before axis/encoding selectors.
    from molecular_calculator.ui.components.formula_builder import render_formula_builder
    df = render_formula_builder(df, page_key=key_prefix)

    # Get column types
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    all_cols = df.columns.tolist()
    categorical_cols = [col for col in all_cols
                       if df[col].dtype == 'object' or df[col].nunique() <= 20]

    if len(numeric_cols) < 1:
        st.warning("No numeric columns available for visualization.")
        return

    st.subheader("Interactive Data Visualization")

    # Chart type definitions
    chart_types = {
        'Scatter Plot': {'requires': ['x', 'y'], 'desc': 'Compare two numeric variables with optional trendline'},
        '3D OLS Regression': {'requires': ['x', 'y', 'z'], 'desc': 'Fit 3D plane: Z = b0 + b1*X + b2*Y'},
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
            "Chart Type:",
            options=list(chart_types.keys()),
            key=f"{key_prefix}_chart_type"
        )
    with col2:
        st.caption(f"Info: {chart_types[selected_chart]['desc']}")

    requires = chart_types[selected_chart]['requires']

    # Initialize axis variables
    x_axis = None
    y_axis = None
    z_axis = None
    color_col = None
    group_by_col = None

    # Variable selection based on chart type
    if selected_chart == 'Scatter Plot':
        num_cols = 4
    elif 'z' in requires:
        num_cols = 4
    else:
        num_cols = 3
    input_cols = st.columns(num_cols)

    with input_cols[0]:
        if 'x' in requires:
            x_axis = st.selectbox(
                "X-axis:",
                options=numeric_cols,
                key=f"{key_prefix}_x_axis"
            )
        elif 'cat' in requires:
            x_axis = st.selectbox(
                "Category:",
                options=all_cols,
                key=f"{key_prefix}_cat_axis"
            )

    with input_cols[1]:
        if 'y' in requires:
            y_axis = st.selectbox(
                "Y-axis:",
                options=numeric_cols,
                index=min(1, len(numeric_cols) - 1) if len(numeric_cols) > 1 else 0,
                key=f"{key_prefix}_y_axis"
            )
        elif selected_chart in ['Box Plot', 'Violin Plot']:
            y_axis = st.selectbox(
                "Variable:",
                options=numeric_cols,
                key=f"{key_prefix}_box_var"
            )

    if 'z' in requires:
        with input_cols[2]:
            z_axis = st.selectbox(
                "Z-axis (Dependent):",
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
                "Color by:",
                options=color_options,
                key=f"{key_prefix}_color"
            )
            color_col = None if color_selection == 'None' else color_selection

    # Trendline toggle for Scatter Plot
    show_trendline = False
    if selected_chart == 'Scatter Plot':
        with input_cols[3]:
            st.markdown("<br>", unsafe_allow_html=True)
            show_trendline = st.toggle(
                "Trendline",
                value=False,
                key=f"{key_prefix}_trendline",
                help="Show OLS regression trendline with R-squared and equation"
            )

    # Advanced options based on chart type
    size_col = None
    shape_col = None
    size_max = 20
    marker_size = 8
    nbins = 30
    show_points = False
    color_scale = 'Viridis'

    if selected_chart == 'Scatter Plot':
        with st.expander("Advanced Scatter Options", expanded=False):
            adv_cols = st.columns(3)

            with adv_cols[0]:
                size_options = ['None (Fixed size)'] + numeric_cols
                size_selection = st.selectbox(
                    "Size by:",
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
                    marker_size = st.slider(
                        "Point size:",
                        min_value=3,
                        max_value=20,
                        value=8,
                        key=f"{key_prefix}_marker_size"
                    )

            with adv_cols[2]:
                shape_options = ['None'] + categorical_cols
                shape_selection = st.selectbox(
                    "Shape by:",
                    options=shape_options,
                    key=f"{key_prefix}_shape_col",
                    help="Map marker shape to a categorical column"
                )
                shape_col = None if shape_selection == 'None' else shape_selection

    elif selected_chart == 'Histogram':
        with st.expander("Histogram Options", expanded=False):
            nbins = st.slider(
                "Number of bins:",
                min_value=5,
                max_value=100,
                value=30,
                key=f"{key_prefix}_nbins"
            )

    elif selected_chart in ['Box Plot', 'Violin Plot']:
        with st.expander("Box/Violin Options", expanded=False):
            opt_cols = st.columns(2)
            with opt_cols[0]:
                group_options = ['None'] + categorical_cols
                group_selection = st.selectbox(
                    "Group by:",
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

    # Color scale options
    show_color_options = (
        selected_chart in ['Heatmap', '3D OLS Regression'] or
        (color_col and color_col in numeric_cols) or
        selected_chart == 'Scatter Plot'
    )

    if show_color_options:
        with st.expander("Color Scale Options", expanded=False):
            color_scales = list(COLOR_SCALES.keys())
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

    # Create and display chart using consolidated charts.py functions
    fig = None

    try:
        if selected_chart == 'Scatter Plot' and x_axis and y_axis:
            fig, _ = create_scatter_plot(
                df=df,
                x_col=x_axis,
                y_col=y_axis,
                color_col=color_col,
                size_col=size_col,
                shape_col=shape_col,
                smiles_col=smiles_col,
                name_col=name_col,
                trendline=show_trendline,
                color_scale=color_scale,
                marker_size=marker_size,
                size_max=size_max,
                key_prefix=key_prefix
            )

        elif selected_chart == '3D OLS Regression' and x_axis and y_axis and z_axis:
            _render_3d_ols_regression(
                df, x_axis, y_axis, z_axis, color_scale, key_prefix,
                smiles_col=smiles_col, name_col=name_col
            )
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
            fig = create_histogram(
                df=df,
                column=x_axis,
                bins=nbins,
                color_col=color_col
            )

        elif selected_chart == 'Box Plot' and y_axis:
            fig = create_box_plot(
                df=df,
                y_col=y_axis,
                x_col=group_by_col,
                color_col=color_col,
                smiles_col=smiles_col,
                name_col=name_col,
                show_points=show_points
            )

        elif selected_chart == 'Violin Plot' and y_axis:
            fig = create_violin_plot(
                df=df,
                y_col=y_axis,
                x_col=group_by_col,
                color_col=color_col,
                smiles_col=smiles_col,
                name_col=name_col,
                show_points=show_points
            )

        elif selected_chart == 'Bar Chart' and x_axis:
            fig = create_bar_chart(
                df=df,
                x_col=x_axis
            )

        elif selected_chart == 'Line Plot' and x_axis and y_axis:
            fig = create_line_plot(
                df=df,
                x_col=x_axis,
                y_col=y_axis,
                color_col=color_col,
                smiles_col=smiles_col,
                name_col=name_col
            )

        elif selected_chart == 'Heatmap':
            fig = create_correlation_heatmap(
                df=df,
                columns=numeric_cols,
                color_scale=color_scale
            )

        if fig:
            # Apply centralized theme
            fig = apply_theme(fig)
            fig.update_layout(height=650)  # Increased height for better visibility
            st.plotly_chart(fig, width='stretch', key=f"{key_prefix}_main_chart")

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
        logger.error(f"Error creating chart: {e}", exc_info=True)
        st.error(f"Error creating chart: {str(e)}")


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
    """Render 3D OLS regression with fitted plane.

    This function is unique to visualization.py as it combines
    regression model fitting with 3D visualization in a specific way.
    """
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

    try:
        ols_model = ThreeDOLSRegression(
            x=plot_data[x_axis].values,
            y=plot_data[y_axis].values,
            z=plot_data[z_axis].values
        )
    except ValueError as e:
        st.error(f"Regression error: {e}")
        return

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

    # Apply 3D theme
    fig = apply_3d_theme(fig)

    st.plotly_chart(fig, width='stretch', config={'scrollZoom': True}, key=f"{key_prefix}_3d_chart")

    # Display regression statistics
    st.markdown("### 3D OLS Regression Statistics")
    stat_cols = st.columns(5)
    with stat_cols[0]:
        st.metric("R-squared", f"{stats['R²']:.4f}")
    with stat_cols[1]:
        st.metric("RMSE", f"{stats['RMSE']:.4f}")
    with stat_cols[2]:
        st.metric("Intercept (b0)", f"{stats['b0']:.4f}")
    with stat_cols[3]:
        st.metric("Coefficient b1", f"{stats['b1']:.4f}")
    with stat_cols[4]:
        st.metric("Coefficient b2", f"{stats['b2']:.4f}")

    st.markdown(f"**Fitted Plane Equation:** {equation}")
    st.caption("Points are colored by their absolute residual. Drag to rotate the 3D plot.")
