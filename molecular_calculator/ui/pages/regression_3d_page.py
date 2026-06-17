"""3D OLS Regression Analysis Page.

Full-featured 3D OLS regression page extracted from the legacy
molecular_properties_app.py with all original features preserved.

Features:
- Full Statistical Summary (statsmodels style)
- Coefficient Table with p-values and confidence intervals
- Key Model Statistics (R^2, Adjusted R^2, RMSE, F-statistic)
- Fitted Plane Equation with Interpretation
- Statistical Significance of Coefficients
- 3D Visualization with options (axis arrows, crosshair, residual vectors, etc.)
- Residual Analysis (histogram, predicted vs actual)
- Diagnostic Tests (Durbin-Watson, Jarque-Bera, Condition Number)
- Export (TXT report, CSV stats, XLSX stats)
- Auto-suggest best predictor pairs
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from io import BytesIO
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Import regression module
from molecular_calculator.models.regression_3d import perform_3d_regression, RegressionSummary, suggest_best_3d_pairs

# Import 3D visualization helpers from shared module
from molecular_calculator.ui.components.plotly_3d_utils import (
    add_axis_arrows,
    add_origin_crosshair,
    add_residual_vectors,
    add_predicted_markers,
)

# Robust CSV reader (handles non-UTF-8 exports from Excel/Windows tools)
from molecular_calculator.ui.components.file_uploader import read_csv_robust
from molecular_calculator.utils.validators import FileValidator, DataFrameValidator


# ============================================================================
# Visualization rendering functions
# ============================================================================

def _render_threejs_visualization(model, summary, x_var: str, y_var: str, z_var: str) -> None:
    """Render Three.js enhanced 3D visualization."""
    import streamlit.components.v1 as components

    try:
        from molecular_calculator.ui.components.threejs_regression import get_threejs_regression_component

        # Get predicted values
        z_predicted = model.predict(model.x, model.y)

        # Generate Three.js component
        html_component = get_threejs_regression_component(
            x=model.x,
            y=model.y,
            z=model.z,
            z_predicted=z_predicted,
            residuals=model.residuals,
            b0=model.b0,
            b1=model.b1,
            b2=model.b2,
            r_squared=summary.r_squared,
            x_label=x_var,
            y_label=y_var,
            z_label=z_var,
            chart_id="reg3d_threejs",
            height=700
        )

        # Render the component
        components.html(html_component, height=750, scrolling=False)

        st.caption("Three.js visualization with interactive controls. Drag to rotate, scroll to zoom, right-click to pan.")

    except ImportError as e:
        st.warning(f"Three.js component not available: {e}")
        st.info("Falling back to Plotly visualization...")
        _render_plotly_visualization(model, summary, x_var, y_var, z_var)
    except Exception as e:
        st.error(f"Error rendering Three.js visualization: {e}")
        st.info("Falling back to Plotly visualization...")
        _render_plotly_visualization(model, summary, x_var, y_var, z_var)


def _render_plotly_visualization(model, summary, x_var: str, y_var: str, z_var: str) -> None:
    """Render Plotly 3D visualization (original implementation)."""
    # Get plane mesh
    X_mesh, Y_mesh, Z_mesh = model.get_plane_mesh(num_points=25)

    # Create 3D scatter plot
    fig = go.Figure()

    # Color by residuals
    residuals_abs = np.abs(model.residuals)

    # Add data points
    fig.add_trace(go.Scatter3d(
        x=model.x,
        y=model.y,
        z=model.z,
        mode='markers',
        marker=dict(
            size=7,
            color=residuals_abs,
            colorscale='Reds',
            showscale=True,
            colorbar=dict(title='|Residual|', x=1.2, xpad=10)
        ),
        name='Data Points',
        text=[f"{x_var}: {x:.3f}<br>{y_var}: {y:.3f}<br>{z_var}: {z:.3f}<br>Residual: {r:.3f}"
              for x, y, z, r in zip(model.x, model.y, model.z, model.residuals)],
        hovertemplate='%{text}<extra></extra>'
    ))

    # Add fitted plane
    fig.add_trace(go.Surface(
        x=X_mesh,
        y=Y_mesh,
        z=Z_mesh,
        opacity=0.7,
        colorscale='Blues',
        showscale=False,
        name='Fitted Plane',
        hovertemplate='Fitted Plane<br>Predicted Value<extra></extra>'
    ))

    # View controls
    ortho = st.checkbox("Orthographic projection", value=False, key="reg3d_ortho_plotly")
    camera_dict = dict(
        eye=dict(x=1.5, y=1.5, z=1.3),
        up=dict(x=0, y=0, z=1),
        projection=dict(type='orthographic' if ortho else 'perspective')
    )

    # Update layout
    fig.update_layout(
        title=f"3D OLS Regression: {z_var} vs {x_var} and {y_var}",
        scene=dict(
            xaxis_title=x_var,
            yaxis_title=y_var,
            zaxis_title=z_var,
            camera=camera_dict,
            aspectmode='data'
        ),
        height=700,
        showlegend=True,
        legend=dict(x=0.01, y=0.99, bgcolor='rgba(0,0,0,0)', borderwidth=0),
        margin=dict(r=120),
        scene_dragmode='orbit',
        uirevision='reg-3d-ols'
    )

    # Compute padded data ranges for visualization elements
    x_min_d, x_max_d = float(np.min(model.x)), float(np.max(model.x))
    y_min_d, y_max_d = float(np.min(model.y)), float(np.max(model.y))
    z_min_d, z_max_d = float(np.min(model.z)), float(np.max(model.z))
    x_pad = 0.1 * (x_max_d - x_min_d) if x_max_d > x_min_d else 1.0
    y_pad = 0.1 * (y_max_d - y_min_d) if y_max_d > y_min_d else 1.0
    x0, x1 = x_min_d - x_pad, x_max_d + x_pad
    y0, y1 = y_min_d - y_pad, y_max_d + y_pad
    z_plane_min, z_plane_max = float(np.min(Z_mesh)), float(np.max(Z_mesh))
    z0, z1 = min(z_min_d, z_plane_min), max(z_max_d, z_plane_max)

    # Optional visualization elements
    show_axis_arrows_reg = st.checkbox("Show 3D axis arrows", value=True, key="reg_axis_arrows_plotly")
    if show_axis_arrows_reg:
        add_axis_arrows(fig, (x0, x1), (y0, y1), (z0, z1), show_labels=True)

    show_crosshair_reg = st.checkbox("Show origin crosshair", value=True, key="reg_crosshair_plotly")
    if show_crosshair_reg:
        add_origin_crosshair(fig, (x0, x1), (y0, y1), (z0, z1), opacity=0.25)

    show_vectors_reg = st.checkbox("Show residual vectors", value=False, key="reg_residual_vectors_plotly")
    show_predicted_reg = st.checkbox("Show predicted markers", value=False, key="reg_predicted_markers_plotly")
    if show_vectors_reg or show_predicted_reg:
        zhat_pts = model.predict(model.x, model.y)
        if show_vectors_reg:
            add_residual_vectors(fig, model.x, model.y, model.z, zhat_pts)
        if show_predicted_reg:
            add_predicted_markers(fig, model.x, model.y, zhat_pts)

    st.plotly_chart(fig, width='stretch')
    st.caption("Points are colored by their absolute residual (distance from fitted plane). Use mouse to rotate the 3D plot.")


def _render_matplotlib_visualization(model, summary, x_var: str, y_var: str, z_var: str) -> None:
    """Render matplotlib static 3D visualization (scipy-lectures style).

    Creates a clean, publication-ready 3D plot similar to:
    https://scipy-lectures.org/packages/statistics/auto_examples/plot_regression_3d.html
    """
    # Create figure with 3D projection
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')

    # Get data
    x = model.x
    y = model.y
    z = model.z

    # Create mesh for fitted plane
    x_min, x_max = float(np.min(x)), float(np.max(x))
    y_min, y_max = float(np.min(y)), float(np.max(y))
    x_pad = 0.1 * (x_max - x_min) if x_max > x_min else 0.5
    y_pad = 0.1 * (y_max - y_min) if y_max > y_min else 0.5

    xx, yy = np.meshgrid(
        np.linspace(x_min - x_pad, x_max + x_pad, 30),
        np.linspace(y_min - y_pad, y_max + y_pad, 30)
    )

    # Calculate Z values on the fitted plane: z = b0 + b1*x + b2*y
    zz = model.b0 + model.b1 * xx + model.b2 * yy

    # Plot fitted plane with color gradient based on Z values
    surf = ax.plot_surface(
        xx, yy, zz,
        cmap='coolwarm',
        alpha=0.6,
        linewidth=0,
        antialiased=True,
        rstride=1,
        cstride=1
    )

    # Plot data points
    scatter = ax.scatter(
        x, y, z,
        c=z,  # Color by Z value to match plane colormap
        cmap='coolwarm',
        s=60,
        edgecolors='black',
        linewidths=0.5,
        alpha=0.9,
        depthshade=True
    )

    # Set labels with proper formatting
    ax.set_xlabel(f'\n{x_var}', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'\n{y_var}', fontsize=12, fontweight='bold')
    ax.set_zlabel(f'\n{z_var}', fontsize=12, fontweight='bold')

    # Set title
    ax.set_title(
        f'3D OLS Regression: {z_var} = f({x_var}, {y_var})\n'
        f'R² = {summary.r_squared:.4f}',
        fontsize=14,
        fontweight='bold',
        pad=20
    )

    # Improve viewing angle (similar to scipy example)
    ax.view_init(elev=20, azim=45)

    # Add equation annotation
    equation = f'{z_var} = {model.b0:.3f} + {model.b1:.3f}×{x_var} + {model.b2:.3f}×{y_var}'
    fig.text(0.5, 0.02, equation, ha='center', fontsize=11,
             style='italic', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Style adjustments
    ax.xaxis.pane.fill = True
    ax.yaxis.pane.fill = True
    ax.zaxis.pane.fill = True
    ax.xaxis.pane.set_facecolor((0.95, 0.95, 0.95, 0.8))
    ax.yaxis.pane.set_facecolor((0.90, 0.90, 0.90, 0.8))
    ax.zaxis.pane.set_facecolor((0.95, 0.95, 0.95, 0.8))

    # Grid styling
    ax.xaxis._axinfo['grid']['color'] = (0.7, 0.7, 0.7, 0.5)
    ax.yaxis._axinfo['grid']['color'] = (0.7, 0.7, 0.7, 0.5)
    ax.zaxis._axinfo['grid']['color'] = (0.7, 0.7, 0.7, 0.5)

    plt.tight_layout()

    # Display in Streamlit
    st.pyplot(fig)
    plt.close(fig)

    st.caption("Static matplotlib visualization. Clean, publication-ready 3D plot with fitted regression plane.")

    # Add download button for the figure
    buf = BytesIO()
    fig2 = plt.figure(figsize=(12, 9))
    ax2 = fig2.add_subplot(111, projection='3d')

    # Recreate the plot for saving
    surf2 = ax2.plot_surface(xx, yy, zz, cmap='coolwarm', alpha=0.6, linewidth=0, antialiased=True)
    scatter2 = ax2.scatter(x, y, z, c=z, cmap='coolwarm', s=60, edgecolors='black', linewidths=0.5, alpha=0.9)
    ax2.set_xlabel(f'\n{x_var}', fontsize=12, fontweight='bold')
    ax2.set_ylabel(f'\n{y_var}', fontsize=12, fontweight='bold')
    ax2.set_zlabel(f'\n{z_var}', fontsize=12, fontweight='bold')
    ax2.set_title(f'3D OLS Regression: {z_var} = f({x_var}, {y_var})\nR² = {summary.r_squared:.4f}',
                  fontsize=14, fontweight='bold', pad=20)
    ax2.view_init(elev=20, azim=45)
    fig2.text(0.5, 0.02, equation, ha='center', fontsize=11, style='italic',
              bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax2.xaxis.pane.set_facecolor((0.95, 0.95, 0.95, 0.8))
    ax2.yaxis.pane.set_facecolor((0.90, 0.90, 0.90, 0.8))
    ax2.zaxis.pane.set_facecolor((0.95, 0.95, 0.95, 0.8))
    plt.tight_layout()

    fig2.savefig(buf, format='png', dpi=150, bbox_inches='tight', facecolor='white')
    buf.seek(0)
    plt.close(fig2)

    st.download_button(
        label="📥 Download Plot (PNG)",
        data=buf,
        file_name=f"3d_regression_{x_var}_{y_var}_{z_var}.png",
        mime="image/png",
        key="download_matplotlib_plot"
    )


# ============================================================================
# Main page renderer
# ============================================================================

def render_3d_regression_page():
    """Render the full 3D OLS Regression Analysis page."""
    st.header("3D OLS Regression Analysis")
    st.markdown("Perform 3D Ordinary Least Squares regression: **Z = b_0 + b_1*X + b_2*Y**")
    st.markdown("Perfect for SAR (Structure-Activity Relationship) analysis and multi-variate modeling")

    render_3d_regression_help()

    # File upload
    reg_file = st.file_uploader(
        "Upload CSV file with your data", type=['csv'],
        accept_multiple_files=False,  # explicit: app processes one file at a time
        key="reg_upload",
    )

    if reg_file is not None:
        fres = FileValidator.validate_upload(reg_file)
        if not fres.is_valid:
            st.error(f"❌ {fres.errors[0] if fres.errors else 'Invalid file'}")
            st.stop()
        try:
            reg_df = read_csv_robust(reg_file)
        except Exception as e:
            st.error(f"Error reading file: {str(e)}")
            st.stop()
        dres = DataFrameValidator.validate(reg_df)
        if not dres.is_valid:
            st.error(f"❌ {dres.errors[0]}")
            st.stop()

        # Data preprocessing - convert numeric-like strings to numbers
        for col in reg_df.columns:
            if reg_df[col].dtype == 'object':
                numeric_converted = pd.to_numeric(reg_df[col], errors='coerce')
                non_null_original = reg_df[col].notna().sum()
                non_null_converted = numeric_converted.notna().sum()
                if non_null_converted > 0 and (non_null_converted / non_null_original) >= 0.5:
                    reg_df[col] = numeric_converted

        # Clear cached regression results if file changed
        if st.session_state.get('reg3d_file') != reg_file.name:
            st.session_state['reg3d_file'] = reg_file.name
            st.session_state.pop('reg3d_model', None)
            st.session_state.pop('reg3d_summary', None)
            st.session_state.pop('reg3d_vars', None)

        st.subheader("Data Preview (first 50 rows)")
        st.info(f"Loaded {len(reg_df)} rows and {len(reg_df.columns)} columns from {reg_file.name}")
        st.dataframe(reg_df.head(50), width='stretch')

        # Get numeric columns
        numeric_cols_reg = reg_df.select_dtypes(include=[np.number]).columns.tolist()

        if len(numeric_cols_reg) >= 3:
            st.markdown("---")
            st.subheader("Variable Selection")
            st.markdown("Select three numeric variables for regression: **Z = b_0 + b_1*X + b_2*Y**")

            col1, col2, col3 = st.columns(3)

            with col1:
                x_var = st.selectbox(
                    "Independent Variable 1 (X):",
                    options=numeric_cols_reg,
                    index=0,
                    key="reg_x_var",
                    help="First predictor variable"
                )

            with col2:
                available_y = [col for col in numeric_cols_reg if col != x_var]
                y_var = st.selectbox(
                    "Independent Variable 2 (Y):",
                    options=available_y,
                    index=0 if available_y else 0,
                    key="reg_y_var",
                    help="Second predictor variable"
                )

            with col3:
                available_z = [col for col in numeric_cols_reg if col != x_var and col != y_var]
                z_var = st.selectbox(
                    "Dependent Variable (Z):",
                    options=available_z,
                    index=0 if available_z else 0,
                    key="reg_z_var",
                    help="Response variable to be predicted"
                )

            # Add example suggestion for SAR data
            st.info(f"**Selected model:** {z_var} = b_0 + b_1*{x_var} + b_2*{y_var}")

            # Auto-suggest best (X, Y) pair for chosen Z
            with st.expander("Suggest best predictors (beta)", expanded=False):
                st.caption("Finds the best two predictors for the selected dependent variable using adjusted R^2 and BIC. Provides a ranked list as proof.")

                col_opt1, col_opt2, col_opt3 = st.columns([1, 1, 1])
                with col_opt1:
                    opt_inter = st.checkbox("Include interaction X*Y", value=False, help="Adds cross term for scoring only. 3D plot remains linear.")
                with col_opt2:
                    opt_quad = st.checkbox("Include quadratics X^2,Y^2", value=False, help="Adds squared terms for scoring only. 3D plot remains linear.")
                with col_opt3:
                    kfolds = st.slider("CV folds", min_value=3, max_value=10, value=5, step=1, help="K-fold CV for RMSE proof metric")

                if st.button("Suggest best (X, Y) pair", key="btn_suggest_xy"):
                    proof_df, top = suggest_best_3d_pairs(
                        reg_df, z_var,
                        include_interaction=opt_inter,
                        include_quadratic=opt_quad,
                        cv_folds=kfolds,
                    )
                    if top is None or proof_df.empty:
                        st.warning("No valid predictor pairs found. Ensure there are at least two numeric columns besides the dependent variable with enough data.")
                    else:
                        st.success(f"Suggested: X = {top['x']}, Y = {top['y']}  (Adj R^2 = {top['adj_r2']:.3f}, BIC = {top['bic']:.2f}, RMSE_CV = {top.get('rmse_cv', float('nan')):.3f}, Model = {top.get('model_spec', 'linear')})")
                        st.dataframe(proof_df, width='stretch')

                        # Let user apply the suggestion to the selectboxes
                        if st.button("Use this suggestion", key="btn_apply_suggestion"):
                            st.session_state["reg_x_var"] = top['x']
                            st.session_state["reg_y_var"] = top['y']
                            st.rerun()

            # Clear cached results if variable selection changed
            if st.session_state.get('reg3d_vars') and st.session_state['reg3d_vars'] != (x_var, y_var, z_var):
                st.session_state.pop('reg3d_model', None)
                st.session_state.pop('reg3d_summary', None)
                st.session_state.pop('reg3d_vars', None)

            # Perform regression button
            run_reg = st.button("Perform 3D OLS Regression", type="primary")
            if run_reg:
                try:
                    with st.spinner("Fitting 3D OLS regression model..."):
                        # Perform regression
                        model, summary = perform_3d_regression(reg_df, x_var, y_var, z_var)

                    # Cache results for persistence across widget changes
                    st.session_state['reg3d_model'] = model
                    st.session_state['reg3d_summary'] = summary
                    st.session_state['reg3d_vars'] = (x_var, y_var, z_var)

                    st.success("Regression analysis complete!")

                except Exception as e:
                    st.error(f"Error performing regression: {str(e)}")
                    import traceback
                    with st.expander("Debug Information"):
                        st.code(traceback.format_exc())

            # Display results if available (from button click or cache)
            if st.session_state.get('reg3d_model') is not None and \
               st.session_state.get('reg3d_vars') == (x_var, y_var, z_var):
                try:
                    model = st.session_state['reg3d_model']
                    summary = st.session_state['reg3d_summary']

                    _display_regression_results(model, summary, x_var, y_var, z_var, reg_file.name)

                except Exception as e:
                    st.error(f"Error displaying regression results: {str(e)}")
                    import traceback
                    with st.expander("Debug Information"):
                        st.code(traceback.format_exc())

        else:
            st.warning(f"Need at least 3 numeric columns for 3D regression. Found {len(numeric_cols_reg)} numeric columns.")
            st.info("Please upload a CSV file with at least 3 numeric columns (2 independent variables and 1 dependent variable).")

    else:
        st.info("Upload a CSV file to get started with 3D regression analysis")
        _display_help_info()


def _display_regression_results(model, summary, x_var: str, y_var: str, z_var: str, filename: str):
    """Display the full regression results with all statistics and visualizations."""
    # Display statsmodels-style summary
    st.markdown("---")
    st.subheader("OLS Regression Results")

    # Text summary
    with st.expander("Full Statistical Summary (statsmodels style)", expanded=True):
        st.code(summary.get_summary_text(), language=None)

    # Coefficient table
    st.markdown("### Coefficient Table")
    coef_df = summary.get_summary_dataframe()
    st.dataframe(coef_df, width='stretch')

    # Key metrics in columns
    st.markdown("### Key Model Statistics")
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("R^2 (Goodness of Fit)", f"{summary.r_squared:.4f}")
        quality = "Excellent" if summary.r_squared >= 0.9 else "Good" if summary.r_squared >= 0.7 else "Moderate" if summary.r_squared >= 0.5 else "Poor"
        st.caption(f"Quality: {quality}")
    with col2:
        st.metric("Adjusted R^2", f"{summary.adj_r_squared:.4f}")
        st.caption("Adjusted for # of predictors")
    with col3:
        st.metric("RMSE", f"{summary.model.rmse:.4f}")
        st.caption("Root Mean Squared Error")
    with col4:
        st.metric("F-statistic", f"{summary.f_statistic:.2f}")
        st.caption(f"p-value: {summary.f_pvalue:.2e}")

    # Equation display
    st.markdown("### Fitted Plane Equation")
    equation = model.get_equation_string(decimals=3)
    equation_display = equation.replace('Z', z_var).replace('X', x_var).replace('Y', y_var)
    st.code(equation_display, language=None)

    # Interpretation
    st.markdown("### Interpretation")
    st.markdown(f"""
    **Model Equation:** `{equation_display}`

    **What this means:**
    - **Intercept (b_0 = {model.b0:.3f})**: Baseline value of {z_var} when both {x_var} and {y_var} are zero
    - **{x_var} coefficient (b_1 = {model.b1:.3f})**: For each 1-unit increase in {x_var}, {z_var} changes by {model.b1:.3f} (holding {y_var} constant)
    - **{y_var} coefficient (b_2 = {model.b2:.3f})**: For each 1-unit increase in {y_var}, {z_var} changes by {model.b2:.3f} (holding {x_var} constant)
    - **R^2 = {summary.r_squared:.3f}**: This model explains {summary.r_squared*100:.1f}% of the variance in {z_var}
    """)

    # Coefficient significance
    st.markdown("### Statistical Significance of Coefficients")
    sig_data = []
    for i, (var_name, coef, se, t, p) in enumerate([
        ('Intercept', model.b0, summary.se_b0, summary.t_b0, summary.p_b0),
        (x_var, model.b1, summary.se_b1, summary.t_b1, summary.p_b1),
        (y_var, model.b2, summary.se_b2, summary.t_b2, summary.p_b2)
    ]):
        significance = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        sig_label = "Highly significant" if p < 0.001 else "Significant" if p < 0.01 else "Moderately significant" if p < 0.05 else "Not significant"
        sig_data.append({
            'Variable': var_name,
            'Coefficient': f"{coef:.4f}",
            'p-value': f"{p:.4e}",
            'Significance': significance,
            'Interpretation': sig_label
        })

    sig_df = pd.DataFrame(sig_data)
    st.dataframe(sig_df, width='stretch')
    st.caption("Significance codes: *** p<0.001, ** p<0.01, * p<0.05, ns: not significant")

    # 3D Visualization
    st.markdown("---")
    st.markdown("### 3D Visualization")

    # Visualization engine selector - Three.js is default (best quality)
    viz_engine = st.radio(
        "Visualization Engine:",
        ["Three.js (WebGL)", "Plotly (Interactive)", "Matplotlib (Static)"],
        horizontal=True,
        key="reg3d_viz_engine",
        help="Three.js: smooth WebGL rendering (recommended) | Plotly: interactive with controls | Matplotlib: publication-quality static plot"
    )

    if viz_engine == "Three.js (WebGL)":
        # Use Three.js visualization (default - best quality)
        _render_threejs_visualization(model, summary, x_var, y_var, z_var)
    elif viz_engine == "Plotly (Interactive)":
        # Use Plotly visualization
        _render_plotly_visualization(model, summary, x_var, y_var, z_var)
    else:
        # Use matplotlib static visualization
        _render_matplotlib_visualization(model, summary, x_var, y_var, z_var)

    # Residual Analysis
    st.markdown("---")
    st.markdown("### Residual Analysis")

    col1, col2 = st.columns(2)

    with col1:
        # Residual distribution
        fig_res_hist = px.histogram(
            x=model.residuals,
            nbins=20,
            title="Residual Distribution",
            labels={'x': 'Residual'},
            color_discrete_sequence=['#4ECDC4']
        )
        fig_res_hist.add_vline(x=0, line_dash="dash", line_color="red", annotation_text="Zero")
        fig_res_hist.update_layout(height=400, showlegend=False)
        st.plotly_chart(fig_res_hist, width='stretch')
        st.caption("Residuals should be normally distributed around zero")

    with col2:
        # Predicted vs Actual
        predicted = model.predict(model.x, model.y)
        fig_pred = px.scatter(
            x=predicted,
            y=model.z,
            title="Predicted vs Actual Values",
            labels={'x': f'Predicted {z_var}', 'y': f'Actual {z_var}'},
            color_discrete_sequence=['#FF6B6B']
        )
        # Add perfect prediction line
        min_val = min(predicted.min(), model.z.min())
        max_val = max(predicted.max(), model.z.max())
        fig_pred.add_trace(go.Scatter(
            x=[min_val, max_val],
            y=[min_val, max_val],
            mode='lines',
            line=dict(color='gray', dash='dash'),
            name='Perfect Prediction',
            showlegend=True
        ))
        fig_pred.update_layout(height=400)
        st.plotly_chart(fig_pred, width='stretch')
        st.caption("Points should lie close to the diagonal line for good fit")

    # Diagnostic tests summary
    st.markdown("### Diagnostic Tests")
    diag_col1, diag_col2, diag_col3 = st.columns(3)

    with diag_col1:
        st.metric("Durbin-Watson", f"{summary.durbin_watson:.3f}")
        dw_interpretation = "No autocorrelation" if 1.5 < summary.durbin_watson < 2.5 else "Possible autocorrelation"
        st.caption(f"{dw_interpretation} (ideal: ~2.0)")

    with diag_col2:
        st.metric("Jarque-Bera (JB)", f"{summary.jb_statistic:.3f}")
        jb_interpretation = "Normal residuals" if summary.jb_pvalue > 0.05 else "Non-normal residuals"
        st.caption(f"{jb_interpretation} (p={summary.jb_pvalue:.3f})")

    with diag_col3:
        st.metric("Condition Number", f"{summary.condition_number:.1f}")
        cn_interpretation = "Low multicollinearity" if summary.condition_number < 30 else "Moderate" if summary.condition_number < 100 else "High multicollinearity"
        st.caption(cn_interpretation)

    # Download results
    st.markdown("---")
    st.markdown("### Export Results")

    # Create summary report
    report_lines = [
        "=" * 80,
        "3D OLS REGRESSION ANALYSIS REPORT",
        "=" * 80,
        f"\nFile: {filename}",
        f"Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"\nModel: {z_var} = b0 + b1*{x_var} + b2*{y_var}",
        "\n" + summary.get_summary_text(),
        "\n" + "=" * 80,
        "END OF REPORT",
        "=" * 80
    ]

    # Append fitted plane equation and coefficient significance before summary
    _eq = model.get_equation_string(decimals=4)
    _eq_disp = _eq.replace('Z', z_var).replace('X', x_var).replace('Y', y_var)
    _sig_header = f"{'Variable':>12s} {'Coef':>12s} {'StdErr':>12s} {'t':>8s} {'P>|t|':>10s} {'Sig':>6s} {'[0.025':>10s} {'0.975]':>10s}"
    _sig_rows = []
    for _name, _coef, _se, _t, _p, _ci in [
        ('Intercept', model.b0, summary.se_b0, summary.t_b0, summary.p_b0, summary.ci_b0),
        (x_var, model.b1, summary.se_b1, summary.t_b1, summary.p_b1, summary.ci_b1),
        (y_var, model.b2, summary.se_b2, summary.t_b2, summary.p_b2, summary.ci_b2),
    ]:
        _sig = '***' if _p < 0.001 else '**' if _p < 0.01 else '*' if _p < 0.05 else 'ns'
        _sig_rows.append(f"{_name:>12s} {(_coef):>12.4f} {(_se):>12.4f} {(_t):>8.3f} {(_p):>10.3e} {(_sig):>6s} {(_ci[0]):>10.3f} {(_ci[1]):>10.3f}")

    _section_lines = [
        "\nFITTED PLANE EQUATION",
        "-" * 80,
        _eq_disp,
        "\nSTATISTICAL SIGNIFICANCE OF COEFFICIENTS",
        "-" * 80,
        _sig_header,
    ] + _sig_rows

    # Insert the new section before the long stats summary (at index 6)
    for _line in reversed(_section_lines):
        report_lines.insert(6, _line)

    report_text = "\n".join(report_lines)

    col1, col2, col3 = st.columns(3)
    with col1:
        st.download_button(
            label="Download Report (TXT)",
            data=report_text,
            file_name=f"3D_OLS_Regression_Report_{z_var}.txt",
            mime="text/plain"
        )

    # Export coefficients and stats
    stats_dict = summary.get_statistics_dict()
    stats_df = pd.DataFrame([stats_dict])
    with col2:
        st.download_button(
            label="Download Stats (CSV)",
            # UTF-8 BOM (utf-8-sig) so the file opens correctly in Excel on Windows.
            data=stats_df.to_csv(index=False).encode("utf-8-sig"),
            file_name=f"3D_OLS_Statistics_{z_var}.csv",
            mime="text/csv"
        )
    with col3:
        xlsx_buffer = BytesIO()
        stats_df.to_excel(xlsx_buffer, index=False, engine='openpyxl')
        xlsx_buffer.seek(0)
        st.download_button(
            label="Download Stats (XLSX)",
            data=xlsx_buffer.getvalue(),
            file_name=f"3D_OLS_Statistics_{z_var}.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )


def _display_help_info():
    """Display help information when no file is uploaded."""
    st.markdown("""
    ### What is 3D OLS Regression?

    3D Ordinary Least Squares regression fits a **plane** to your data in 3-dimensional space:

    **Z = b_0 + b_1*X + b_2*Y**

    Where:
    - **Z** is the dependent variable (what you want to predict)
    - **X** and **Y** are independent variables (predictors)
    - **b_0, b_1, b_2** are coefficients calculated to minimize prediction error

    ### Use Cases:

    **SAR Analysis (Structure-Activity Relationship)**:
    - Predict **biological activity (pKi)** from **LogP** and **TPSA**
    - Model **IC50** values from molecular descriptors

    **General Applications**:
    - Multi-variate modeling with two predictors
    - Understanding combined effects of two variables
    - Scientific data analysis requiring 3D relationships

    ### Output Includes:

    - Complete statistical summary (like statsmodels)
    - Coefficient table with p-values and confidence intervals
    - R^2, Adjusted R^2, RMSE, F-statistic
    - Interactive 3D visualization
    - Residual analysis plots
    - Diagnostic tests (Durbin-Watson, Jarque-Bera)
    - Downloadable report and statistics

    ### Example Data Structure:

    Your CSV should have at least 3 numeric columns:

    | LogP | TPSA | pKi  |
    |------|------|------|
    | 2.3  | 45.2 | 7.1  |
    | 3.1  | 62.8 | 6.8  |
    | 1.9  | 38.5 | 7.5  |
    | ...  | ...  | ...  |

    Upload your data above to start!
    """)


def render_3d_regression_help():
    """Render the help documentation for 3D regression."""
    with st.expander("3D Regression Analysis - Help & Documentation"):
        st.markdown("""
        ### 3D Ordinary Least Squares (OLS) Regression

        **Mathematical Model:**
        ```
        Z = b_0 + b_1*X + b_2*Y + e
        ```

        Where:
        - **Z** = Dependent variable (response, what you want to predict)
        - **X, Y** = Independent variables (predictors, features)
        - **b_0** = Intercept (baseline value when X=0, Y=0)
        - **b_1** = Coefficient for X (effect of X on Z, holding Y constant)
        - **b_2** = Coefficient for Y (effect of Y on Z, holding X constant)
        - **e** = Error term (residual)

        ---

        ### Statistical Output Explained

        #### **Model Quality Metrics:**
        - **R^2 (R-squared)**: Proportion of variance explained (0-1, higher is better)
          - 0.9+ = Excellent fit
          - 0.7-0.9 = Good fit
          - 0.5-0.7 = Moderate fit
          - <0.5 = Poor fit

        - **Adjusted R^2**: R^2 adjusted for number of predictors (penalizes overfitting)

        - **RMSE**: Root Mean Squared Error (average prediction error in original units)

        - **F-statistic**: Tests if the model is significantly better than a null model

        #### **Coefficient Interpretation:**
        - **b_0 (Intercept)**: Expected Z when X=0 and Y=0
        - **b_1 (X coefficient)**: Change in Z for each 1-unit increase in X (holding Y constant)
        - **b_2 (Y coefficient)**: Change in Z for each 1-unit increase in Y (holding X constant)

        #### **Diagnostic Tests:**
        - **Durbin-Watson**: Tests for autocorrelation in residuals (ideal: ~2.0)
        - **Jarque-Bera**: Tests if residuals are normally distributed
        - **Condition Number**: Measures multicollinearity (<30 is good)

        ---

        ### Interpretation Guide

        1. **Check R^2**: How well does the model explain variance?
        2. **Check p-values**: Are coefficients statistically significant?
        3. **Check residuals**: Are they normally distributed around zero?
        4. **Check diagnostics**: Any issues with autocorrelation or multicollinearity?
        """)
