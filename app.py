"""Molecular Properties Calculator - Streamlit Application.

A production-ready application for calculating and visualizing
molecular properties from SMILES, InChI, and InChI Key formats.

This is the refactored entry point that uses the modular architecture.
For the legacy monolithic app, see molecular_properties_app.py.

Usage:
    streamlit run app.py

Developed by: Yashwanth Reddy for ITR-UIC
Part of: Chemo-Informatics Toolkit
"""

import streamlit as st
import warnings
import numpy as np

# Suppress numpy runtime warnings for division by zero
warnings.filterwarnings('ignore', category=RuntimeWarning, module='numpy')
np.seterr(divide='ignore', invalid='ignore')

# Import from refactored package
from molecular_calculator import MolecularCalculator, __version__
from molecular_calculator.config.settings import config
from molecular_calculator.utils.session_state import SessionState

# Import UI components
from molecular_calculator.ui.pages import (
    render_single_molecule_page,
    render_batch_processing_page,
    render_3d_regression_page,
    render_3d_regression_help,
)


def main():
    """Main application entry point."""
    # Page configuration
    st.set_page_config(
        page_title=config.APP_NAME,
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )

    # Initialize session state
    SessionState.init_defaults()

    # Title
    st.title(f"🧬 {config.APP_NAME}")
    st.caption(f"Version {__version__}")

    # Sidebar settings
    _render_sidebar()

    # Get current mode
    input_mode = st.session_state.get('input_mode', 'Single Molecule')
    enable_online_lookup = st.session_state.get('enable_online_lookup', True)

    # Render appropriate page
    if input_mode == "Single Molecule":
        render_single_molecule_page(enable_online_lookup=enable_online_lookup)

    elif input_mode == "Batch Processing":
        render_batch_processing_page(enable_online_lookup=enable_online_lookup)

    elif input_mode == "Data Visualization":
        _render_visualization_page()

    elif input_mode == "3D Regression Analysis":
        render_3d_regression_page()

    # Footer
    _render_footer()


def _render_sidebar():
    """Render the sidebar with settings and navigation."""
    st.sidebar.title("Navigation")

    # Settings section
    st.sidebar.subheader("⚙️ Settings")

    suppress_warnings = st.sidebar.checkbox(
        "Suppress RDKit warnings",
        value=True,
        help="Hide stereochemistry conflict warnings",
        key="suppress_warnings"
    )

    enable_online_lookup = st.sidebar.checkbox(
        "Enable InChI Key conversion",
        value=True,
        help="Convert InChI Keys using online databases (NIH CIR, PubChem)",
        key="enable_online_lookup"
    )

    # Apply settings
    MolecularCalculator.suppress_rdkit_warnings(suppress_warnings)

    # Mode selection
    st.sidebar.header("📊 Input Options")

    input_mode = st.sidebar.radio(
        "Select input mode:",
        [
            "Single Molecule",
            "Batch Processing",
            "Data Visualization",
            "3D Regression Analysis"
        ],
        key="input_mode"
    )

    # About section
    st.sidebar.markdown("---")
    st.sidebar.markdown("""
    ### About

    **Molecular Properties Calculator** calculates physicochemical
    properties from molecular structures using RDKit.

    **Supported Formats:**
    - SMILES
    - InChI
    - InChI Key (online lookup)

    **Developed by:**
    Yashwanth Reddy for ITR-UIC
    """)


def _render_visualization_page():
    """Render the data visualization page.

    Uses the shared visualization component for consistency with batch processing.
    """
    st.header("Data Visualization")

    st.info("""
    📊 **Visualization Mode**

    Upload a dataset with calculated properties to create interactive visualizations.

    **Available Charts:**
    - Scatter plots with trendlines and regression statistics
    - 3D OLS Regression with fitted plane
    - Histograms with mean/median lines
    - Box plots and violin plots with grouping
    - Correlation heatmaps
    - Line plots and bar charts

    **Features:**
    - Size by, Shape by, Color by options
    - Trendlines with R² and equation
    - Color scale selection with reverse option
    - Export charts as images
    """)

    # Import and use visualization components
    from molecular_calculator.ui.components import (
        render_file_upload_section,
        render_distribution_plots,
        render_interactive_visualization,
    )

    # File upload
    df, smiles_col, name_col = render_file_upload_section(key_prefix="dataviz")

    if df is None:
        st.info("👆 Upload a file to get started with visualization")
        return

    # Quick Distribution Analysis
    render_distribution_plots(df, key_prefix="dataviz_dist")

    # Interactive Visualization - using shared component with unique key prefix
    render_interactive_visualization(
        df,
        key_prefix="dataviz_viz",
        smiles_col=smiles_col,
        name_col=name_col
    )


def _render_footer():
    """Render the application footer."""
    st.markdown("---")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("**Developed by:** Yashwanth Reddy for ITR-UIC")

    with col2:
        st.markdown(f"**Version:** {__version__}")

    with col3:
        st.markdown("**Powered by:** RDKit, Streamlit, Plotly")


if __name__ == "__main__":
    main()
