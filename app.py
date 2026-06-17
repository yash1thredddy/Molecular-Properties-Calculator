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
from molecular_calculator.config import setup_logging

# Initialize logging (INFO level by default, can be configured via environment)
import os
_log_level = os.environ.get('LOG_LEVEL', 'INFO')
setup_logging(level=_log_level)
from molecular_calculator.utils.session_state import SessionState

# Import UI components
from molecular_calculator.ui.pages import (
    render_single_molecule_page,
    render_batch_processing_page,
    render_3d_regression_page,
    render_gmm_page,
)

# Plain-language documentation site (GitHub Pages).
DOCS_URL = "https://yash1thredddy.github.io/molecular-properties-calculator-docs/"


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

    # Suppress RDKit warnings by default
    MolecularCalculator.suppress_rdkit_warnings(True)

    # Title
    st.title(f"🧬 {config.APP_NAME}")

    # About section
    with st.expander("About", expanded=False):
        st.markdown(f"""
        **Version:** {__version__}

        Calculate physicochemical properties from molecular structures using RDKit.

        **Supported Input Formats:**
        - SMILES (Simplified Molecular Input Line Entry System)
        - InChI (International Chemical Identifier)
        - InChI Key (with automatic online lookup via NIH/PubChem)

        **Features:**
        - Single molecule analysis with detailed property breakdown
        - Batch processing for CSV/Excel files
        - Interactive data visualization with multiple chart types
        - 3D OLS regression analysis for structure-activity relationships
        - Gaussian Mixture Model (GMM) grouping with plain-language results

        **Developed by:** Yashwanth Reddy for ITR-UIC
        """)

    # Navigation — grouped, icon'd pages via Streamlit's multipage API.
    # Each st.Page wraps a page-render callable; the menu renders in the sidebar.
    pg = st.navigation({
        "Analyze": [
            st.Page(_single_molecule_page, title="Single Molecule",
                    icon=":material/science:", default=True),
            st.Page(_batch_processing_page, title="Batch Processing",
                    icon=":material/dataset:"),
            st.Page(_render_visualization_page, title="Data Visualization",
                    icon=":material/insights:"),
        ],
        "Modeling": [
            st.Page(render_3d_regression_page, title="3D Regression Analysis",
                    icon=":material/deployed_code:"),
            st.Page(render_gmm_page, title="GMM Analysis",
                    icon=":material/bubble_chart:"),
        ],
    })

    # Supplementary info below the nav menu.
    _render_sidebar()

    # Render the selected page, then the shared footer.
    pg.run()
    _render_footer()


def _single_molecule_page():
    """Page wrapper: single-molecule analysis with online lookup enabled."""
    render_single_molecule_page(enable_online_lookup=True)


def _batch_processing_page():
    """Page wrapper: batch processing with online lookup enabled."""
    render_batch_processing_page(enable_online_lookup=True)


def _render_sidebar():
    """Render supplementary info beneath the navigation menu."""
    with st.sidebar:
        # Documentation — placed directly under the navigation menu (filling the
        # empty gap above "Supported Formats") and styled primary so it stands out.
        # Opens the plain-language docs site in a new tab.
        st.link_button(
            "Documentation",
            DOCS_URL,
            icon=":material/menu_book:",
            type="primary",
            width="stretch",
        )
        st.markdown("---")
        st.markdown("""
        **Supported Formats:**
        - SMILES
        - InChI
        - InChI Key (online lookup)
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
