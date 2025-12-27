"""Molecule viewer component.

This module provides components for viewing molecular structures,
wrapping the existing structure_viewer_component.
"""

import streamlit as st
import streamlit.components.v1 as components
from typing import Optional, Dict, Any

# Import the existing structure viewer
try:
    from molecular_calculator.ui.components.structure_viewer import (
        get_structure_viewer_component,
        get_structure_viewer_hint
    )
    STRUCTURE_VIEWER_AVAILABLE = True
except ImportError:
    STRUCTURE_VIEWER_AVAILABLE = False


def render_structure_viewer_hint() -> None:
    """Render a hint message about clicking points to view structures."""
    if not STRUCTURE_VIEWER_AVAILABLE:
        return

    hint_html = get_structure_viewer_hint()
    st.markdown(hint_html, unsafe_allow_html=True)


def embed_structure_viewer(
    chart_id: str = "plotly_chart",
    x_col: Optional[str] = None,
    y_col: Optional[str] = None,
    z_col: Optional[str] = None,
    name_col: Optional[str] = None,
    height: int = 0
) -> None:
    """Embed the structure viewer component.

    This component attaches to Plotly charts and shows molecular structures
    when points are clicked.

    Args:
        chart_id: Unique identifier for the chart
        x_col: Name of X-axis column
        y_col: Name of Y-axis column
        z_col: Name of Z-axis column (for 3D charts)
        name_col: Name of the molecule name/ID column
        height: Height of the component (0 for minimal)
    """
    if not STRUCTURE_VIEWER_AVAILABLE:
        st.warning("Structure viewer component not available")
        return

    viewer_html = get_structure_viewer_component(
        chart_id=chart_id,
        x_col=x_col,
        y_col=y_col,
        z_col=z_col,
        name_col=name_col
    )

    components.html(viewer_html, height=height)


def render_smiles_input(
    label: str = "Enter molecular structure",
    placeholder: str = "SMILES, InChI, or InChI Key",
    help_text: str = "Paste your molecular structure here",
    key: str = "molecule_input"
) -> str:
    """Render a SMILES/molecular structure input field.

    Args:
        label: Label for the input
        placeholder: Placeholder text
        help_text: Help text
        key: Unique key for the widget

    Returns:
        Input string
    """
    return st.text_area(
        label,
        placeholder=placeholder,
        help=help_text,
        key=key
    )


def render_format_selector(
    detected_format: str = "smiles",
    key: str = "format_selector"
) -> str:
    """Render a format selector dropdown.

    Args:
        detected_format: Auto-detected format
        key: Unique key for the widget

    Returns:
        Selected format string
    """
    formats = ["smiles", "inchi", "inchi_key"]

    # Determine default index
    default_idx = 0
    if detected_format in formats:
        default_idx = formats.index(detected_format)

    return st.selectbox(
        "Confirm input format:",
        options=formats,
        index=default_idx,
        key=key
    )


def render_molecule_info(
    smiles: str,
    detected_format: str,
    original_input: Optional[str] = None
) -> None:
    """Render molecule information display.

    Args:
        smiles: SMILES string
        detected_format: Detected/selected format
        original_input: Original input if converted
    """
    st.info(f"Detected format: {detected_format.upper()}")

    if smiles:
        st.success(f"SMILES: {smiles}")

        if original_input and original_input != smiles:
            with st.expander("View original input"):
                st.code(original_input)


def prepare_chart_customdata(
    df,
    smiles_col: str,
    name_col: Optional[str] = None
) -> tuple:
    """Prepare customdata for Plotly charts with structure viewer.

    The structure viewer expects customdata in format:
    - [SMILES, name, index] if name_col is provided
    - [SMILES, index] otherwise

    Args:
        df: DataFrame with data
        smiles_col: Name of SMILES column
        name_col: Optional name/ID column

    Returns:
        Tuple of (modified_df, customdata_columns)
    """
    import pandas as pd

    # Add row index if not present
    if '_row_index' not in df.columns:
        df = df.copy()
        df['_row_index'] = range(len(df))

    # Build customdata columns list
    if name_col and name_col in df.columns:
        customdata_cols = [smiles_col, name_col, '_row_index']
    else:
        customdata_cols = [smiles_col, '_row_index']

    return df, customdata_cols


def get_structure_viewer_guide() -> str:
    """Get the structure viewer usage guide.

    Returns:
        Markdown string with usage guide
    """
    return """
### Interactive Structure Viewer Guide

**How to Use:**
1. Create any scatter plot, 3D scatter, or other point-based visualization
2. Click on any data point in the chart
3. A side panel will slide in showing:
   - 2D molecular structure (rendered by SmilesDrawer)
   - Molecule name/ID (if available)
   - SMILES string
   - Point coordinates

**Requirements:**
- Your data must have a SMILES column
- The chart must be created with Plotly
- customdata must include SMILES as the first element

**Troubleshooting:**
- If structures don't appear, check that SMILES values are valid
- Some complex structures may take a moment to render
- If the panel doesn't open, try refreshing the page

**Browser Compatibility:**
- Works in Chrome, Firefox, Edge, and Safari
- Some tracking prevention settings may block the SmilesDrawer library
"""
