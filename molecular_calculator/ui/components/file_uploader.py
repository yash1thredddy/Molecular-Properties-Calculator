"""File upload component for molecular data.

This module provides a reusable file upload component with:
- CSV and XLSX support
- File validation
- SMILES column detection
- Column mapping
"""

import streamlit as st
import pandas as pd
from typing import Optional, Tuple, List
from io import BytesIO

from molecular_calculator.utils.validators import FileValidator, DataFrameValidator
from molecular_calculator.utils.session_state import SessionState
from molecular_calculator.config.settings import config


def render_file_uploader(
    key: str = "file_uploader",
    help_text: str = "Upload CSV or XLSX file with molecular structures"
) -> Optional[st.runtime.uploaded_file_manager.UploadedFile]:
    """Render a file upload widget.

    Args:
        key: Unique key for the widget
        help_text: Help text to display

    Returns:
        Uploaded file object or None
    """
    return st.file_uploader(
        "Upload molecular data file",
        type=['csv', 'xlsx'],
        accept_multiple_files=False,  # explicit: app processes one file at a time
        help=help_text,
        key=key
    )


def validate_uploaded_file(
    uploaded_file
) -> Tuple[bool, Optional[str], List[str]]:
    """Validate an uploaded file.

    Args:
        uploaded_file: Streamlit uploaded file object

    Returns:
        Tuple of (is_valid, error_message, warnings)
    """
    result = FileValidator.validate_upload(uploaded_file)

    if not result.is_valid:
        return False, result.errors[0] if result.errors else "Invalid file", result.warnings

    return True, None, result.warnings


def read_uploaded_file(
    uploaded_file,
    show_errors: bool = True
) -> Optional[pd.DataFrame]:
    """Read an uploaded file into a DataFrame.

    Args:
        uploaded_file: Streamlit uploaded file object
        show_errors: Whether to show error messages in Streamlit

    Returns:
        DataFrame or None if reading fails
    """
    if uploaded_file is None:
        return None

    try:
        filename = uploaded_file.name.lower()

        if filename.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
        elif filename.endswith('.xlsx'):
            df = pd.read_excel(uploaded_file, engine='openpyxl')
        else:
            if show_errors:
                st.error(f"Unsupported file type: {filename}")
            return None

        return df

    except Exception as e:
        if show_errors:
            st.error(f"Error reading file: {str(e)}")
        return None


def detect_smiles_column(df: pd.DataFrame) -> Optional[str]:
    """Auto-detect SMILES column in DataFrame.

    Args:
        df: DataFrame to analyze

    Returns:
        Column name or None if not found
    """
    possible_names = [
        'smiles', 'SMILES', 'Smiles', 'smi', 'SMI',
        'canonical_smiles', 'CANONICAL_SMILES',
        'CanonicalSMILES', 'Canonical_SMILES'
    ]

    for col in df.columns:
        if col in possible_names:
            return col
        if col.lower() in [name.lower() for name in possible_names]:
            return col

    return None


def detect_name_column(df: pd.DataFrame, smiles_col: str) -> Optional[str]:
    """Auto-detect name/ID column in DataFrame.

    Args:
        df: DataFrame to analyze
        smiles_col: SMILES column to exclude

    Returns:
        Column name or None if not found
    """
    name_patterns = [
        'Id', 'ID', 'id', 'Name', 'name', 'NAME',
        'Compound', 'compound', 'COMPOUND',
        'Molecule', 'molecule', 'MOLECULE',
        'Title', 'title', 'TITLE'
    ]

    for col_candidate in name_patterns:
        if col_candidate in df.columns and col_candidate != smiles_col:
            return col_candidate

    # Fallback to first non-numeric object column
    for col in df.columns:
        if col != smiles_col and df[col].dtype == 'object':
            return col

    return None


def render_column_selector(
    df: pd.DataFrame,
    label: str = "Select SMILES column",
    detected_col: Optional[str] = None,
    key: str = "smiles_col"
) -> Optional[str]:
    """Render a column selector dropdown.

    Args:
        df: DataFrame with columns to select from
        label: Label for the dropdown
        detected_col: Pre-detected column name
        key: Unique key for the widget

    Returns:
        Selected column name
    """
    columns = list(df.columns)

    # Set default index
    default_idx = 0
    if detected_col and detected_col in columns:
        default_idx = columns.index(detected_col)

    return st.selectbox(
        label,
        options=columns,
        index=default_idx,
        key=key
    )


def render_file_upload_section(
    key_prefix: str = "batch"
) -> Tuple[Optional[pd.DataFrame], Optional[str], Optional[str]]:
    """Render complete file upload section with validation.

    This is a convenience function that combines file upload,
    validation, reading, and column detection.

    Args:
        key_prefix: Prefix for widget keys

    Returns:
        Tuple of (DataFrame, smiles_column, name_column)
    """
    uploaded_file = render_file_uploader(key=f"{key_prefix}_uploader")

    if uploaded_file is None:
        return None, None, None

    # Validate file
    is_valid, error, warnings = validate_uploaded_file(uploaded_file)

    if not is_valid:
        st.error(f"❌ {error}")
        return None, None, None

    for warning in warnings:
        st.warning(warning)

    # Check if file changed - only clear results if truly a new file
    if SessionState.file_changed(uploaded_file):
        # Clear cached results - but do NOT clear uploaded_file_hash (we just set it)
        # Only clear specific result keys, not all batch state
        SessionState.clear('batch_results_df')
        SessionState.clear('batch_calculated_columns')
        SessionState.clear('batch_results_smiles_col')
        SessionState.clear('batch_results_name_col')

    # Read file
    df = read_uploaded_file(uploaded_file)

    if df is None:
        return None, None, None

    # Validate DataFrame
    df_result = DataFrameValidator.validate(df)

    if not df_result.is_valid:
        st.error(f"❌ {df_result.errors[0]}")
        return None, None, None

    for warning in df_result.warnings:
        st.warning(warning)

    # Show file info
    st.success(f"✅ Loaded {len(df):,} rows, {len(df.columns)} columns")

    # Show data preview
    with st.expander("📋 Data Preview (first 50 rows)", expanded=True):
        st.dataframe(df.head(50), width='stretch')

    # Detect and select columns
    detected_smiles = detect_smiles_column(df)
    detected_name = None

    col1, col2 = st.columns(2)

    with col1:
        smiles_col = render_column_selector(
            df,
            label="SMILES Column",
            detected_col=detected_smiles,
            key=f"{key_prefix}_smiles_col"
        )

    with col2:
        if smiles_col:
            detected_name = detect_name_column(df, smiles_col)

        name_col = render_column_selector(
            df,
            label="Name/ID Column (optional)",
            detected_col=detected_name,
            key=f"{key_prefix}_name_col"
        )

    return df, smiles_col, name_col


def create_download_button(
    df: pd.DataFrame,
    filename: str = "results.csv",
    label: str = "Download CSV",
    key: str = "download_csv"
) -> None:
    """Create a download button for DataFrame.

    Args:
        df: DataFrame to download
        filename: Default filename
        label: Button label
        key: Unique key for the widget
    """
    csv = df.to_csv(index=False)
    st.download_button(
        label=label,
        data=csv,
        file_name=filename,
        mime="text/csv",
        key=key
    )


def create_excel_download_button(
    df: pd.DataFrame,
    filename: str = "results.xlsx",
    label: str = "Download Excel",
    key: str = "download_xlsx"
) -> None:
    """Create a download button for DataFrame as Excel.

    Args:
        df: DataFrame to download
        filename: Default filename
        label: Button label
        key: Unique key for the widget
    """
    buffer = BytesIO()
    df.to_excel(buffer, index=False, engine='openpyxl')
    buffer.seek(0)

    st.download_button(
        label=label,
        data=buffer,
        file_name=filename,
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key=key
    )
