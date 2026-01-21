"""
Assay Interference Flags Display Component

This module provides UI components for displaying assay interference detection
results in a visually appealing format matching the design specification.
"""

import streamlit as st
import pandas as pd
from typing import Dict, List, Optional, Any

from molecular_calculator.services.assay_interference import (
    InterferenceFlags,
    FLAG_DESCRIPTIONS,
    get_interference_flags_from_smiles,
)


def render_interference_header() -> None:
    """Render the header section for assay interference flags."""
    st.markdown("### Assay Interference Flags")
    st.markdown(
        "*Detection of compounds with known assay interference mechanisms*",
        help="Based on Bisson et al. (2016) and Baell & Holloway (2010)"
    )


def render_flag_metric(
    flag_name: str,
    count: int,
    total: int,
    is_flagged: bool = False
) -> None:
    """
    Render a single flag metric with status indicator.

    Args:
        flag_name: Name of the flag (PAINS, Aggregator, etc.)
        count: Number of compounds with this flag
        total: Total number of compounds
        is_flagged: Whether this flag is raised (for single molecule)
    """
    flag_info = FLAG_DESCRIPTIONS.get(flag_name, {})
    full_name = flag_info.get('full_name', flag_name)
    description = flag_info.get('description', '')

    # Calculate percentage
    percentage = (count / total * 100) if total > 0 else 0

    # Determine status
    if count == 0:
        status_html = '<span style="color: #51cf66;">&#8593; Clean</span>'
    else:
        status_html = f'<span style="color: #ff6b6b;">&#8593; {percentage:.0f}%</span>'

    # Render metric
    st.markdown(f"**{flag_name}**")
    st.markdown(f"<h2 style='margin: 0; padding: 0;'>{count}</h2>", unsafe_allow_html=True)
    st.markdown(status_html, unsafe_allow_html=True)


def render_interference_metrics(flags: InterferenceFlags) -> None:
    """
    Render interference flag metrics for a single molecule.

    Args:
        flags: InterferenceFlags object with detection results
    """
    cols = st.columns(5)

    flag_names = ['PAINS', 'Aggregator', 'Redox', 'Fluorescence', 'Thiol']
    flag_values = [flags.pains, flags.aggregator, flags.redox, flags.fluorescence, flags.thiol]

    for col, name, is_flagged in zip(cols, flag_names, flag_values):
        with col:
            flag_info = FLAG_DESCRIPTIONS.get(name, {})
            full_name = flag_info.get('full_name', name)

            # Status indicator
            if is_flagged:
                status = f'<span style="color: #ff6b6b;">&#9679; Flagged</span>'
                value = "1"
            else:
                status = f'<span style="color: #51cf66;">&#9679; Clean</span>'
                value = "0"

            st.markdown(f"**{name}**")
            st.markdown(f"<h2 style='margin: 0; padding: 0;'>{value}</h2>", unsafe_allow_html=True)
            st.markdown(status, unsafe_allow_html=True)


def render_batch_interference_metrics(
    flags_df: pd.DataFrame,
    total_compounds: int
) -> None:
    """
    Render interference flag metrics for batch processing results.

    Args:
        flags_df: DataFrame with interference flag columns
        total_compounds: Total number of compounds processed
    """
    cols = st.columns(5)

    flag_columns = ['PAINS', 'Aggregator', 'Redox', 'Fluorescence', 'Thiol']

    for col, flag_name in zip(cols, flag_columns):
        with col:
            if flag_name in flags_df.columns:
                count = int(flags_df[flag_name].sum())
            else:
                count = 0

            percentage = (count / total_compounds * 100) if total_compounds > 0 else 0

            # Status indicator
            if count == 0:
                status = '<span style="color: #51cf66;">&#8593; Clean</span>'
            else:
                status = f'<span style="color: #ff6b6b;">&#8593; {percentage:.0f}%</span>'

            st.markdown(f"**{flag_name}**")
            st.markdown(f"<h2 style='margin: 0; padding: 0;'>{count}</h2>", unsafe_allow_html=True)
            st.markdown(status, unsafe_allow_html=True)


def render_flag_summary_table(flags_df: pd.DataFrame, total_compounds: int) -> None:
    """
    Render a summary table of flag counts and percentages.

    Args:
        flags_df: DataFrame with interference flag columns
        total_compounds: Total number of compounds
    """
    st.markdown("#### Flag Summary")

    flag_data = []
    flag_columns = ['PAINS', 'Aggregator', 'Redox', 'Fluorescence', 'Thiol']
    descriptions = [
        'Pan-Assay Interference',
        'Colloidal Aggregation',
        'Redox Cycling',
        'Fluorescence Interference',
        'Thiol Reactivity'
    ]
    colors = ['#ff6b6b', '#ffa94d', '#ffd43b', '#74c0fc', '#b197fc']

    for flag_name, desc, color in zip(flag_columns, descriptions, colors):
        if flag_name in flags_df.columns:
            count = int(flags_df[flag_name].sum())
        else:
            count = 0

        percentage = (count / total_compounds * 100) if total_compounds > 0 else 0

        flag_data.append({
            'Flag': flag_name,
            'Count': count,
            '%': f"{percentage:.0f}%",
            'Description': desc
        })

    summary_df = pd.DataFrame(flag_data)

    # Style the dataframe
    st.dataframe(
        summary_df,
        use_container_width=True,
        hide_index=True,
        column_config={
            'Flag': st.column_config.TextColumn('Flag', width='small'),
            'Count': st.column_config.NumberColumn('Count', width='small'),
            '%': st.column_config.TextColumn('%', width='small'),
            'Description': st.column_config.TextColumn('Description', width='large'),
        }
    )


def render_flagged_compounds_table(
    df: pd.DataFrame,
    id_column: str = 'ChEMBL_ID',
    name_column: Optional[str] = None,
    max_rows: int = 100
) -> None:
    """
    Render a table showing flagged compounds with their flags.

    Args:
        df: DataFrame with compound data and flag columns
        id_column: Column name for compound ID
        name_column: Optional column name for compound name
        max_rows: Maximum rows to display
    """
    st.markdown("#### Flagged Compounds")

    flag_columns = ['PAINS', 'Aggregator', 'Redox', 'Fluorescence', 'Thiol']
    available_flags = [col for col in flag_columns if col in df.columns]

    if not available_flags:
        st.info("No interference flag data available.")
        return

    # Filter to only flagged compounds
    mask = df[available_flags].sum(axis=1) > 0
    flagged_df = df[mask].copy()

    if len(flagged_df) == 0:
        st.success("No flagged compounds detected!")
        return

    # Determine which flag(s) each compound has
    def get_flags(row):
        flags = []
        for flag in available_flags:
            if row.get(flag, 0) == 1:
                flags.append(flag)
        return ', '.join(flags)

    flagged_df['Flag'] = flagged_df.apply(get_flags, axis=1)

    # Select columns to display
    display_cols = []
    if id_column in flagged_df.columns:
        display_cols.append(id_column)
    display_cols.append('Flag')
    if name_column and name_column in flagged_df.columns:
        display_cols.append(name_column)

    # Rename for display
    display_df = flagged_df[display_cols].head(max_rows)

    column_config = {}
    if id_column in display_cols:
        column_config[id_column] = st.column_config.TextColumn(id_column, width='medium')
    column_config['Flag'] = st.column_config.TextColumn('Flag', width='medium')
    if name_column and name_column in display_cols:
        column_config[name_column] = st.column_config.TextColumn('Molecule', width='medium')

    st.dataframe(
        display_df,
        use_container_width=True,
        hide_index=True,
        column_config=column_config
    )

    if len(flagged_df) > max_rows:
        st.caption(f"Showing {max_rows} of {len(flagged_df)} flagged compounds")


def render_interference_section(
    flags: Optional[InterferenceFlags] = None,
    flags_df: Optional[pd.DataFrame] = None,
    total_compounds: int = 1,
    df: Optional[pd.DataFrame] = None,
    id_column: str = 'ID',
    name_column: Optional[str] = None,
    show_details: bool = True
) -> None:
    """
    Render the complete assay interference flags section.

    Can be used for both single molecule and batch processing results.

    Args:
        flags: InterferenceFlags object for single molecule
        flags_df: DataFrame with flag columns for batch results
        total_compounds: Total number of compounds
        df: Full DataFrame with compound data (for flagged compounds table)
        id_column: Column name for compound ID
        name_column: Column name for compound name
        show_details: Whether to show detailed tables
    """
    render_interference_header()

    # Single molecule mode
    if flags is not None:
        render_interference_metrics(flags)

        if show_details and flags.total_flags > 0:
            with st.expander("View Details", expanded=False):
                details = []
                if flags.pains and flags.pains_details:
                    details.append(f"**PAINS:** {', '.join(flags.pains_details)}")
                if flags.aggregator and flags.aggregator_reason:
                    details.append(f"**Aggregator:** {flags.aggregator_reason}")
                if flags.redox and flags.redox_groups:
                    details.append(f"**Redox:** {', '.join(flags.redox_groups)}")
                if flags.fluorescence and flags.fluorescence_scaffolds:
                    details.append(f"**Fluorescence:** {', '.join(flags.fluorescence_scaffolds)}")
                if flags.thiol and flags.thiol_electrophiles:
                    details.append(f"**Thiol:** {', '.join(flags.thiol_electrophiles)}")

                for detail in details:
                    st.markdown(detail)

    # Batch mode
    elif flags_df is not None:
        render_batch_interference_metrics(flags_df, total_compounds)

        if show_details:
            col1, col2 = st.columns(2)

            with col1:
                render_flag_summary_table(flags_df, total_compounds)

            with col2:
                if df is not None:
                    render_flagged_compounds_table(df, id_column, name_column)


def calculate_batch_interference_flags(
    df: pd.DataFrame,
    smiles_column: str
) -> pd.DataFrame:
    """
    Calculate interference flags for all compounds in a DataFrame.

    Args:
        df: DataFrame with SMILES column
        smiles_column: Name of the column containing SMILES

    Returns:
        DataFrame with added interference flag columns
    """
    result_df = df.copy()

    # Initialize flag columns
    flag_columns = ['PAINS', 'Aggregator', 'Redox', 'Fluorescence', 'Thiol']
    for col in flag_columns:
        result_df[col] = 0

    # Calculate flags for each compound
    for idx, row in result_df.iterrows():
        smiles = row.get(smiles_column, '')
        if smiles and pd.notna(smiles):
            flags = get_interference_flags_from_smiles(str(smiles))
            result_df.at[idx, 'PAINS'] = 1 if flags.pains else 0
            result_df.at[idx, 'Aggregator'] = 1 if flags.aggregator else 0
            result_df.at[idx, 'Redox'] = 1 if flags.redox else 0
            result_df.at[idx, 'Fluorescence'] = 1 if flags.fluorescence else 0
            result_df.at[idx, 'Thiol'] = 1 if flags.thiol else 0

    return result_df
