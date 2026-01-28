"""Batch processing page.

This module provides the batch processing page for the Streamlit app.
"""

import streamlit as st
import pandas as pd
import numpy as np
from typing import Optional, Set, Dict, List

from molecular_calculator.core import MolecularCalculator
from molecular_calculator.ui.components import (
    render_file_upload_section,
    render_property_selector,
    render_batch_results,
    render_rule_compliance_summary,
    create_download_button,
    create_excel_download_button,
    render_property_explanations,
    render_distribution_plots,
    render_interactive_visualization,
    render_interference_section,
)
from molecular_calculator.utils.session_state import SessionState
from molecular_calculator.services.assay_interference import get_interference_flags_from_smiles

# Import LEI functionality
from molecular_calculator.services.ligand_efficiency import (
    DependencyChecker,
    LigandEfficiencyCalculator,
    get_lei_descriptions,
)


def render_batch_processing_page(
    enable_online_lookup: bool = True
) -> None:
    """Render the batch processing page.

    Args:
        enable_online_lookup: Whether to enable InChI Key online lookup
    """
    st.header("Batch Processing")

    st.write("""
    Upload a CSV or Excel file containing molecular structures (SMILES, InChI, or InChI Key).
    Select the properties you want to calculate and process all molecules at once.
    """)

    # File upload section
    df, smiles_col, name_col = render_file_upload_section(key_prefix="batch")

    if df is None or smiles_col is None:
        st.info("👆 Upload a file to get started")
        render_property_explanations()
        return

    # Initialize LEI session state
    if 'batch_selected_leis' not in st.session_state:
        st.session_state.batch_selected_leis = set()

    # Property selection
    st.subheader("Select Properties to Calculate")

    selected_properties = render_property_selector(
        key_prefix="batch",
        include_lei=False,
        default_expanded=["Basic Properties", "Lipinski Properties"]
    )

    # LEI Section
    st.markdown("---")
    selected_leis, pki_col = _render_lei_section(df, smiles_col)

    # Process button
    st.markdown("---")
    st.subheader("🚀 Ready to Process?")

    has_selections = len(selected_properties) > 0 or len(selected_leis) > 0

    # Show selection summary
    if has_selections:
        items = []
        if selected_properties:
            items.append(f"{len(selected_properties)} properties")
        if selected_leis:
            items.append(f"{len(selected_leis)} LEIs")
        st.info(f"📊 **Ready to calculate:** {' + '.join(items)}")
    else:
        st.warning("Please select at least one property or LEI to calculate.")

    col1, col2, col3 = st.columns([1, 1, 2])

    with col1:
        process_btn = st.button(
            "🚀 Calculate Properties",
            type="primary",
            disabled=not has_selections,
            key="batch_process_btn"
        )

    with col2:
        use_parallel = st.checkbox(
            "⚡ Parallel Processing",
            value=len(df) > 50,
            key="batch_use_parallel",
            help="Use parallel processing for faster calculation (recommended for >50 molecules)"
        )

    # Process batch
    if process_btn and has_selections:
        results_df = _process_batch(
            df,
            smiles_col,
            selected_properties,
            selected_leis,
            pki_col,
            enable_online_lookup,
            use_parallel=use_parallel
        )

        if results_df is not None:
            SessionState.set('batch_results_df', results_df)
            # Store which columns were calculated
            all_calculated = set(selected_properties) | set(selected_leis)
            SessionState.set('batch_calculated_columns', all_calculated)
            # Store column names for display persistence (use different key names to avoid widget conflict)
            SessionState.set('batch_results_smiles_col', smiles_col)
            SessionState.set('batch_results_name_col', name_col)

    # Display results if available
    results_df = SessionState.get('batch_results_df')

    if results_df is not None:
        calculated_columns = SessionState.get('batch_calculated_columns', set())
        # Use stored column names to ensure consistency across reruns
        stored_smiles_col = SessionState.get('batch_results_smiles_col', smiles_col)
        stored_name_col = SessionState.get('batch_results_name_col', name_col)
        _display_batch_results(results_df, calculated_columns, stored_smiles_col, stored_name_col)


def _render_lei_section(df: pd.DataFrame, smiles_col: str) -> tuple:
    """Render the LEI (Ligand Efficiency Indices) section.

    Args:
        df: The uploaded DataFrame
        smiles_col: Name of the SMILES column

    Returns:
        Tuple of (Set of selected LEI names, pKi column name or None)
    """
    st.subheader("⚗️ Ligand Efficiency Indices (LEI)")
    st.markdown("*Based on AtlasCBS methodology (Cele Abad-Zapatero, A. Cortes-Cabrera 2013) - requires pKi column*")

    # Auto-detect pKi column
    pki_detected = DependencyChecker.detect_column(df, 'pki')

    # Get all numeric columns for manual selection
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

    # pKi column selection
    col1, col2 = st.columns([2, 1])

    with col1:
        # Build options list with auto-detected at top
        pki_options = ['None (No LEI calculations)']
        if pki_detected:
            pki_options.append(f"{pki_detected} (Auto-detected)")
        pki_options.extend([c for c in numeric_cols if c != pki_detected])

        pki_selection = st.selectbox(
            "📊 Select pKi/pIC50 Column:",
            options=pki_options,
            index=1 if pki_detected else 0,
            key="batch_pki_col",
            help="Select the column containing pKi or pIC50 values for LEI calculations"
        )

    with col2:
        st.markdown("**Supported names:**")
        st.caption("pKi, pki, pIC50, pic50, p_Ki")

    # Parse selection
    if pki_selection == 'None (No LEI calculations)':
        st.info("ℹ️ **No pKi column selected** - Select a pKi column to enable LEI calculations")
        return set(), None

    # Extract actual column name from selection
    if '(Auto-detected)' in pki_selection:
        pki_col = pki_detected
    else:
        pki_col = pki_selection

    st.success(f"✅ pKi column selected: '{pki_col}'")

    # Get detected columns
    detected_cols = DependencyChecker.detect_all_columns(df)
    # Override with user-selected columns
    detected_cols['smiles'] = smiles_col
    detected_cols['pki'] = pki_col

    has_smiles = smiles_col is not None
    has_pki = pki_col is not None
    can_calculate_all = has_smiles and has_pki

    # Show column detection status
    with st.expander("📊 Column Detection Status", expanded=False):
        col1, col2 = st.columns(2)

        with col1:
            st.markdown("**Available Columns:**")
            for dep_name, col_name in detected_cols.items():
                if col_name:
                    st.write(f"✓ {dep_name.upper()}: `{col_name}`")

        with col2:
            st.markdown("**Missing Columns:**")
            missing = [dep for dep, col in detected_cols.items() if col is None]
            if missing:
                for dep in missing:
                    if dep == 'smiles':
                        st.write(f"✗ {dep.upper()} (required for calculating missing properties)")
                    else:
                        st.write(f"✗ {dep.upper()}")
            else:
                st.write("All standard columns detected!")

    # LEI property selection
    st.markdown("**Select LEI Properties to Calculate:**")

    lei_descriptions = get_lei_descriptions()
    all_leis = list(lei_descriptions.keys())

    # Check which LEIs can be calculated
    lei_check = DependencyChecker.check_lei_dependencies(df, all_leis)

    if can_calculate_all:
        can_calculate_leis = set(all_leis)
        st.info("✨ SMILES and pKi detected - all LEIs enabled (missing properties will be calculated from SMILES)")
    else:
        can_calculate_leis = set(lei_check.get('can_calculate', []))

    # Create three columns for LEI selection
    lei_col1, lei_col2, lei_col3 = st.columns(3)

    for idx, lei in enumerate(all_leis):
        col = [lei_col1, lei_col2, lei_col3][idx % 3]

        with col:
            can_calculate = lei in can_calculate_leis
            status_icon = "✅" if can_calculate else "⚠️"

            lei_selected = lei in st.session_state.batch_selected_leis
            checkbox_label = f"{status_icon} {lei}"

            if st.checkbox(
                checkbox_label,
                value=lei_selected,
                key=f"lei_{lei}",
                help=lei_descriptions[lei],
                disabled=not can_calculate
            ):
                st.session_state.batch_selected_leis.add(lei)
            else:
                st.session_state.batch_selected_leis.discard(lei)

    # Show calculation plan
    if st.session_state.batch_selected_leis:
        selected_lei_check = DependencyChecker.check_lei_dependencies(
            df, list(st.session_state.batch_selected_leis)
        )

        with st.expander("ℹ️ Calculation Plan", expanded=False):
            status_msg = DependencyChecker.generate_status_message(selected_lei_check)
            st.markdown(status_msg)

        st.info(f"📋 Selected LEIs: {len(st.session_state.batch_selected_leis)}")

    return st.session_state.batch_selected_leis.copy(), pki_col


def _process_batch(
    df: pd.DataFrame,
    smiles_col: str,
    selected_properties: Set[str],
    selected_leis: Set[str],
    pki_col: Optional[str],
    enable_online_lookup: bool,
    use_parallel: bool = True
) -> Optional[pd.DataFrame]:
    """Process a batch of molecules.

    Args:
        df: Input DataFrame
        smiles_col: Name of SMILES column
        selected_properties: Set of properties to calculate
        selected_leis: Set of LEI names to calculate
        pki_col: Name of pKi column for LEI calculations
        enable_online_lookup: Whether to enable InChI Key lookup
        use_parallel: Whether to use parallel processing (faster for large datasets)

    Returns:
        DataFrame with calculated properties
    """
    progress_bar = st.progress(0, text="Processing molecules...")
    status_text = st.empty()

    try:
        total = len(df)

        # Calculate molecular properties
        if selected_properties:
            status_text.text("Calculating molecular properties...")

            if use_parallel and total > 50:
                # Use parallel processing for larger datasets
                def progress_callback(completed, total_count):
                    progress = completed / total_count
                    progress_bar.progress(progress, text=f"Processing molecule {completed}/{total_count}...")

                final_df = MolecularCalculator.process_batch_parallel(
                    df,
                    smiles_col,
                    selected_properties=selected_properties,
                    enable_online_lookup=enable_online_lookup,
                    progress_callback=progress_callback
                )
            else:
                # Use sequential processing for small datasets
                results = []
                for i, (idx, row) in enumerate(df.iterrows()):
                    smiles = row[smiles_col]

                    if pd.isna(smiles):
                        results.append({})
                        continue

                    # Auto-detect and convert
                    smiles_str = str(smiles)
                    props = MolecularCalculator.calculate_molecular_properties(smiles_str)

                    # Filter to selected properties
                    if props:
                        props = {k: v for k, v in props.items() if k in selected_properties}

                    results.append(props)

                    # Update progress
                    progress = (i + 1) / total
                    progress_bar.progress(progress, text=f"Processing molecule {i + 1}/{total}...")

                # Create results DataFrame
                results_df = pd.DataFrame(results)

                # Combine with original DataFrame
                final_df = pd.concat([df.reset_index(drop=True), results_df], axis=1)
        else:
            final_df = df.copy()

        # Calculate LEIs if selected
        if selected_leis and pki_col:
            status_text.text("Calculating Ligand Efficiency Indices...")
            progress_bar.progress(0.9, text="Calculating LEIs...")

            # Build manual mappings using user-selected columns
            manual_mappings = {
                'smiles': smiles_col,
                'pki': pki_col
            }

            lei_result_df, lei_status = LigandEfficiencyCalculator.process_batch(
                final_df,
                list(selected_leis),
                show_errors=False,
                manual_mappings=manual_mappings
            )

            if lei_status['success']:
                final_df = lei_result_df
                status_text.text(f"LEI calculation complete! Calculated: {', '.join(lei_status['calculated_leis'])}")
            else:
                st.warning(f"⚠️ LEI calculation issue: {lei_status.get('message', 'Unknown error')}")

        # Calculate interference flags if any are selected
        interference_props = {'PAINS', 'Aggregator', 'Redox', 'Fluorescence', 'Thiol'}
        selected_interference = selected_properties & interference_props

        if selected_interference and smiles_col in final_df.columns:
            status_text.text("Calculating assay interference flags...")
            progress_bar.progress(0.95, text="Calculating assay interference flags...")

            # Calculate flags for each molecule (use detailed dict for pattern info)
            flag_results = []
            empty_result = {
                'PAINS': 0, 'PAINS_Details': '',
                'Aggregator': 0, 'Aggregator_Details': '',
                'Redox': 0, 'Redox_Details': '',
                'Fluorescence': 0, 'Fluorescence_Details': '',
                'Thiol': 0, 'Thiol_Details': ''
            }
            for smiles in final_df[smiles_col]:
                if pd.isna(smiles):
                    flag_results.append(empty_result.copy())
                else:
                    flags = get_interference_flags_from_smiles(str(smiles))
                    flag_results.append(flags.to_detailed_dict())

            # Add flag columns and detail columns to dataframe
            flags_df = pd.DataFrame(flag_results)
            for col in flags_df.columns:
                # Add flag columns if selected
                if col in selected_interference:
                    final_df[col] = flags_df[col]
                # Also add detail columns for selected flags
                elif col.endswith('_Details'):
                    base_flag = col.replace('_Details', '')
                    if base_flag in selected_interference:
                        final_df[col] = flags_df[col]

        progress_bar.empty()
        status_text.empty()

        st.success(f"✅ Successfully processed {total} molecules")

        return final_df

    except Exception as e:
        progress_bar.empty()
        status_text.empty()
        st.error(f"Error processing batch: {str(e)}")
        return None


def _display_batch_results(
    df: pd.DataFrame,
    calculated_columns: Set[str],
    smiles_col: Optional[str] = None,
    name_col: Optional[str] = None
) -> None:
    """Display batch processing results.

    Args:
        df: Results DataFrame
        calculated_columns: Set of calculated column names
        smiles_col: Optional SMILES column for structure viewer
        name_col: Optional name/ID column for display
    """
    st.subheader("Results")

    # Filter to columns that exist
    present_cols = [c for c in calculated_columns if c in df.columns]

    # Show compliance summary
    if 'Lipinski_Violations' in df.columns or 'Veber_Violations' in df.columns:
        render_rule_compliance_summary(df)

    # Check if any interference properties were selected
    interference_props = {'PAINS', 'Aggregator', 'Redox', 'Fluorescence', 'Thiol'}
    selected_interference = calculated_columns & interference_props

    # Show Assay Interference Flags section only if selected
    # Note: Flags are already calculated during _process_batch(), no need to recalculate
    if selected_interference:
        st.markdown("---")
        flag_columns = ['PAINS', 'Aggregator', 'Redox', 'Fluorescence', 'Thiol']
        available_flag_columns = [col for col in flag_columns if col in df.columns]

        if available_flag_columns:
            # Use pre-calculated flags from batch processing
            render_interference_section(
                flags_df=df[available_flag_columns] if len(available_flag_columns) == len(flag_columns) else df,
                total_compounds=len(df),
                df=df,
                id_column=name_col if name_col else (smiles_col if smiles_col else df.columns[0]),
                name_column=name_col,
                show_details=True
            )
        else:
            st.warning("Interference flags were selected but not found in results. Please re-run the calculation.")

    # Show FULL results table FIRST (most important for users)
    st.subheader("📋 Results Table")
    st.write(f"Showing all **{len(df):,}** molecules with **{len(present_cols)}** calculated properties")
    st.dataframe(df, width='stretch', height=400)

    # Download options right after the table
    st.subheader("📥 Export Results")

    col1, col2 = st.columns(2)

    with col1:
        create_download_button(
            df,
            filename="batch_results.csv",
            label="📥 Download CSV",
            key="batch_download_csv"
        )

    with col2:
        create_excel_download_button(
            df,
            filename="batch_results.xlsx",
            label="📥 Download Excel",
            key="batch_download_xlsx"
        )

    # Show statistics
    if present_cols:
        st.subheader("📊 Statistics Summary")
        # Only use numeric columns for statistics
        numeric_cols = [c for c in present_cols if df[c].dtype in ['float64', 'int64', 'float32', 'int32']]
        if numeric_cols:
            stats = df[numeric_cols].describe()
            st.dataframe(stats.round(3), width='stretch')

    # Quick Distribution Analysis
    st.markdown("---")
    render_distribution_plots(df, calculated_columns, key_prefix="batch_dist")

    # Interactive Visualization
    st.markdown("---")
    render_interactive_visualization(
        df,
        key_prefix="batch_viz",
        smiles_col=smiles_col,
        name_col=name_col
    )

    # Property explanations
    st.markdown("---")
    render_property_explanations()
