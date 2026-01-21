"""Single molecule analysis page.

This module provides the single molecule analysis page for the Streamlit app.
"""

import streamlit as st
from typing import Optional

from molecular_calculator.core import MolecularCalculator
from molecular_calculator.models import InputFormat
from molecular_calculator.ui.components import (
    render_smiles_input,
    render_format_selector,
    render_molecule_info,
    render_property_selector,
    render_properties_by_group,
    render_property_explanations,
    create_download_button,
    render_interference_section,
)
from molecular_calculator.services.assay_interference import (
    get_interference_flags_from_smiles,
)


def render_single_molecule_page(
    enable_online_lookup: bool = True
) -> None:
    """Render the single molecule analysis page.

    Args:
        enable_online_lookup: Whether to enable InChI Key online lookup
    """
    st.header("Single Molecule Analysis")

    # Initialize calculator
    calculator = MolecularCalculator()

    # Initialize session state
    if 'single_molecule_analyzed' not in st.session_state:
        st.session_state.single_molecule_analyzed = False
    if 'single_current_smiles' not in st.session_state:
        st.session_state.single_current_smiles = None
    if 'single_selected_properties' not in st.session_state:
        st.session_state.single_selected_properties = set()

    # Input molecule
    molecule_input = render_smiles_input(
        label="Enter molecular structure:",
        placeholder="SMILES, InChI, or InChI Key",
        help_text="Paste your molecular structure here - analysis will happen automatically",
        key="single_molecule_input"
    )

    # Reset button
    if st.session_state.single_molecule_analyzed:
        col1, col2 = st.columns([3, 1])
        with col2:
            if st.button("🔄 Reset", help="Clear analysis and start over"):
                st.session_state.single_molecule_analyzed = False
                st.session_state.single_current_smiles = None
                st.session_state.single_selected_properties = set()
                st.rerun()

    # Process input
    if molecule_input and molecule_input.strip():
        input_text = molecule_input.strip()

        # Auto-detect format
        detected_format = calculator.conversion_service.detect_format(input_text)
        st.info(f"Detected format: {detected_format.value.upper()}")

        # Manual format selection
        format_options = ["smiles", "inchi", "inchi_key"]
        format_map = {
            InputFormat.SMILES: 0,
            InputFormat.INCHI: 1,
            InputFormat.INCHI_KEY: 2,
            InputFormat.UNKNOWN: 0
        }

        input_format = st.selectbox(
            "Confirm input format:",
            options=format_options,
            index=format_map.get(detected_format, 0),
            key="single_format_select"
        )

        # Convert to SMILES if needed
        smiles = None

        if input_format == 'inchi_key' and not enable_online_lookup:
            st.warning(
                "InChI Key conversion is disabled. "
                "Please enable 'InChI Key conversion' in Settings."
            )
        else:
            # Convert
            if input_format == 'inchi_key':
                with st.spinner("Converting InChI Key using online databases..."):
                    smiles = calculator.convert_to_smiles(
                        input_text, input_format, enable_online_lookup
                    )
            else:
                smiles = calculator.convert_to_smiles(
                    input_text, input_format, enable_online_lookup
                )

        # Handle conversion result
        if smiles is None and input_format == 'inchi_key' and enable_online_lookup:
            st.error(
                "❌ Could not resolve this InChI Key. "
                "Please verify the format or check your internet connection."
            )

        if smiles:
            st.success(f"✅ SMILES: {smiles}")
            st.session_state.single_molecule_analyzed = True
            st.session_state.single_current_smiles = smiles
        else:
            st.session_state.single_molecule_analyzed = False
            st.session_state.single_current_smiles = None

    elif molecule_input and not molecule_input.strip():
        st.warning("Please enter a valid molecular structure.")

    # Show property selection if molecule is analyzed
    if st.session_state.single_molecule_analyzed and st.session_state.single_current_smiles:
        st.subheader("Select Properties to Calculate")

        # Property selection
        selected_properties = render_property_selector(
            key_prefix="single",
            include_lei=False,
            default_expanded=["Basic Properties", "Lipinski Properties", "Drug-likeness"]
        )

        st.session_state.single_selected_properties = selected_properties

        # Calculate button
        if st.button("Calculate Properties", type="primary", key="single_calculate_btn"):
            if not selected_properties:
                st.warning("Please select at least one property to calculate.")
            else:
                _calculate_and_display_properties(
                    st.session_state.single_current_smiles,
                    selected_properties
                )

    # Property explanations
    render_property_explanations()


def _calculate_and_display_properties(
    smiles: str,
    selected_properties: set
) -> None:
    """Calculate and display properties for a molecule.

    Args:
        smiles: SMILES string
        selected_properties: Set of property names to calculate
    """
    calculator = MolecularCalculator()

    with st.spinner("Calculating properties..."):
        result = calculator.calculate(smiles)

    if not result.success:
        st.error(f"Calculation failed: {result.error}")
        return

    # Filter to selected properties
    all_props = result.properties.to_dict()
    filtered_props = {
        k: v for k, v in all_props.items()
        if k in selected_properties
    }

    if not filtered_props:
        st.warning("No properties were calculated for this molecule.")
        return

    # Display results
    st.subheader("Calculated Properties")
    render_properties_by_group(filtered_props)

    # Calculate and display interference flags
    st.markdown("---")
    interference_flags = get_interference_flags_from_smiles(smiles)
    render_interference_section(flags=interference_flags, show_details=True)

    # Download options
    st.markdown("---")
    st.subheader("Export Results")

    import pandas as pd
    results_df = pd.DataFrame([{
        'SMILES': smiles,
        **filtered_props,
        **interference_flags.to_dict()
    }])

    col1, col2 = st.columns(2)

    with col1:
        create_download_button(
            results_df,
            filename="molecular_properties.csv",
            label="📥 Download CSV",
            key="single_download_csv"
        )

    with col2:
        # Create Excel download
        from molecular_calculator.ui.components import create_excel_download_button
        create_excel_download_button(
            results_df,
            filename="molecular_properties.xlsx",
            label="📥 Download Excel",
            key="single_download_xlsx"
        )
