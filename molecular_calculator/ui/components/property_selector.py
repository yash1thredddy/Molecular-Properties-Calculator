"""Property selection component.

This module provides UI components for selecting molecular properties
to calculate, organized by property groups.
"""

import streamlit as st
from typing import Set, Dict, List, Optional

from molecular_calculator.models import PROPERTY_GROUPS, LEI_PROPERTY_GROUP


def get_all_property_groups(include_lei: bool = False) -> Dict[str, List[str]]:
    """Get all property groups.

    Args:
        include_lei: Whether to include LEI properties

    Returns:
        Dictionary mapping group names to property lists
    """
    groups = dict(PROPERTY_GROUPS)
    if include_lei:
        groups.update(LEI_PROPERTY_GROUP)
    return groups


def get_all_properties(include_lei: bool = False) -> Set[str]:
    """Get set of all available properties.

    Args:
        include_lei: Whether to include LEI properties

    Returns:
        Set of all property names
    """
    all_props = set()
    for props in get_all_property_groups(include_lei).values():
        all_props.update(props)
    return all_props


def render_property_selector(
    key_prefix: str = "props",
    include_lei: bool = False,
    default_expanded: Optional[List[str]] = None
) -> Set[str]:
    """Render property selection UI with grouped checkboxes.

    Args:
        key_prefix: Prefix for widget keys
        include_lei: Whether to include LEI properties
        default_expanded: List of group names to expand by default

    Returns:
        Set of selected property names
    """
    if default_expanded is None:
        default_expanded = ["Basic Properties", "Lipinski Properties", "Drug-likeness"]

    # Initialize session state for selections
    state_key = f"{key_prefix}_selected"
    calc_all_key = f"{key_prefix}_calc_all"
    prev_calc_all_key = f"{key_prefix}_prev_calc_all"

    if state_key not in st.session_state:
        st.session_state[state_key] = set()
    if prev_calc_all_key not in st.session_state:
        st.session_state[prev_calc_all_key] = False

    property_groups = get_all_property_groups(include_lei)
    all_properties = get_all_properties(include_lei)

    # Calculate All checkbox
    calc_all = st.checkbox(
        "Calculate All Properties",
        value=st.session_state.get(prev_calc_all_key, False),
        key=calc_all_key
    )

    # Handle toggle changes
    prev_calc_all = st.session_state.get(prev_calc_all_key, False)

    if calc_all and not prev_calc_all:
        # Just checked - select all properties
        st.session_state[state_key] = all_properties.copy()
    elif not calc_all and prev_calc_all:
        # Just unchecked - clear all properties
        st.session_state[state_key] = set()

    # Update previous state
    st.session_state[prev_calc_all_key] = calc_all

    if calc_all:
        # When "Calculate All" is checked, show message and return all
        st.info(f"✅ All {len(all_properties)} properties will be calculated")
        return all_properties.copy()

    # Create expandable sections for each property group
    for group_name, properties in property_groups.items():
        is_expanded = group_name in default_expanded

        with st.expander(f"📊 {group_name}", expanded=is_expanded):
            # Group checkbox to select/deselect all in group
            group_selected = all(
                prop in st.session_state[state_key]
                for prop in properties
            )

            group_check = st.checkbox(
                f"Select All {group_name}",
                value=group_selected,
                key=f"{key_prefix}_group_{group_name}"
            )

            if group_check and not group_selected:
                st.session_state[state_key].update(properties)
            elif not group_check and group_selected:
                st.session_state[state_key] -= set(properties)

            # Individual property checkboxes
            cols = st.columns(2 if len(properties) > 4 else 1)

            for i, prop in enumerate(properties):
                col_idx = i % 2 if len(properties) > 4 else 0

                with cols[col_idx]:
                    checked = st.checkbox(
                        prop.replace('_', ' '),
                        value=prop in st.session_state[state_key],
                        key=f"{key_prefix}_prop_{prop}"
                    )

                    if checked:
                        st.session_state[state_key].add(prop)
                    else:
                        st.session_state[state_key].discard(prop)

    return st.session_state[state_key].copy()


def render_compact_property_selector(
    key_prefix: str = "compact_props",
    include_lei: bool = False
) -> Set[str]:
    """Render a compact property selector using multiselect.

    Args:
        key_prefix: Prefix for widget keys
        include_lei: Whether to include LEI properties

    Returns:
        Set of selected property names
    """
    all_props = sorted(get_all_properties(include_lei))

    selected = st.multiselect(
        "Select Properties to Calculate",
        options=all_props,
        default=[],
        key=f"{key_prefix}_multiselect",
        help="Choose which molecular properties to calculate"
    )

    return set(selected)


def render_property_group_selector(
    key_prefix: str = "group_select",
    include_lei: bool = False
) -> Set[str]:
    """Render a property selector organized by groups using multiselect.

    Args:
        key_prefix: Prefix for widget keys
        include_lei: Whether to include LEI properties

    Returns:
        Set of selected property names
    """
    property_groups = get_all_property_groups(include_lei)

    selected_props = set()

    # Quick select all
    if st.checkbox("Select All Properties", key=f"{key_prefix}_all"):
        return get_all_properties(include_lei)

    for group_name, properties in property_groups.items():
        with st.expander(group_name, expanded=False):
            selected = st.multiselect(
                "Properties",
                options=properties,
                default=[],
                key=f"{key_prefix}_{group_name}",
                label_visibility="collapsed"
            )
            selected_props.update(selected)

    return selected_props


def render_property_summary(selected_properties: Set[str]) -> None:
    """Render a summary of selected properties.

    Args:
        selected_properties: Set of selected property names
    """
    if not selected_properties:
        st.info("No properties selected. Please select at least one property to calculate.")
        return

    st.write(f"**Selected Properties:** {len(selected_properties)}")

    # Group by category
    property_groups = get_all_property_groups(include_lei=True)

    for group_name, properties in property_groups.items():
        group_selected = [p for p in properties if p in selected_properties]
        if group_selected:
            st.caption(f"{group_name}: {', '.join(group_selected)}")


def get_property_description(property_name: str) -> str:
    """Get description for a property.

    Args:
        property_name: Property name

    Returns:
        Description string
    """
    descriptions = {
        'Molecular_Weight': 'Molecular weight in Daltons (Da)',
        'Heavy_Atom_Count': 'Number of non-hydrogen atoms',
        'Atom_Count': 'Total number of atoms including hydrogens',
        'Bond_Count': 'Total number of bonds',
        'Formal_Charge': 'Net formal charge of the molecule',
        'LogP': 'Partition coefficient (lipophilicity), range -2 to 6 typical',
        'HB_Donors': 'Hydrogen bond donors (≤5 for drug-likeness)',
        'HB_Acceptors': 'Hydrogen bond acceptors (≤10 for drug-likeness)',
        'TPSA': 'Topological polar surface area in Ų (≤140 for oral bioavailability)',
        '10xPSA_MW': 'PSA/MW ratio scaled by 10 (membrane permeability indicator)',
        'Rotatable_Bonds': 'Number of rotatable bonds (≤10 for oral bioavailability)',
        'QED': 'Quantitative Estimate of Drug-likeness (0-1 scale)',
        'Aromatic_Rings': 'Number of aromatic rings',
        'Aliphatic_Rings': 'Number of non-aromatic rings',
        'Saturated_Rings': 'Number of saturated rings',
        'Ring_Count': 'Total number of rings',
        'Heteroatoms': 'Number of non-carbon, non-hydrogen atoms',
        'BertzCT': 'Bertz complexity index',
        'Chi0': 'Chi connectivity index 0',
        'Chi1': 'Chi connectivity index 1',
        'CrippenLogP': "Crippen's LogP calculation",
        'CrippenMR': "Crippen's molar refractivity",
        'LabuteASA': "Labute's approximate surface area",
        'Lipinski_Violations': 'Rule of Five violations (0=pass, 1=fail)',
        'Veber_Violations': 'Veber rule violations (0=pass, 1=fail)',
        # LEI properties
        'NSEI': 'Normalized Surface Efficiency Index (pKi / Polar Atoms)',
        'NBEI': 'Normalized Binding Efficiency Index (pKi / Heavy Atoms)',
        'BEI': 'Binding Efficiency Index (pKi / (MW/1000))',
        'SEI': 'Surface Efficiency Index (pKi / (TPSA/100))',
        'nBEI': 'Alternative Binding Efficiency (-log10(Ki / Heavy Atoms))',
        'mBEI': 'Molecular Binding Efficiency (-log10(Ki / MW))',
        'LEH': 'Ligand Efficiency Hopkins (-ΔG / Heavy Atoms)',
        'LEP': 'Ligand Efficiency Polar (-ΔG / Polar Atoms)',
    }

    return descriptions.get(property_name, property_name)
