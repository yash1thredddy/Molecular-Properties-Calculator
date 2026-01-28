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

    property_groups = get_all_property_groups(include_lei)
    all_properties = get_all_properties(include_lei)

    # Session state key for tracking selected properties
    state_key = f"{key_prefix}_selected"

    # Initialize session state for selected properties
    if state_key not in st.session_state:
        st.session_state[state_key] = set()

    # Initialize widget state keys on first run only (not on every render)
    # This prevents the "value set via both default and Session State API" warning
    calc_all_widget_key = f"{key_prefix}_calc_all_widget"
    if calc_all_widget_key not in st.session_state:
        st.session_state[calc_all_widget_key] = False

    # Initialize group and property widget keys
    for group_name, properties in property_groups.items():
        group_widget_key = f"{key_prefix}_group_{group_name}_widget"
        if group_widget_key not in st.session_state:
            st.session_state[group_widget_key] = False
        for prop in properties:
            prop_widget_key = f"{key_prefix}_prop_{prop}_widget"
            if prop_widget_key not in st.session_state:
                st.session_state[prop_widget_key] = False

    # Callback for "Calculate All" checkbox
    def on_calc_all_change():
        if st.session_state[calc_all_widget_key]:
            st.session_state[state_key] = all_properties.copy()
            # Update all group and individual checkbox widget states
            for group_name, props in property_groups.items():
                group_widget_key = f"{key_prefix}_group_{group_name}_widget"
                st.session_state[group_widget_key] = True
                for prop in props:
                    prop_widget_key = f"{key_prefix}_prop_{prop}_widget"
                    st.session_state[prop_widget_key] = True
        else:
            st.session_state[state_key] = set()
            # Update all group and individual checkbox widget states
            for group_name, props in property_groups.items():
                group_widget_key = f"{key_prefix}_group_{group_name}_widget"
                st.session_state[group_widget_key] = False
                for prop in props:
                    prop_widget_key = f"{key_prefix}_prop_{prop}_widget"
                    st.session_state[prop_widget_key] = False

    # Calculate All checkbox - use session state key only (no value= parameter)
    calc_all = st.checkbox(
        "Calculate All Properties",
        key=calc_all_widget_key,
        on_change=on_calc_all_change
    )

    if calc_all:
        # When "Calculate All" is checked, show message and return all
        st.info(f"All {len(all_properties)} properties will be calculated")
        return all_properties.copy()

    # Create callbacks for group checkboxes
    def make_group_callback(group_name: str, props: List[str]):
        def callback():
            widget_key = f"{key_prefix}_group_{group_name}_widget"
            if st.session_state[widget_key]:
                # Add all properties in group
                st.session_state[state_key].update(props)
                # Also update individual checkbox widget states
                for prop in props:
                    prop_widget_key = f"{key_prefix}_prop_{prop}_widget"
                    st.session_state[prop_widget_key] = True
            else:
                # Remove all properties in group
                st.session_state[state_key] -= set(props)
                # Also update individual checkbox widget states
                for prop in props:
                    prop_widget_key = f"{key_prefix}_prop_{prop}_widget"
                    st.session_state[prop_widget_key] = False
        return callback

    # Create callbacks for individual property checkboxes
    def make_prop_callback(prop: str, group_name: str, group_props: List[str]):
        def callback():
            widget_key = f"{key_prefix}_prop_{prop}_widget"
            if st.session_state[widget_key]:
                st.session_state[state_key].add(prop)
            else:
                st.session_state[state_key].discard(prop)
            # Update group checkbox state based on whether all props in group are selected
            all_in_group_selected = all(
                p in st.session_state[state_key] for p in group_props
            )
            group_widget_key = f"{key_prefix}_group_{group_name}_widget"
            st.session_state[group_widget_key] = all_in_group_selected
        return callback

    # Create expandable sections for each property group
    for group_name, properties in property_groups.items():
        is_expanded = group_name in default_expanded

        with st.expander(f"{group_name}", expanded=is_expanded):
            # Sync group widget state with actual selection state
            group_all_selected = all(
                prop in st.session_state[state_key]
                for prop in properties
            )
            group_widget_key = f"{key_prefix}_group_{group_name}_widget"
            # Only update if out of sync (avoid triggering rerun)
            if st.session_state[group_widget_key] != group_all_selected:
                st.session_state[group_widget_key] = group_all_selected

            # Group "Select All" checkbox - use session state key only
            st.checkbox(
                f"Select All {group_name}",
                key=group_widget_key,
                on_change=make_group_callback(group_name, properties)
            )

            # Individual property checkboxes
            cols = st.columns(2 if len(properties) > 4 else 1)

            for i, prop in enumerate(properties):
                col_idx = i % 2 if len(properties) > 4 else 0
                prop_widget_key = f"{key_prefix}_prop_{prop}_widget"

                # Sync widget state with selection state
                is_selected = prop in st.session_state[state_key]
                if st.session_state[prop_widget_key] != is_selected:
                    st.session_state[prop_widget_key] = is_selected

                with cols[col_idx]:
                    # Use session state key only (no value= parameter)
                    st.checkbox(
                        prop.replace('_', ' '),
                        key=prop_widget_key,
                        on_change=make_prop_callback(prop, group_name, properties)
                    )

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
        'NPOLoNHA': 'Polar atoms / Heavy atoms ratio (polarity indicator)',
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
        # Assay Interference Flags
        'PAINS': 'Pan-Assay Interference Substructures (0=Clean, 1=Flagged)',
        'Aggregator': 'Colloidal aggregation risk (0=Clean, 1=Flagged)',
        'Redox': 'Redox-active functional groups (0=Clean, 1=Flagged)',
        'Fluorescence': 'Autofluorescent scaffolds (0=Clean, 1=Flagged)',
        'Thiol': 'Thiol-reactive electrophiles (0=Clean, 1=Flagged)',
    }

    return descriptions.get(property_name, property_name)
