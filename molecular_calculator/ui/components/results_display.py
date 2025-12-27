"""Results display components.

This module provides components for displaying calculation results,
including property tables, summaries, and statistics.
"""

import streamlit as st
import pandas as pd
from typing import Dict, Any, List, Optional, Set

from molecular_calculator.models import (
    MolecularProperties,
    CalculationResult,
    PROPERTY_GROUPS,
)
from molecular_calculator.ui.components.property_selector import get_property_description


def render_property_card(
    name: str,
    value: Any,
    description: Optional[str] = None
) -> None:
    """Render a single property as a metric card.

    Args:
        name: Property name
        value: Property value
        description: Optional description
    """
    display_name = name.replace('_', ' ')

    if isinstance(value, float):
        formatted_value = f"{value:.3f}"
    else:
        formatted_value = str(value)

    st.metric(
        label=display_name,
        value=formatted_value,
        help=description or get_property_description(name)
    )


def render_properties_grid(
    properties: Dict[str, Any],
    columns: int = 3
) -> None:
    """Render properties in a grid layout.

    Args:
        properties: Dictionary of property names to values
        columns: Number of columns in the grid
    """
    if not properties:
        st.info("No properties to display")
        return

    prop_list = list(properties.items())
    cols = st.columns(columns)

    for i, (name, value) in enumerate(prop_list):
        col_idx = i % columns
        with cols[col_idx]:
            render_property_card(name, value)


def render_properties_by_group(
    properties: Dict[str, Any],
    expanded_groups: Optional[List[str]] = None
) -> None:
    """Render properties organized by group.

    Args:
        properties: Dictionary of property names to values
        expanded_groups: List of group names to expand by default
    """
    if expanded_groups is None:
        expanded_groups = ["Basic Properties", "Lipinski Properties"]

    for group_name, group_props in PROPERTY_GROUPS.items():
        # Filter to properties that are present
        present_props = {
            p: properties[p]
            for p in group_props
            if p in properties
        }

        if not present_props:
            continue

        is_expanded = group_name in expanded_groups

        with st.expander(f"📊 {group_name}", expanded=is_expanded):
            render_properties_grid(present_props, columns=3)


def render_calculation_result(
    result: CalculationResult,
    show_smiles: bool = True
) -> None:
    """Render a calculation result.

    Args:
        result: CalculationResult object
        show_smiles: Whether to show the SMILES string
    """
    if not result.success:
        st.error(f"Calculation failed: {result.error}")
        return

    if show_smiles and result.smiles:
        st.success(f"SMILES: {result.smiles}")

    if result.warnings:
        for warning in result.warnings:
            st.warning(warning)

    if result.properties:
        render_properties_by_group(result.properties.to_dict())


def render_dataframe_preview(
    df: pd.DataFrame,
    max_rows: int = 10,
    title: str = "Data Preview"
) -> None:
    """Render a preview of a DataFrame.

    Args:
        df: DataFrame to preview
        max_rows: Maximum rows to show
        title: Title for the preview
    """
    st.subheader(title)
    st.write(f"Showing {min(len(df), max_rows)} of {len(df):,} rows")
    st.dataframe(df.head(max_rows), use_container_width=True)


def render_batch_results(
    df: pd.DataFrame,
    calculated_columns: List[str],
    show_stats: bool = True
) -> None:
    """Render batch processing results.

    Args:
        df: DataFrame with results
        calculated_columns: List of columns that were calculated
        show_stats: Whether to show statistics
    """
    st.success(f"✅ Processed {len(df):,} molecules")

    # Show calculated columns
    if calculated_columns:
        st.write(f"**Calculated Properties:** {', '.join(calculated_columns)}")

    # Show statistics
    if show_stats and calculated_columns:
        render_batch_statistics(df, calculated_columns)

    # Show data preview
    st.subheader("Results Preview")
    st.dataframe(df.head(20), use_container_width=True)


def render_batch_statistics(
    df: pd.DataFrame,
    numeric_columns: List[str]
) -> None:
    """Render statistics for batch results.

    Args:
        df: DataFrame with results
        numeric_columns: Columns to show statistics for
    """
    # Filter to numeric columns that exist
    valid_cols = [c for c in numeric_columns if c in df.columns]

    if not valid_cols:
        return

    with st.expander("📈 Statistics Summary", expanded=False):
        # Calculate statistics
        stats_df = df[valid_cols].describe()
        st.dataframe(stats_df.round(3), use_container_width=True)

        # Show column-wise stats in metrics
        cols = st.columns(min(4, len(valid_cols)))

        for i, col_name in enumerate(valid_cols[:4]):
            with cols[i % 4]:
                mean_val = df[col_name].mean()
                std_val = df[col_name].std()

                st.metric(
                    label=col_name.replace('_', ' '),
                    value=f"{mean_val:.3f}" if pd.notna(mean_val) else "N/A",
                    delta=f"σ = {std_val:.3f}" if pd.notna(std_val) else None
                )


def render_rule_compliance_summary(
    df: pd.DataFrame,
    lipinski_col: str = "Lipinski_Violations",
    veber_col: str = "Veber_Violations"
) -> None:
    """Render a summary of drug-likeness rule compliance.

    Args:
        df: DataFrame with violation columns
        lipinski_col: Name of Lipinski violations column
        veber_col: Name of Veber violations column
    """
    cols = st.columns(2)

    if lipinski_col in df.columns:
        with cols[0]:
            compliant = (df[lipinski_col] == 0).sum()
            total = df[lipinski_col].notna().sum()
            pct = (compliant / total * 100) if total > 0 else 0

            st.metric(
                label="Lipinski Compliant",
                value=f"{compliant}/{total}",
                delta=f"{pct:.1f}%"
            )

    if veber_col in df.columns:
        with cols[1]:
            compliant = (df[veber_col] == 0).sum()
            total = df[veber_col].notna().sum()
            pct = (compliant / total * 100) if total > 0 else 0

            st.metric(
                label="Veber Compliant",
                value=f"{compliant}/{total}",
                delta=f"{pct:.1f}%"
            )


def render_processing_progress(
    current: int,
    total: int,
    message: str = "Processing molecules..."
) -> None:
    """Render a progress bar for batch processing.

    Args:
        current: Current item number
        total: Total number of items
        message: Progress message
    """
    progress = current / total if total > 0 else 0
    st.progress(progress, text=f"{message} ({current}/{total})")


def render_error_summary(
    errors: List[str],
    max_display: int = 5
) -> None:
    """Render a summary of processing errors.

    Args:
        errors: List of error messages
        max_display: Maximum number of errors to display
    """
    if not errors:
        return

    with st.expander(f"⚠️ {len(errors)} errors occurred", expanded=False):
        for error in errors[:max_display]:
            st.error(error)

        if len(errors) > max_display:
            st.warning(f"... and {len(errors) - max_display} more errors")


def render_property_explanations() -> None:
    """Render property explanations in an expander."""
    with st.expander("ℹ️ Property Explanations"):
        st.markdown("""
### Supported Input Formats
- **SMILES**: Simplified Molecular Input Line Entry System
- **InChI**: International Chemical Identifier
- **InChI Key**: Hashed version of InChI (database lookup required)

### Property Groups

#### Basic Properties
- **Molecular Weight**: Molecular weight in Daltons (Da)
- **Heavy Atom Count**: Number of non-hydrogen atoms
- **Atom Count**: Total number of atoms (including hydrogens)
- **Bond Count**: Total number of bonds
- **Formal Charge**: Net formal charge of the molecule

#### Lipinski Properties (Rule of Five)
- **LogP**: Partition coefficient (lipophilicity, -2 to 6 typical range)
- **HB Donors**: Hydrogen bond donors (≤5 for drug-likeness)
- **HB Acceptors**: Hydrogen bond acceptors (≤10 for drug-likeness)
- **TPSA**: Topological polar surface area in Ų (≤140 for oral bioavailability)
- **10xPSA/MW**: PSA/MW ratio scaled by 10 (membrane permeability indicator)
- **Rotatable Bonds**: Number of rotatable bonds (≤10 for oral bioavailability)

#### Drug-likeness
- **QED**: Quantitative Estimate of Drug-likeness (0-1 scale, higher is better)

#### Rule Violations (Binary: 0=Compliant, 1=Violates)
- **Lipinski Violations**: Passes Lipinski Rule of Five (MW≤500, LogP≤5, HBD≤5, HBA≤10)
- **Veber Violations**: Passes Veber Rule (TPSA≤140, RotBonds≤10)

#### Ring Properties
- **Aromatic Rings**: Number of aromatic rings
- **Aliphatic Rings**: Number of non-aromatic rings
- **Saturated Rings**: Number of saturated rings
- **Ring Count**: Total number of rings
- **Heteroatoms**: Number of non-carbon, non-hydrogen atoms

#### Complexity
- **BertzCT**: Bertz complexity index (higher = more complex)
- **Chi0/Chi1**: Chi connectivity indices (molecular connectivity)

#### Additional Descriptors
- **CrippenLogP**: Crippen's LogP calculation method
- **CrippenMR**: Crippen's molar refractivity
- **LabuteASA**: Labute's approximate surface area
        """)
