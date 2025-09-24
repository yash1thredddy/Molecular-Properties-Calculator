import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski, rdMolDescriptors
import plotly.graph_objects as go
import plotly.express as px
from io import StringIO
import re

def convert_to_smiles(input_text, input_type):
    """Convert InChI or InChI Key to SMILES"""
    try:
        if input_type.lower() == 'smiles':
            return input_text
        elif input_type.lower() == 'inchi':
            mol = Chem.MolFromInchi(input_text)
            return Chem.MolToSmiles(mol) if mol else None
        elif input_type.lower() == 'inchi_key':
            # InChI Key cannot be directly converted to SMILES without additional data
            st.error("Direct conversion from InChI Key to SMILES is not possible without a database lookup.")
            return None
    except:
        return None

def calculate_molecular_properties(smiles):
    """Calculate various molecular properties from SMILES"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}

        properties = {}

        # Basic properties
        properties['Molecular_Weight'] = round(Descriptors.MolWt(mol), 3)
        properties['Heavy_Atom_Count'] = Descriptors.HeavyAtomCount(mol)
        properties['Atom_Count'] = mol.GetNumAtoms()
        properties['Bond_Count'] = mol.GetNumBonds()

        # Lipinski Rule of Five
        properties['LogP'] = round(Descriptors.MolLogP(mol), 3)
        properties['HB_Donors'] = Descriptors.NumHDonors(mol)
        properties['HB_Acceptors'] = Descriptors.NumHAcceptors(mol)
        properties['TPSA'] = round(Descriptors.TPSA(mol), 2)
        properties['Rotatable_Bonds'] = Descriptors.NumRotatableBonds(mol)

        # Drug-likeness
        properties['QED'] = round(QED.qed(mol), 4)

        # Ring properties
        properties['Aromatic_Rings'] = Descriptors.NumAromaticRings(mol)
        properties['Aliphatic_Rings'] = Descriptors.NumAliphaticRings(mol)
        properties['Saturated_Rings'] = Descriptors.NumSaturatedRings(mol)
        properties['Heteroatoms'] = Descriptors.NumHeteroatoms(mol)

        # Crippen descriptors
        properties['CrippenLogP'] = round(Crippen.MolLogP(mol), 3)
        properties['CrippenMR'] = round(Crippen.MolMR(mol), 3)

        # Try to calculate additional descriptors with error handling
        try:
            properties['BertzCT'] = round(Descriptors.BertzCT(mol), 3)
        except:
            properties['BertzCT'] = None

        try:
            properties['LabuteASA'] = round(Descriptors.LabuteASA(mol), 3)
        except:
            properties['LabuteASA'] = None

        # Try rdMolDescriptors with fallback
        try:
            from rdkit.Chem import rdMolDescriptors
            properties['Chi0'] = round(rdMolDescriptors.Chi0(mol), 3)
            properties['Chi1'] = round(rdMolDescriptors.Chi1(mol), 3)
        except:
            properties['Chi0'] = None
            properties['Chi1'] = None

        # Additional simple descriptors
        try:
            properties['Formal_Charge'] = Chem.rdmolops.GetFormalCharge(mol)
        except:
            properties['Formal_Charge'] = 0

        # Fragment-based descriptors
        try:
            properties['Ring_Count'] = rdMolDescriptors.CalcNumRings(mol) if hasattr(rdMolDescriptors, 'CalcNumRings') else None
        except:
            properties['Ring_Count'] = None

        # Lipinski Rule compliance (0 = compliant, 1 = violates rule)
        lipinski_violations = sum([
            properties['Molecular_Weight'] > 500,
            properties['LogP'] > 5,
            properties['HB_Donors'] > 5,
            properties['HB_Acceptors'] > 10
        ])
        properties['Lipinski_Violations'] = 1 if lipinski_violations > 0 else 0

        # Veber Rule compliance (0 = compliant, 1 = violates rule)
        veber_violations = sum([
            properties['TPSA'] > 140,
            properties['Rotatable_Bonds'] > 10
        ])
        properties['Veber_Violations'] = 1 if veber_violations > 0 else 0

        # Filter out None values
        properties = {k: v for k, v in properties.items() if v is not None}

        return properties
    except Exception as e:
        st.error(f"Error calculating properties: {str(e)}")
        return {}

def detect_smiles_column(df):
    """Detect SMILES column with case-insensitive matching"""
    possible_names = ['smiles', 'SMILES', 'Smiles', 'smi', 'SMI', 'canonical_smiles', 'CANONICAL_SMILES']
    for col in df.columns:
        if col in possible_names or col.lower() in [name.lower() for name in possible_names]:
            return col
    return None

def detect_input_format(input_text):
    """Auto-detect input format"""
    if input_text.startswith('InChI='):
        return 'inchi'
    elif re.match(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$', input_text):
        return 'inchi_key'
    else:
        return 'smiles'

# Streamlit App
st.set_page_config(
    page_title="Molecular Properties Calculator",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("🧪 Molecular Properties Calculator")
st.markdown("Calculate various chemical properties from molecular structures")

# Sidebar for input options
st.sidebar.header("Input Options")
input_mode = st.sidebar.radio("Select input mode:", ["Single Molecule", "Batch Processing"])

if input_mode == "Single Molecule":
    st.header("Single Molecule Analysis")

    # Input molecule
    molecule_input = st.text_area("Enter molecular structure:", placeholder="SMILES, InChI, or InChI Key")

    if molecule_input:
        # Auto-detect format
        detected_format = detect_input_format(molecule_input.strip())
        st.info(f"Detected format: {detected_format.upper()}")

        # Manual format selection
        input_format = st.selectbox("Confirm input format:",
                                   ["smiles", "inchi", "inchi_key"],
                                   index=0 if detected_format == "smiles" else 1 if detected_format == "inchi" else 2)

        # Convert to SMILES if needed
        smiles = convert_to_smiles(molecule_input.strip(), input_format)

        if smiles:
            st.success(f"SMILES: {smiles}")

            # Property selection
            st.subheader("Select Properties to Calculate")

            # Define all properties organized by groups
            property_groups = {
                "Basic Properties": ['Molecular_Weight', 'Heavy_Atom_Count', 'Atom_Count', 'Bond_Count', 'Formal_Charge'],
                "Lipinski Properties": ['LogP', 'HB_Donors', 'HB_Acceptors', 'TPSA', 'Rotatable_Bonds'],
                "Drug-likeness": ['QED'],
                "Rule Violations": ['Lipinski_Violations', 'Veber_Violations'],
                "Ring Properties": ['Aromatic_Rings', 'Aliphatic_Rings', 'Saturated_Rings', 'Ring_Count', 'Heteroatoms'],
                "Complexity": ['BertzCT', 'Chi0', 'Chi1'],
                "Additional": ['CrippenLogP', 'CrippenMR', 'LabuteASA']
            }

            # Initialize session state for property selection
            if 'selected_properties' not in st.session_state:
                st.session_state.selected_properties = set()

            # Calculate All checkbox
            calc_all = st.checkbox("Calculate All Properties", value=False)

            if calc_all:
                st.session_state.selected_properties = set()
                for props in property_groups.values():
                    st.session_state.selected_properties.update(props)

            if not calc_all:
                # Create expandable sections for each property group
                for group_name, properties in property_groups.items():
                    with st.expander(f"📊 {group_name}", expanded=True if group_name in ["Basic Properties", "Lipinski Properties", "Drug-likeness"] else False):

                        # Group checkbox to select/deselect all in group
                        group_selected = all(prop in st.session_state.selected_properties for prop in properties)
                        group_check = st.checkbox(f"Select All {group_name}", value=group_selected, key=f"group_{group_name}")

                        if group_check and not group_selected:
                            # If group checkbox is checked, add all properties
                            st.session_state.selected_properties.update(properties)
                        elif not group_check and group_selected:
                            # If group checkbox is unchecked, remove all properties
                            st.session_state.selected_properties.difference_update(properties)

                        # Individual property checkboxes
                        cols = st.columns(2 if len(properties) > 4 else 1)
                        for i, prop in enumerate(properties):
                            col_idx = i % len(cols)
                            with cols[col_idx]:
                                # Format property name for display
                                display_name = prop.replace('_', ' ').title()
                                prop_selected = prop in st.session_state.selected_properties

                                if st.checkbox(display_name, value=prop_selected, key=f"prop_{prop}"):
                                    st.session_state.selected_properties.add(prop)
                                else:
                                    st.session_state.selected_properties.discard(prop)

            # Show selected properties count
            if not calc_all:
                if len(st.session_state.selected_properties) == 0:
                    st.warning("No properties selected. Please select at least one property or check 'Calculate All Properties'.")
                else:
                    st.info(f"Selected properties: {len(st.session_state.selected_properties)}")
            else:
                total_props = sum(len(props) for props in property_groups.values())
                st.info(f"All properties selected: {total_props}")

            if st.button("Calculate Properties", type="primary"):
                if calc_all or st.session_state.selected_properties:
                    properties = calculate_molecular_properties(smiles)

                    if properties:
                        # Display results based on selection
                        results_df = pd.DataFrame([properties]).T
                        results_df.columns = ['Value']

                        # Filter properties based on selection
                        if not calc_all and st.session_state.selected_properties:
                            # Only show selected properties that actually exist in the results
                            available_selected = [prop for prop in st.session_state.selected_properties if prop in results_df.index]
                            if available_selected:
                                results_df = results_df.loc[available_selected]

                        st.subheader("Calculated Properties")
                        st.dataframe(results_df, use_container_width=True)

                        # Create visualization for numeric properties
                        numeric_props = results_df.select_dtypes(include=[np.number])
                        if not numeric_props.empty:
                            fig = px.bar(
                                x=numeric_props.index,
                                y=numeric_props.iloc[:, 0],
                                title="Molecular Properties"
                            )
                            fig.update_layout(xaxis_title="Property", yaxis_title="Value")
                            st.plotly_chart(fig, use_container_width=True)

                        # Download results
                        csv = results_df.to_csv()
                        st.download_button(
                            label="Download Results as CSV",
                            data=csv,
                            file_name=f"molecular_properties_{smiles.replace('/', '_')[:20]}.csv",
                            mime="text/csv"
                        )
                    else:
                        st.error("Could not calculate properties for this molecule.")
                else:
                    st.warning("Please select at least one property to calculate.")
        else:
            st.error("Could not process the input. Please check the format.")

else:  # Batch Processing
    st.header("Batch Processing")

    uploaded_file = st.file_uploader("Upload file", type=['csv', 'xlsx'])

    if uploaded_file is not None:
        # Handle different file formats
        file_extension = uploaded_file.name.split('.')[-1].lower()

        try:
            if file_extension == 'csv':
                df = pd.read_csv(uploaded_file)
            elif file_extension == 'xlsx':
                df = pd.read_excel(uploaded_file)
            else:
                st.error("Unsupported file format. Please upload CSV or XLSX files.")
                st.stop()
        except Exception as e:
            st.error(f"Error reading file: {str(e)}")
            st.stop()

        st.subheader("Data Preview")
        st.info(f"Loaded {len(df)} rows from {uploaded_file.name}")
        st.dataframe(df.head(), use_container_width=True)

        # Detect SMILES column
        smiles_col = detect_smiles_column(df)

        if smiles_col:
            st.success(f"Detected SMILES column: '{smiles_col}'")
        else:
            st.warning("Could not detect SMILES column automatically")
            smiles_col = st.selectbox("Select SMILES column:", df.columns)

        # Property selection for batch
        st.subheader("Select Properties to Calculate")

        # Use same property groups as single molecule
        property_groups_batch = {
            "Basic Properties": ['Molecular_Weight', 'Heavy_Atom_Count', 'Atom_Count', 'Bond_Count', 'Formal_Charge'],
            "Lipinski Properties": ['LogP', 'HB_Donors', 'HB_Acceptors', 'TPSA', 'Rotatable_Bonds'],
            "Drug-likeness": ['QED'],
            "Rule Violations": ['Lipinski_Violations', 'Veber_Violations'],
            "Ring Properties": ['Aromatic_Rings', 'Aliphatic_Rings', 'Saturated_Rings', 'Ring_Count', 'Heteroatoms'],
            "Complexity": ['BertzCT', 'Chi0', 'Chi1'],
            "Additional": ['CrippenLogP', 'CrippenMR', 'LabuteASA']
        }

        # Initialize batch session state
        if 'batch_selected_properties' not in st.session_state:
            st.session_state.batch_selected_properties = set()

        # Calculate All checkbox for batch
        calc_all_batch = st.checkbox("Calculate All Properties (Batch)", value=False)

        if calc_all_batch:
            st.session_state.batch_selected_properties = set()
            for props in property_groups_batch.values():
                st.session_state.batch_selected_properties.update(props)

        if not calc_all_batch:
            # Create expandable sections for batch processing
            for group_name, properties in property_groups_batch.items():
                with st.expander(f"📊 {group_name}", expanded=True if group_name in ["Basic Properties", "Lipinski Properties", "Drug-likeness"] else False):

                    # Group checkbox
                    group_selected = all(prop in st.session_state.batch_selected_properties for prop in properties)
                    group_check = st.checkbox(f"Select All {group_name}", value=group_selected, key=f"batch_group_{group_name}")

                    if group_check and not group_selected:
                        st.session_state.batch_selected_properties.update(properties)
                    elif not group_check and group_selected:
                        st.session_state.batch_selected_properties.difference_update(properties)

                    # Individual property checkboxes
                    cols = st.columns(2 if len(properties) > 4 else 1)
                    for i, prop in enumerate(properties):
                        col_idx = i % len(cols)
                        with cols[col_idx]:
                            display_name = prop.replace('_', ' ').title()
                            prop_selected = prop in st.session_state.batch_selected_properties

                            if st.checkbox(display_name, value=prop_selected, key=f"batch_prop_{prop}"):
                                st.session_state.batch_selected_properties.add(prop)
                            else:
                                st.session_state.batch_selected_properties.discard(prop)

        # Show selected properties count for batch
        if not calc_all_batch:
            if len(st.session_state.batch_selected_properties) == 0:
                st.warning("No properties selected. Please select at least one property or check 'Calculate All Properties (Batch)'.")
            else:
                st.info(f"Selected properties: {len(st.session_state.batch_selected_properties)}")
        else:
            total_props = sum(len(props) for props in property_groups_batch.values())
            st.info(f"All properties selected: {total_props}")

        if st.button("Process Batch", type="primary"):
            if smiles_col and (calc_all_batch or st.session_state.batch_selected_properties):
                progress_bar = st.progress(0)
                status_text = st.empty()

                results = []
                total_rows = len(df)

                for idx, row in df.iterrows():
                    smiles = row[smiles_col]

                    # Auto-detect and convert if needed
                    input_format = detect_input_format(str(smiles)) if pd.notna(smiles) else 'smiles'
                    if input_format != 'smiles':
                        smiles = convert_to_smiles(str(smiles), input_format)

                    properties = calculate_molecular_properties(smiles) if pd.notna(smiles) else {}

                    # Filter properties based on selection - only if properties is not None/empty
                    if properties and not calc_all_batch and st.session_state.batch_selected_properties:
                        properties = {k: v for k, v in properties.items() if k in st.session_state.batch_selected_properties}

                    results.append(properties)

                    # Update progress
                    progress = (idx + 1) / total_rows
                    progress_bar.progress(progress)
                    status_text.text(f"Processing: {idx + 1}/{total_rows}")

                # Create results DataFrame
                results_df = pd.DataFrame(results)
                final_df = pd.concat([df, results_df], axis=1)

                status_text.text("Processing complete!")

                st.subheader("Results")
                st.dataframe(final_df, use_container_width=True)

                # Summary statistics
                if not results_df.empty:
                    st.subheader("Summary Statistics")
                    summary = results_df.describe()
                    st.dataframe(summary, use_container_width=True)

                # Download results
                csv = final_df.to_csv(index=False)
                st.download_button(
                    label="Download Results as CSV",
                    data=csv,
                    file_name="molecular_properties_batch_results.csv",
                    mime="text/csv"
                )
            else:
                st.error("Please select a SMILES column and at least one property.")

# Information section
with st.expander("ℹ️ Information & Property Explanations"):
    st.markdown("""
    ### Supported Input Formats
    - **SMILES**: Simplified Molecular Input Line Entry System
    - **InChI**: International Chemical Identifier
    - **InChI Key**: Hashed version of InChI (database lookup required)

    ### File Formats
    - **CSV**: Comma-separated values
    - **XLSX**: Excel spreadsheet format

    ### Property Explanations

    #### Basic Properties
    - **Molecular_Weight**: Molecular weight in Daltons (Da)
    - **Heavy_Atom_Count**: Number of non-hydrogen atoms
    - **Atom_Count**: Total number of atoms (including hydrogens)
    - **Bond_Count**: Total number of bonds
    - **Formal_Charge**: Net formal charge of the molecule

    #### Lipinski Properties (Rule of Five)
    - **LogP**: Partition coefficient (lipophilicity, -2 to 6 typical range)
    - **HB_Donors**: Hydrogen bond donors (≤5 for drug-likeness)
    - **HB_Acceptors**: Hydrogen bond acceptors (≤10 for drug-likeness)
    - **TPSA**: Topological polar surface area in Ų (≤140 for oral bioavailability)
    - **Rotatable_Bonds**: Number of rotatable bonds (≤10 for oral bioavailability)

    #### Drug-likeness
    - **QED**: Quantitative Estimate of Drug-likeness (0-1 scale, higher is better)

    #### Rule Violations (Binary: 0=Compliant, 1=Violates)
    - **Lipinski_Violations**:
      - 0 = Passes Lipinski Rule of Five (MW≤500, LogP≤5, HBD≤5, HBA≤10)
      - 1 = Violates at least one Lipinski rule
    - **Veber_Violations**:
      - 0 = Passes Veber Rule (TPSA≤140, RotBonds≤10)
      - 1 = Violates at least one Veber rule

    #### Ring Properties
    - **Aromatic_Rings**: Number of aromatic rings
    - **Aliphatic_Rings**: Number of non-aromatic rings
    - **Saturated_Rings**: Number of saturated rings
    - **Ring_Count**: Total number of rings
    - **Heteroatoms**: Number of non-carbon, non-hydrogen atoms

    #### Complexity
    - **BertzCT**: Bertz complexity index (higher = more complex)
    - **Chi0**: Chi connectivity index 0 (molecular connectivity)
    - **Chi1**: Chi connectivity index 1 (molecular connectivity)

    #### Additional Descriptors
    - **CrippenLogP**: Crippen's LogP calculation method
    - **CrippenMR**: Crippen's molar refractivity
    - **LabuteASA**: Labute's approximate surface area

    ### Usage Tips
    - No properties are selected by default - choose what you need
    - Perfect for files with existing calculations - add only missing properties
    - Common SMILES column names are automatically detected
    - Invalid molecules will result in empty property values
    - All calculations use RDKit library
    - Results are compatible with StarDrop, Pipeline Pilot, and other software
    """)