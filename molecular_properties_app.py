import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from io import StringIO
from molecular_calculator import MolecularCalculator, PropertyExplanations


# Streamlit App
st.set_page_config(
    page_title="ITR - Molecular Properties Calculator",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("🧪 ITR - Molecular Properties Calculator")
st.markdown("Calculate various chemical properties from molecular structures")

# Developer signature
st.markdown("---")
st.markdown("**Developed by:** Yashwanth Reddy for ITR-UIC | **Part of:** Chemo-Informatics Toolkit as part of Dr. Guido Pauli's Team")
st.markdown("---")

# Settings
st.sidebar.subheader("Settings")
suppress_warnings = st.sidebar.checkbox("Suppress RDKit warnings", value=True, help="Hide stereochemistry conflict warnings")
enable_online_lookup = st.sidebar.checkbox("Enable InChI Key conversion", value=True, help="Convert InChI Keys using online databases (NIH CIR, PubChem)")

MolecularCalculator.suppress_rdkit_warnings(suppress_warnings)

# Sidebar for input options
st.sidebar.header("Input Options")
input_mode = st.sidebar.radio("Select input mode:", ["Single Molecule", "Batch Processing"])

if input_mode == "Single Molecule":
    st.header("Single Molecule Analysis")

    # Input molecule
    molecule_input = st.text_area("Enter molecular structure:", placeholder="SMILES, InChI, or InChI Key", help="Paste your molecular structure here and click 'Analyze Molecule' to proceed")

    # Add analyze button for better UX
    col1, col2 = st.columns([3, 1])
    with col1:
        analyze_button = st.button("🔬 Analyze Molecule", type="primary", help="Click to process the molecular structure", disabled=not molecule_input)
    with col2:
        if st.session_state.get('molecule_analyzed', False):
            if st.button("🔄 Reset", help="Clear analysis and start over"):
                st.session_state.molecule_analyzed = False
                st.session_state.current_smiles = None
                st.rerun()

    # Initialize session state for molecule analysis
    if 'molecule_analyzed' not in st.session_state:
        st.session_state.molecule_analyzed = False
    if 'current_smiles' not in st.session_state:
        st.session_state.current_smiles = None

    if molecule_input and analyze_button:
        # Auto-detect format
        detected_format = MolecularCalculator.detect_input_format(molecule_input.strip())
        st.info(f"Detected format: {detected_format.upper()}")

        # Manual format selection
        input_format = st.selectbox("Confirm input format:",
                                   ["smiles", "inchi", "inchi_key"],
                                   index=0 if detected_format == "smiles" else 1 if detected_format == "inchi" else 2)

        # Convert to SMILES if needed
        if input_format == 'inchi_key' and not enable_online_lookup:
            st.warning("InChI Key conversion is disabled. Please enable 'InChI Key conversion' in Settings.")
            smiles = None
        else:
            with st.spinner("Converting InChI Key using online databases..." if input_format == 'inchi_key' else None):
                smiles = MolecularCalculator.convert_to_smiles(molecule_input.strip(), input_format, enable_online_lookup)

        if smiles is None and input_format == 'inchi_key' and enable_online_lookup:
            st.error("❌ Could not resolve this InChI Key. Please verify the format or check your internet connection.")

        if smiles:
            st.success(f"✅ SMILES: {smiles}")
            # Store in session state
            st.session_state.molecule_analyzed = True
            st.session_state.current_smiles = smiles
        else:
            st.session_state.molecule_analyzed = False
            st.session_state.current_smiles = None

    elif not molecule_input and analyze_button:
        st.warning("Please enter a molecular structure first.")

    # Show property selection if we have analyzed a molecule
    if st.session_state.molecule_analyzed and st.session_state.current_smiles:
        # Show current molecule
        #st.success(f"✅ Analyzed molecule: {st.session_state.current_smiles}")

        # Property selection
        st.subheader("Select Properties to Calculate")

        # Get property groups from backend
        property_groups = MolecularCalculator.get_property_groups()

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
                properties = MolecularCalculator.calculate_molecular_properties(st.session_state.current_smiles)

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
                        file_name=f"molecular_properties_{st.session_state.current_smiles.replace('/', '_')[:20]}.csv",
                        mime="text/csv"
                    )
                else:
                    st.error("Could not calculate properties for this molecule.")
            else:
                st.warning("Please select at least one property to calculate.")

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
        smiles_col = MolecularCalculator.detect_smiles_column(df)

        if smiles_col:
            st.success(f"Detected SMILES column: '{smiles_col}'")
        else:
            st.warning("Could not detect SMILES column automatically")
            smiles_col = st.selectbox("Select SMILES column:", df.columns)

        # Property selection for batch
        st.subheader("Select Properties to Calculate")

        # Use same property groups as single molecule
        property_groups_batch = MolecularCalculator.get_property_groups()

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
                    input_format = MolecularCalculator.detect_input_format(str(smiles)) if pd.notna(smiles) else 'smiles'
                    if input_format != 'smiles':
                        smiles = MolecularCalculator.convert_to_smiles(str(smiles), input_format, enable_online_lookup)

                    properties = MolecularCalculator.calculate_molecular_properties(smiles) if pd.notna(smiles) else {}

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
    st.markdown(PropertyExplanations.get_explanations())