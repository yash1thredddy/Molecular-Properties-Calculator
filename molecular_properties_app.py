import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from io import StringIO
from molecular_calculator import MolecularCalculator, PropertyExplanations, ThreeDOLSRegression
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# Import Ligand Efficiency functionality
from ligand_efficiency import (
    DependencyChecker,
    LigandEfficiencyCalculator,
    get_lei_descriptions
)


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
st.markdown("**Developed by:** Yashwanth Reddy for ITR-UIC | **Part of:** Chemo-Informatics Toolkit")
st.markdown("---")

# Settings
st.sidebar.subheader("Settings")
suppress_warnings = st.sidebar.checkbox("Suppress RDKit warnings", value=True, help="Hide stereochemistry conflict warnings")
enable_online_lookup = st.sidebar.checkbox("Enable InChI Key conversion", value=True, help="Convert InChI Keys using online databases (NIH CIR, PubChem)")

MolecularCalculator.suppress_rdkit_warnings(suppress_warnings)

# Sidebar for input options
st.sidebar.header("Input Options")
input_mode = st.sidebar.radio("Select input mode:", ["Single Molecule", "Batch Processing", "Data Visualization", "3D Regression Analysis"])

if input_mode == "Single Molecule":
    st.header("Single Molecule Analysis")

    # Input molecule
    molecule_input = st.text_area("Enter molecular structure:", placeholder="SMILES, InChI, or InChI Key", help="Paste your molecular structure here - analysis will happen automatically")

    # Add reset button only if molecule is analyzed
    if st.session_state.get('molecule_analyzed', False):
        col1, col2 = st.columns([3, 1])
        with col2:
            if st.button("🔄 Reset", help="Clear analysis and start over"):
                st.session_state.molecule_analyzed = False
                st.session_state.current_smiles = None
                st.rerun()

    # Initialize session state for molecule analysis
    if 'molecule_analyzed' not in st.session_state:
        st.session_state.molecule_analyzed = False
    if 'current_smiles' not in st.session_state:
        st.session_state.current_smiles = None

    # Auto-analyze when input changes
    if molecule_input and molecule_input.strip():
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

    elif molecule_input and not molecule_input.strip():
        st.warning("Please enter a valid molecular structure.")

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
                            if prop == 'QED':
                                display_name = 'QED'
                            else:
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

elif input_mode == "Batch Processing":
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

        # Detect SMILES column first
        smiles_col = MolecularCalculator.detect_smiles_column(df)

        # Data preprocessing for batch mode (avoid converting SMILES column)
        original_dtypes_batch = df.dtypes.to_dict()
        converted_columns_batch = []

        # Try to convert string columns that might contain numbers
        for col in df.columns:
            if df[col].dtype == 'object' and col != smiles_col:  # Don't convert SMILES column
                numeric_converted = pd.to_numeric(df[col], errors='coerce')
                non_null_original = df[col].notna().sum()
                non_null_converted = numeric_converted.notna().sum()

                if non_null_converted > 0 and (non_null_converted / non_null_original) >= 0.5:
                    df[col] = numeric_converted
                    converted_columns_batch.append(col)

        if converted_columns_batch:
            st.info(f"✅ Auto-converted {len(converted_columns_batch)} columns to numeric: {', '.join(converted_columns_batch)}")

        # Clear previous results when new file is uploaded (DO THIS FIRST!)
        if 'current_batch_file' not in st.session_state:
            st.session_state.current_batch_file = None

        # Check if this is a new file
        if st.session_state.current_batch_file != uploaded_file.name:
            # New file detected - clear all previous selections
            st.session_state.current_batch_file = uploaded_file.name
            st.session_state.batch_results_df = None
            st.session_state.batch_final_df = None
            st.session_state.batch_selected_properties = set()  # Clear property selections
            st.session_state.batch_selected_leis = set()  # Clear LEI selections
            st.session_state.column_mappings = {}  # Clear column mappings
            st.session_state.mapping_applied = False  # Reset mapping status
            st.session_state.file_just_changed = True  # Flag to trigger rerun

        # If file just changed, show message and trigger rerun immediately
        if st.session_state.get('file_just_changed', False):
            st.session_state.file_just_changed = False  # Reset flag
            st.info("🔄 New file detected - clearing previous selections...")
            st.rerun()

        st.subheader("Data Preview")
        st.info(f"Loaded {len(df)} rows from {uploaded_file.name}")
        st.dataframe(df.head(), use_container_width=True)

        # ==================== COLUMN MAPPING SECTION ====================
        st.markdown("---")
        st.subheader("🔗 Column Mapping")

        # Initialize column mapping state
        if 'column_mappings' not in st.session_state:
            st.session_state.column_mappings = {}

        # Check for manual mapping toggle
        show_manual_mapping = st.checkbox(
            "📝 Manually map columns to standard names",
            value=False,
            help="Use this if your columns have non-standard names (e.g., 'molecular' instead of 'SMILES', 'MICMtb(μM)' instead of 'pKi')"
        )

        if show_manual_mapping:
            with st.expander("🗺️ Column Mapping Configuration", expanded=True):
                st.markdown("**Map your file columns to standard property names:**")

                # Create mapping UI
                map_col1, map_col2 = st.columns(2)

                with map_col1:
                    st.markdown("**Essential Columns:**")

                    # SMILES mapping
                    smiles_mapping = st.selectbox(
                        "SMILES column:",
                        options=['[Not mapped]'] + list(df.columns),
                        index=list(df.columns).index(smiles_col) + 1 if smiles_col else 0,
                        key="map_smiles",
                        help="Select the column containing molecular structures (SMILES/InChI)"
                    )
                    if smiles_mapping != '[Not mapped]':
                        st.session_state.column_mappings['smiles'] = smiles_mapping
                        smiles_col = smiles_mapping  # Update smiles_col
                    else:
                        smiles_col = None

                    # pKi mapping
                    pki_detected_auto = DependencyChecker.detect_column(df, 'pki')
                    pki_mapping = st.selectbox(
                        "pKi column:",
                        options=['[Not mapped]'] + list(df.columns),
                        index=list(df.columns).index(pki_detected_auto) + 1 if pki_detected_auto else 0,
                        key="map_pki",
                        help="Select the column containing pKi values (required for LEI calculations)"
                    )
                    if pki_mapping != '[Not mapped]':
                        st.session_state.column_mappings['pki'] = pki_mapping

                with map_col2:
                    st.markdown("**Optional Property Columns:**")
                    st.markdown("*These will be calculated from SMILES if not mapped*")

                    # MW mapping
                    mw_detected_auto = DependencyChecker.detect_column(df, 'mw')
                    mw_mapping = st.selectbox(
                        "Molecular Weight (MW):",
                        options=['[Not mapped]'] + list(df.columns),
                        index=list(df.columns).index(mw_detected_auto) + 1 if mw_detected_auto else 0,
                        key="map_mw"
                    )
                    if mw_mapping != '[Not mapped]':
                        st.session_state.column_mappings['mw'] = mw_mapping

                    # TPSA mapping
                    tpsa_detected_auto = DependencyChecker.detect_column(df, 'tpsa')
                    tpsa_mapping = st.selectbox(
                        "TPSA:",
                        options=['[Not mapped]'] + list(df.columns),
                        index=list(df.columns).index(tpsa_detected_auto) + 1 if tpsa_detected_auto else 0,
                        key="map_tpsa"
                    )
                    if tpsa_mapping != '[Not mapped]':
                        st.session_state.column_mappings['tpsa'] = tpsa_mapping

                    # Heavy Atoms mapping
                    heavy_detected_auto = DependencyChecker.detect_column(df, 'heavy_atoms')
                    heavy_mapping = st.selectbox(
                        "Heavy Atom Count:",
                        options=['[Not mapped]'] + list(df.columns),
                        index=list(df.columns).index(heavy_detected_auto) + 1 if heavy_detected_auto else 0,
                        key="map_heavy"
                    )
                    if heavy_mapping != '[Not mapped]':
                        st.session_state.column_mappings['heavy_atoms'] = heavy_mapping

                # Show current mappings
                if st.session_state.column_mappings:
                    st.success("✅ **Current Mappings:**")
                    for std_name, col_name in st.session_state.column_mappings.items():
                        st.write(f"  • `{col_name}` → **{std_name.upper()}**")

                # Apply Mapping Button
                st.markdown("---")
                if st.button("✅ Apply Mapping", type="primary", key="apply_mapping_btn"):
                    st.session_state.mapping_applied = True
                    st.success("🎉 Column mapping applied! You can now proceed with property and LEI selection.")
                    st.rerun()

            # Initialize mapping_applied state
            if 'mapping_applied' not in st.session_state:
                st.session_state.mapping_applied = False

            # Show status if mapping not applied yet
            if show_manual_mapping and not st.session_state.mapping_applied:
                st.warning("⚠️ Please click 'Apply Mapping' button to confirm your column mappings before proceeding.")
        else:
            # Auto-detect SMILES column (no manual mapping)
            st.session_state.mapping_applied = True  # Auto-detection is instant
            if smiles_col:
                st.success(f"✅ Auto-detected SMILES column: '{smiles_col}'")
            else:
                st.warning("⚠️ Could not detect SMILES column automatically")
                st.info("💡 Please select your SMILES/InChI/InChIKey column from the dropdown, or use 'Manually map columns' option above")
                selected_smiles = st.selectbox(
                    "Select SMILES column:",
                    options=['[Not selected]'] + list(df.columns),
                    key="manual_smiles_select",
                    help="Choose the column containing molecular structures (SMILES/InChI/InChIKey)"
                )
                if selected_smiles != '[Not selected]':
                    smiles_col = selected_smiles
                    st.success(f"✅ Selected SMILES column: '{smiles_col}'")
                else:
                    smiles_col = None
        # ==================== END COLUMN MAPPING SECTION ====================

        # Property selection for batch
        st.subheader("Select Properties to Calculate")

        # Use same property groups as single molecule
        property_groups_batch = MolecularCalculator.get_property_groups()

        # Initialize batch session state
        if 'batch_selected_properties' not in st.session_state:
            st.session_state.batch_selected_properties = set()

        # Calculate All checkbox for batch
        # Track previous state to detect when user unchecks "Calculate All"
        if 'prev_calc_all_batch' not in st.session_state:
            st.session_state.prev_calc_all_batch = False

        calc_all_batch = st.checkbox("Calculate All Properties (Batch)", value=False)

        if calc_all_batch:
            # Select all properties when checked
            st.session_state.batch_selected_properties = set()
            for props in property_groups_batch.values():
                st.session_state.batch_selected_properties.update(props)
            # Update tracking flag
            st.session_state.prev_calc_all_batch = True
        else:
            # If previously was "calc all" and now unchecked, clear all selections
            if st.session_state.prev_calc_all_batch and not calc_all_batch:
                # User just unchecked "Calculate All" - clear selections
                st.session_state.batch_selected_properties = set()

            # Update the previous state
            st.session_state.prev_calc_all_batch = False

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
                            if prop == 'QED':
                                display_name = 'QED'
                            else:
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

        # ==================== LEI CALCULATIONS SECTION ====================
        # Initialize LEI session state FIRST (before anything else)
        if 'batch_selected_leis' not in st.session_state:
            st.session_state.batch_selected_leis = set()

        # Only show LEI section if mapping is applied
        if st.session_state.get('mapping_applied', False):
            st.markdown("---")
            st.subheader("⚗️ Ligand Efficiency Indices (LEI)")
            st.markdown("*Based on AtlasCBS methodology (Cele Abad-Zapatero, A. Cortes-Cabrera 2013) - requires pKi column*")

            # Check if pKi column exists (use manual mapping if available)
            if 'pki' in st.session_state.column_mappings:
                pki_detected = st.session_state.column_mappings['pki']
            else:
                pki_detected = DependencyChecker.detect_column(df, 'pki')

            if pki_detected:
                st.success(f"✅ pKi column detected: '{pki_detected}'")

                # Override detected_cols with manual mappings if they exist
                detected_cols = DependencyChecker.detect_all_columns(df)
                # Use manual mappings to override auto-detection
                for std_name, col_name in st.session_state.column_mappings.items():
                    detected_cols[std_name] = col_name

                # Check if both SMILES and pKi are available (required for LEI calculations)
                has_smiles = detected_cols.get('smiles') is not None
                has_pki = pki_detected is not None
                can_calculate_all_leis = has_smiles and has_pki

                # Show detected columns status
                with st.expander("📊 Column Detection Status", expanded=True):
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

                # Check which LEIs can be calculated (pass manual mappings)
                lei_check = DependencyChecker.check_lei_dependencies(
                    df, all_leis, manual_mappings=st.session_state.column_mappings
                )

                # If both SMILES and pKi are available, all LEIs can be calculated
                if can_calculate_all_leis:
                    # Enable all LEIs since we have SMILES (to calculate missing properties) and pKi
                    can_calculate_leis = set(all_leis)
                    st.info("✨ SMILES and pKi detected - all LEIs enabled (missing properties will be calculated from SMILES)")
                else:
                    # Use normal dependency checking
                    can_calculate_leis = set(lei_check.get('can_calculate', []))

                # Create three columns for LEI selection
                lei_col1, lei_col2, lei_col3 = st.columns(3)

                for idx, lei in enumerate(all_leis):
                    col = [lei_col1, lei_col2, lei_col3][idx % 3]

                    with col:
                        # Check if this LEI can be calculated
                        can_calculate = lei in can_calculate_leis
                        status_icon = "✅" if can_calculate else "⚠️"

                    lei_selected = lei in st.session_state.batch_selected_leis
                    checkbox_label = f"{status_icon} {lei}"

                    if st.checkbox(checkbox_label, value=lei_selected, key=f"lei_{lei}",
                                 help=lei_descriptions[lei], disabled=not can_calculate):
                        st.session_state.batch_selected_leis.add(lei)

                        # Auto-select required molecular properties for this LEI
                        if lei in lei_check['status_by_lei']:
                            needs_calc = lei_check['status_by_lei'][lei].get('needs_calc', [])

                            # Auto-select MW, TPSA, Heavy_Atom_Count if needed and not present in file
                            for dep in needs_calc:
                                if dep == 'heavy_atoms' and 'Heavy_Atom_Count' not in [detected_cols.get('heavy_atoms')]:
                                    st.session_state.batch_selected_properties.add('Heavy_Atom_Count')
                                elif dep == 'mw' and detected_cols.get('mw') is None:
                                    st.session_state.batch_selected_properties.add('Molecular_Weight')
                                elif dep == 'tpsa' and detected_cols.get('tpsa') is None:
                                    st.session_state.batch_selected_properties.add('TPSA')
                    else:
                        st.session_state.batch_selected_leis.discard(lei)

                # Show what will be calculated
                if st.session_state.batch_selected_leis:
                    selected_lei_check = DependencyChecker.check_lei_dependencies(
                        df, list(st.session_state.batch_selected_leis),
                        manual_mappings=st.session_state.column_mappings
                    )

                    # Check if any properties were auto-selected
                    auto_selected = []
                    needs_calc = selected_lei_check.get('needs_calculation', [])
                    if 'heavy_atoms' in needs_calc or 'polar_atoms' in needs_calc:
                        if 'Heavy_Atom_Count' in st.session_state.batch_selected_properties:
                            auto_selected.append('Heavy_Atom_Count')
                    if 'mw' in needs_calc and detected_cols.get('mw') is None:
                        if 'Molecular_Weight' in st.session_state.batch_selected_properties:
                            auto_selected.append('Molecular_Weight')
                    if 'tpsa' in needs_calc and detected_cols.get('tpsa') is None:
                        if 'TPSA' in st.session_state.batch_selected_properties:
                            auto_selected.append('TPSA')

                    if auto_selected:
                        st.success(f"✨ Auto-selected required properties: {', '.join(auto_selected)}")

                    with st.expander("ℹ️ Calculation Plan", expanded=False):
                        status_msg = DependencyChecker.generate_status_message(selected_lei_check)
                        st.markdown(status_msg)

                    st.info(f"📋 Selected LEIs: {len(st.session_state.batch_selected_leis)}")

            else:
                st.info("ℹ️ **No pKi column detected** - LEI calculations require a pKi column")
                st.markdown("**Supported pKi column names:**")
                st.code("pKi, pki, PKI, pIC50, pic50, p_Ki, p_ki")
                st.markdown("Upload a file with a pKi column to enable LEI calculations.")
        # ==================== END LEI SECTION ====================

        # Initialize batch results session state
        if 'batch_results_df' not in st.session_state:
            st.session_state.batch_results_df = None
        if 'batch_final_df' not in st.session_state:
            st.session_state.batch_final_df = None

        # ==================== VALIDATION BEFORE PROCESSING ====================
        st.markdown("---")
        st.subheader("🚀 Ready to Process?")

        # Check if SMILES column is detected
        if smiles_col:
            st.success(f"✅ **Molecular structure column detected:** `{smiles_col}`")
        else:
            st.error("❌ **No molecular structure column detected**")
            st.warning("⚠️ **Cannot proceed without SMILES, InChI, or InChIKey column**")
            st.info("💡 **How to fix this:**\n" +
                   "- Make sure your file has a column with molecular structures (SMILES/InChI/InChIKey)\n" +
                   "- Supported column names: `SMILES`, `smiles`, `InChI`, `inchi`, `InChIKey`, `inchikey`\n" +
                   "- OR use 'Manually map columns' option above to map your custom column name")

        # Check if any properties or LEIs are selected
        has_selections = calc_all_batch or st.session_state.batch_selected_properties or st.session_state.batch_selected_leis
        if not has_selections and smiles_col:
            st.warning("⚠️ **No properties or LEIs selected** - Please select at least one property or LEI to calculate")
        elif has_selections and smiles_col:
            selection_count = 0
            if calc_all_batch:
                selection_count = sum(len(props) for props in property_groups_batch.values())
                st.info(f"📊 **Ready to calculate:** All properties ({selection_count} total)")
            else:
                if st.session_state.batch_selected_properties:
                    selection_count += len(st.session_state.batch_selected_properties)
                if st.session_state.batch_selected_leis:
                    selection_count += len(st.session_state.batch_selected_leis)
                items = []
                if st.session_state.batch_selected_properties:
                    items.append(f"{len(st.session_state.batch_selected_properties)} properties")
                if st.session_state.batch_selected_leis:
                    items.append(f"{len(st.session_state.batch_selected_leis)} LEIs")
                st.info(f"📊 **Ready to calculate:** {' + '.join(items)}")

        # Show the Process Batch button with appropriate state
        can_process = smiles_col and has_selections

        if st.button("Process Batch", type="primary", disabled=not can_process):
            if smiles_col and (calc_all_batch or st.session_state.batch_selected_properties or st.session_state.batch_selected_leis):
                progress_bar = st.progress(0)
                status_text = st.empty()

                results = []
                total_rows = len(df)

                # Calculate regular molecular properties
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
                    progress_bar.progress(progress * 0.7)  # First 70% for regular properties
                    status_text.text(f"Calculating properties: {idx + 1}/{total_rows}")

                # Create results DataFrame
                results_df = pd.DataFrame(results)
                final_df = pd.concat([df, results_df], axis=1)

                # Calculate LEIs if selected
                if st.session_state.batch_selected_leis:
                    status_text.text("Calculating Ligand Efficiency Indices...")
                    progress_bar.progress(0.75)

                    lei_result_df, lei_status = LigandEfficiencyCalculator.process_batch(
                        df=final_df,
                        selected_leis=list(st.session_state.batch_selected_leis),
                        show_errors=False,
                        manual_mappings=st.session_state.get('column_mappings', {})
                    )

                    if lei_status['success']:
                        final_df = lei_result_df
                        status_text.text(f"LEI calculation complete! Calculated: {', '.join(lei_status['calculated_leis'])}")
                    else:
                        st.warning(f"⚠️ LEI calculation issue: {lei_status.get('message', 'Unknown error')}")

                progress_bar.progress(1.0)

                # Store results in session state
                st.session_state.batch_results_df = results_df
                st.session_state.batch_final_df = final_df

                status_text.text("✅ Processing complete!")
            else:
                st.error("Please select a SMILES column and at least one property or LEI.")

        # Show helpful message when no processing has been done yet
        if (st.session_state.batch_results_df is None and smiles_col and
            (calc_all_batch or st.session_state.batch_selected_properties)):
            st.info("👆 Click 'Process Batch' above to calculate properties for your molecules.")

        # Display results if they exist in session state and are for the current file
        if (st.session_state.batch_results_df is not None and
            st.session_state.batch_final_df is not None and
            st.session_state.current_batch_file == uploaded_file.name):
            st.subheader("Results")
            st.dataframe(st.session_state.batch_final_df, use_container_width=True)

            # Summary statistics
            if not st.session_state.batch_results_df.empty:
                st.subheader("Summary Statistics")
                summary = st.session_state.batch_results_df.describe()
                st.dataframe(summary, use_container_width=True)

                # Visualizations
                st.subheader("📊 Property Visualizations")

                # Create plots for key properties if they exist
                plot_properties = ['QED', 'LogP', 'Molecular_Weight', 'TPSA']
                available_props = [prop for prop in plot_properties if prop in st.session_state.batch_results_df.columns]

                if available_props:
                    st.subheader("🔍 Distribution Analysis")
                    cols = st.columns(2)

                    # Color palette for histograms
                    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4']

                    for i, prop in enumerate(available_props[:4]):  # Limit to 4 plots max
                        with cols[i % 2]:
                            # Create enhanced histogram with reference lines
                            fig = px.histogram(
                                st.session_state.batch_results_df,
                                x=prop,
                                title=f"{prop.replace('_', ' ')} Distribution",
                                nbins=20,
                                color_discrete_sequence=[colors[i % len(colors)]]
                            )

                            # Add mean line
                            mean_val = st.session_state.batch_results_df[prop].mean()
                            fig.add_vline(
                                x=mean_val,
                                line_dash="dash",
                                line_color="red",
                                annotation_text=f"Mean: {mean_val:.2f}",
                                annotation_position="top"
                            )

                            # Add reference lines for drug-like properties
                            if prop == 'QED':
                                fig.add_vline(x=0.5, line_dash="dot", line_color="orange",
                                            annotation_text="QED > 0.5 (Drug-like)", annotation_position="bottom")
                            elif prop == 'LogP':
                                fig.add_vline(x=5, line_dash="dot", line_color="orange",
                                            annotation_text="LogP < 5 (Lipinski)", annotation_position="bottom")
                            elif prop == 'Molecular_Weight':
                                fig.add_vline(x=500, line_dash="dot", line_color="orange",
                                            annotation_text="MW < 500 (Lipinski)", annotation_position="bottom")
                            elif prop == 'TPSA':
                                fig.add_vline(x=140, line_dash="dot", line_color="orange",
                                            annotation_text="TPSA < 140 (Drug-like)", annotation_position="bottom")

                            fig.update_layout(
                                xaxis_title=prop.replace('_', ' '),
                                yaxis_title="Count",
                                height=400,
                                showlegend=False,
                                plot_bgcolor='rgba(0,0,0,0)',
                                paper_bgcolor='rgba(0,0,0,0)'
                            )
                            st.plotly_chart(fig, use_container_width=True)

                # Enhanced Interactive Visualization System
                # Include ALL numeric columns from final_df (original + calculated)
                numeric_cols = st.session_state.batch_final_df.select_dtypes(include=[np.number]).columns.tolist()
                all_cols = st.session_state.batch_final_df.columns.tolist()

                if len(numeric_cols) >= 1:
                    st.subheader("🎨 Interactive Data Visualization")

                    # Chart type selection
                    chart_types = {
                        'Scatter Plot': {'requires': ['x', 'y'], 'desc': 'Compare two numeric variables'},
                        '3D OLS Regression': {'requires': ['x', 'y', 'z'], 'desc': 'Fit 3D plane: Z = b0 + b1·X + b2·Y'},
                        'Histogram': {'requires': ['x'], 'desc': 'Distribution of a single variable'},
                        'Box Plot': {'requires': ['y'], 'desc': 'Statistical summary with quartiles'},
                        'Violin Plot': {'requires': ['y'], 'desc': 'Distribution shape with density'},
                        'Line Plot': {'requires': ['x', 'y'], 'desc': 'Trends over continuous data'},
                        'Heatmap': {'requires': [], 'desc': 'Correlation matrix of all numeric variables'}
                    }

                    col1, col2 = st.columns([1, 2])
                    with col1:
                        selected_chart = st.selectbox(
                            "📊 Chart Type:",
                            options=list(chart_types.keys()),
                            key="chart_type_select"
                        )
                    with col2:
                        st.info(f"💡 {chart_types[selected_chart]['desc']}")

                    # Dynamic input controls based on chart type
                    requires = chart_types[selected_chart]['requires']

                    # Create input columns
                    input_cols = st.columns(4)

                    # Initialize variables
                    x_axis, y_axis, z_axis, color_col, size_col = None, None, None, None, None

                    with input_cols[0]:
                        if 'x' in requires:
                            x_axis = st.selectbox(
                                "🔢 X-axis:",
                                options=numeric_cols,
                                index=numeric_cols.index('Molecular_Weight') if 'Molecular_Weight' in numeric_cols else 0,
                                key="viz_x_axis"
                            )

                    with input_cols[1]:
                        if 'y' in requires:
                            available_y = [col for col in numeric_cols if col != x_axis] if x_axis else numeric_cols
                            y_axis = st.selectbox(
                                "📈 Y-axis:",
                                options=available_y,
                                index=available_y.index('LogP') if 'LogP' in available_y else 0,
                                key="viz_y_axis"
                            )
                        elif selected_chart in ['Histogram', 'Box Plot', 'Violin Plot']:
                            y_axis = st.selectbox(
                                "📊 Variable:",
                                options=numeric_cols,
                                index=numeric_cols.index('QED') if 'QED' in numeric_cols else 0,
                                key="viz_single_var"
                            )

                    # Z-axis for 3D OLS Regression
                    if 'z' in requires:
                        available_z = [col for col in numeric_cols if col != x_axis and col != y_axis]
                        if available_z:
                            z_axis = st.selectbox(
                                "📊 Z-axis (dependent):",
                                options=available_z,
                                index=available_z.index('TPSA') if 'TPSA' in available_z else 0,
                                key="viz_z_axis"
                            )
                        else:
                            st.warning("Need at least 3 numeric columns for 3D regression")

                    with input_cols[2]:
                        # Color mapping options
                        color_options = ['None'] + all_cols
                        color_col = st.selectbox(
                            "🎨 Color by:",
                            options=color_options,
                            key="viz_color_col"
                        )
                        if color_col == 'None':
                            color_col = None

                    with input_cols[3]:
                        # Additional options based on chart type
                        if selected_chart == 'Scatter Plot':
                            show_trendline = st.checkbox("📈 Trend Line", value=False, key="viz_trendline")
                        elif selected_chart == 'Histogram':
                            bins = st.slider("📊 Bins:", min_value=5, max_value=50, value=20, key="viz_bins")
                        elif selected_chart in ['Box Plot', 'Violin Plot']:
                            show_points = st.checkbox("🔹 Show Points", value=True, key="viz_points")

                    # Color scale options (if color column is selected)
                    color_scale = 'viridis'  # Default
                    if color_col and color_col in numeric_cols:
                        st.markdown("**🌈 Color Scale Options:**")
                        color_scale_options = {
                            'Viridis': 'viridis',
                            'Plasma': 'plasma',
                            'Inferno': 'inferno',
                            'Rainbow': 'rainbow',
                            'RdYlBu': 'RdYlBu',
                            'Blue-Red': 'RdBu',
                            'Spectral': 'Spectral'
                        }
                        color_scale_col1, color_scale_col2 = st.columns(2)
                        with color_scale_col1:
                            color_scale = st.selectbox(
                                "Color Scale:",
                                options=list(color_scale_options.values()),
                                format_func=lambda x: [k for k, v in color_scale_options.items() if v == x][0],
                                index=0,
                                key="viz_color_scale"
                            )
                        with color_scale_col2:
                            reverse_scale = st.checkbox("🔄 Reverse Scale", value=False, key="viz_reverse_scale")
                            if reverse_scale:
                                color_scale = color_scale + '_r'

                    # Size parameter for scatter plots
                    size_col = None
                    if selected_chart == 'Scatter Plot':
                        st.markdown("**📏 Size Options:**")
                        size_options = ['None'] + numeric_cols
                        size_col = st.selectbox(
                            "Size by:",
                            options=size_options,
                            key="viz_size_col"
                        )
                        if size_col == 'None':
                            size_col = None

                    # Faceting/Grouping options removed for clarity
                    facet_col = None

                    # Generate visualization automatically when parameters change
                    if True:  # Always generate visualization
                        try:
                            # Prepare data
                            plot_data = st.session_state.batch_final_df.copy()

                            # Faceting validation removed

                            # Handle color and size columns
                            color_param = {}
                            if color_col:
                                if color_col in numeric_cols:
                                    # Continuous color scale
                                    color_param = {
                                        'color': color_col,
                                        'color_continuous_scale': color_scale
                                    }
                                else:
                                    # Discrete color scale
                                    color_param = {'color': color_col}

                            # Add size parameter for scatter plots
                            if size_col and selected_chart == 'Scatter Plot':
                                color_param['size'] = size_col

                            # Faceting parameter removed

                            # Create the selected chart
                            if selected_chart == 'Scatter Plot' and x_axis and y_axis:
                                fig = px.scatter(
                                    plot_data,
                                    x=x_axis,
                                    y=y_axis,
                                    title=f"{x_axis.replace('_', ' ')} vs {y_axis.replace('_', ' ')}",
                                    trendline="ols" if show_trendline else None,
                                    **color_param
                                )

                                if not color_col:
                                    # Add correlation coefficient and regression stats
                                    correlation = plot_data[x_axis].corr(plot_data[y_axis])

                                    # Calculate regression statistics if trendline is shown
                                    if show_trendline:

                                        # Prepare data for regression
                                        X = plot_data[[x_axis]].dropna()
                                        y = plot_data[y_axis].dropna()

                                        # Ensure same indices
                                        common_idx = X.index.intersection(y.index)
                                        X = X.loc[common_idx]
                                        y = y.loc[common_idx]

                                        if len(X) > 1:
                                            # Fit linear regression
                                            reg = LinearRegression()
                                            reg.fit(X, y)

                                            # Calculate R²
                                            y_pred = reg.predict(X)
                                            r2 = r2_score(y, y_pred)

                                            # Get slope and intercept
                                            slope = reg.coef_[0]
                                            intercept = reg.intercept_

                                            # Format equation with actual variable names
                                            sign = "+" if intercept >= 0 else "-"
                                            y_name = y_axis.replace('_', ' ')
                                            x_name = x_axis.replace('_', ' ')
                                            equation = f"{y_name} = {slope:.3f} × {x_name} {sign} {abs(intercept):.3f}"

                                            # Add regression info
                                            fig.add_annotation(
                                                x=0.02, y=0.98,
                                                xref="paper", yref="paper",
                                                text=f"<b>Regression Statistics:</b><br>" +
                                                     f"Equation: {equation}<br>" +
                                                     f"R² = {r2:.3f}<br>" +
                                                     f"Correlation: {correlation:.3f}",
                                                showarrow=False,
                                                bgcolor="rgba(255,255,255,0.9)",
                                                bordercolor="black",
                                                borderwidth=1,
                                                align="left",
                                                font=dict(size=11)
                                            )
                                    else:
                                        # Just show correlation if no trendline
                                        fig.add_annotation(
                                            x=0.02, y=0.98,
                                            xref="paper", yref="paper",
                                            text=f"Correlation: {correlation:.3f}",
                                            showarrow=False,
                                            bgcolor="rgba(255,255,255,0.8)",
                                            bordercolor="black",
                                            borderwidth=1
                                        )

                            elif selected_chart == 'Histogram' and y_axis:
                                fig = px.histogram(
                                    plot_data,
                                    x=y_axis,
                                    title=f"Distribution of {y_axis.replace('_', ' ')}",
                                    nbins=bins,
                                    **color_param
                                )

                                # Add mean and median lines
                                mean_val = plot_data[y_axis].mean()
                                median_val = plot_data[y_axis].median()
                                fig.add_vline(x=mean_val, line_dash="dash", line_color="red",
                                            annotation_text=f"Mean: {mean_val:.2f}")
                                fig.add_vline(x=median_val, line_dash="dot", line_color="blue",
                                            annotation_text=f"Median: {median_val:.2f}")

                            elif selected_chart == 'Box Plot' and y_axis:
                                fig = px.box(
                                    plot_data,
                                    y=y_axis,
                                    title=f"Box Plot of {y_axis.replace('_', ' ')}",
                                    points="all" if show_points else False,
                                    **color_param
                                )

                            elif selected_chart == 'Violin Plot' and y_axis:
                                fig = px.violin(
                                    plot_data,
                                    y=y_axis,
                                    title=f"Violin Plot of {y_axis.replace('_', ' ')}",
                                    points="all" if show_points else False,
                                    **color_param
                                )

                            elif selected_chart == 'Line Plot' and x_axis and y_axis:
                                # Sort data for line plot
                                plot_data_sorted = plot_data.sort_values(x_axis)
                                fig = px.line(
                                    plot_data_sorted,
                                    x=x_axis,
                                    y=y_axis,
                                    title=f"{y_axis.replace('_', ' ')} over {x_axis.replace('_', ' ')}",
                                    **color_param
                                )

                            elif selected_chart == '3D OLS Regression' and x_axis and y_axis and z_axis:
                                # Perform 3D OLS Regression
                                try:
                                    # Create OLS regression model
                                    ols_model = ThreeDOLSRegression(
                                        x=plot_data[x_axis],
                                        y=plot_data[y_axis],
                                        z=plot_data[z_axis]
                                    )

                                    # Get statistics
                                    stats = ols_model.get_statistics()
                                    equation = ols_model.get_equation_string(decimals=3)

                                    # Get plane mesh for visualization
                                    X_mesh, Y_mesh, Z_mesh = ols_model.get_plane_mesh(num_points=20)

                                    # Create 3D scatter plot
                                    fig = go.Figure()

                                    # Add data points
                                    if color_col and color_col in numeric_cols:
                                        # Color by a variable
                                        fig.add_trace(go.Scatter3d(
                                            x=ols_model.x,
                                            y=ols_model.y,
                                            z=ols_model.z,
                                            mode='markers',
                                            marker=dict(
                                                size=5,
                                                color=plot_data.loc[plot_data.index.isin(plot_data.index), color_col],
                                                colorscale=color_scale,
                                                showscale=True,
                                                colorbar=dict(title=color_col.replace('_', ' '))
                                            ),
                                            name='Data Points',
                                            text=[f"{x_axis}: {x:.2f}<br>{y_axis}: {y:.2f}<br>{z_axis}: {z:.2f}"
                                                  for x, y, z in zip(ols_model.x, ols_model.y, ols_model.z)],
                                            hovertemplate='%{text}<extra></extra>'
                                        ))
                                    else:
                                        # Uniform color, can color by residuals
                                        residuals_abs = np.abs(ols_model.residuals)
                                        fig.add_trace(go.Scatter3d(
                                            x=ols_model.x,
                                            y=ols_model.y,
                                            z=ols_model.z,
                                            mode='markers',
                                            marker=dict(
                                                size=5,
                                                color=residuals_abs,
                                                colorscale='Reds',
                                                showscale=True,
                                                colorbar=dict(title='|Residual|')
                                            ),
                                            name='Data Points',
                                            text=[f"{x_axis}: {x:.2f}<br>{y_axis}: {y:.2f}<br>{z_axis}: {z:.2f}<br>Residual: {r:.3f}"
                                                  for x, y, z, r in zip(ols_model.x, ols_model.y, ols_model.z, ols_model.residuals)],
                                            hovertemplate='%{text}<extra></extra>'
                                        ))

                                    # Compute plane corners and add a lightweight Mesh3d plane (faster than Surface)
                                    x_min_d, x_max_d = float(np.min(ols_model.x)), float(np.max(ols_model.x))
                                    y_min_d, y_max_d = float(np.min(ols_model.y)), float(np.max(ols_model.y))
                                    # 10% padding similar to get_plane_mesh
                                    x_pad = 0.1 * (x_max_d - x_min_d) if x_max_d > x_min_d else 1.0
                                    y_pad = 0.1 * (y_max_d - y_min_d) if y_max_d > y_min_d else 1.0
                                    x0, x1 = x_min_d - x_pad, x_max_d + x_pad
                                    y0, y1 = y_min_d - y_pad, y_max_d + y_pad

                                    corners_x = np.array([x0, x1, x1, x0])
                                    corners_y = np.array([y0, y0, y1, y1])
                                    corners_z = ols_model.predict(corners_x, corners_y)

                                    fig.add_trace(go.Mesh3d(
                                        x=corners_x,
                                        y=corners_y,
                                        z=corners_z,
                                        i=[0, 1, 2, 0, 2, 3],
                                        opacity=0.6,
                                        color='royalblue',
                                        lighting=dict(
                                            ambient=0.5,
                                            diffuse=0.7,
                                            specular=0.2,
                                            roughness=0.9,
                                            fresnel=0.2
                                        ),
                                        lightposition=dict(x=200, y=100, z=300),
                                        name='Fitted Plane',
                                        hoverinfo='skip'
                                    ))

                                    # Update layout for 3D plot
                                    # Fix axis ranges and preserve camera to avoid flipping/jumps
                                    z_min_d, z_max_d = float(np.min(ols_model.z)), float(np.max(ols_model.z))
                                    z_plane_min, z_plane_max = float(np.min(Z_mesh)), float(np.max(Z_mesh))
                                    z0, z1 = min(z_min_d, z_plane_min), max(z_max_d, z_plane_max)

                                    fig.update_layout(
                                        title=f"3D OLS Regression: {z_axis} vs {x_axis} and {y_axis}",
                                        uirevision="3d-ols",  # keep view stable across reruns
                                        paper_bgcolor='#0e1117',
                                        font=dict(color='#d9d9d9'),
                                        scene=dict(
                                            bgcolor='#0e1117',
                                            xaxis_title=x_axis.replace('_', ' '),
                                            yaxis_title=y_axis.replace('_', ' '),
                                            zaxis_title=z_axis.replace('_', ' '),
                                            xaxis=dict(
                                                range=[x0, x1], showspikes=False,
                                                showbackground=True, backgroundcolor='#1b1f2a',
                                                gridcolor='#2b3242', zerolinecolor='#444',
                                                tickfont=dict(color='#c2c7cf'), titlefont=dict(color='#c2c7cf')
                                            ),
                                            yaxis=dict(
                                                range=[y0, y1], showspikes=False,
                                                showbackground=True, backgroundcolor='#1b1f2a',
                                                gridcolor='#2b3242', zerolinecolor='#444',
                                                tickfont=dict(color='#c2c7cf'), titlefont=dict(color='#c2c7cf')
                                            ),
                                            zaxis=dict(
                                                range=[z0, z1], showspikes=False,
                                                showbackground=True, backgroundcolor='#1b1f2a',
                                                gridcolor='#2b3242', zerolinecolor='#444',
                                                tickfont=dict(color='#c2c7cf'), titlefont=dict(color='#c2c7cf')
                                            ),
                                            aspectmode='data',
                                            camera=dict(
                                                eye=dict(x=1.8, y=1.8, z=1.2),
                                                up=dict(x=0, y=0, z=1),
                                                projection=dict(type='perspective')
                                            )
                                        ),
                                        scene_dragmode='orbit',
                                        height=700,
                                        showlegend=True
                                    )

                                    # Display the figure with smoother config
                                    st.plotly_chart(
                                        fig,
                                        use_container_width=True,
                                        config={
                                            'displaylogo': False,
                                            'scrollZoom': True,
                                            'responsive': True
                                        }
                                    )

                                    # Display regression statistics
                                    st.markdown("### 📊 3D OLS Regression Statistics")
                                    col1, col2, col3, col4, col5 = st.columns(5)
                                    with col1:
                                        st.metric("R² (Coefficient of Determination)", f"{stats['R²']:.4f}")
                                    with col2:
                                        st.metric("RMSE", f"{stats['RMSE']:.4f}")
                                    with col3:
                                        st.metric("Intercept (b₀)", f"{stats['b0']:.4f}")
                                    with col4:
                                        st.metric("Coefficient b₁", f"{stats['b1']:.4f}")
                                    with col5:
                                        st.metric("Coefficient b₂", f"{stats['b2']:.4f}")

                                    # Display equation
                                    st.markdown(f"**Fitted Plane Equation:**")
                                    # Format equation for LaTeX display (fix escape sequences)
                                    latex_eq = equation.replace('Z', z_axis.replace('_', r'\_')).replace('X', x_axis.replace('_', r'\_')).replace('Y', y_axis.replace('_', r'\_')).replace('·', r'\cdot ')
                                    st.latex(latex_eq)
                                    st.code(equation, language=None)

                                    # Explanation
                                    st.info(f"ℹ️ The R² value of {stats['R²']:.4f} indicates that {stats['R²']*100:.2f}% of the variance in {z_axis} is explained by the linear relationship with {x_axis} and {y_axis}. Points are colored by their absolute residual (distance from the fitted plane).")

                                    # Skip the normal layout update section for 3D plots
                                    layout_height = None  # Flag to skip standard layout updates

                                except Exception as e:
                                    st.error(f"Error performing 3D OLS regression: {str(e)}")
                                    st.info("Make sure you have at least 3 valid data points and that X and Y are not perfectly collinear.")
                                    layout_height = None

                            elif selected_chart == 'Heatmap':
                                # Create correlation heatmap
                                corr_matrix = plot_data[numeric_cols].corr()
                                fig = px.imshow(
                                    corr_matrix,
                                    title="Properties Correlation Heatmap",
                                    color_continuous_scale='RdBu_r',
                                    aspect='auto',
                                    text_auto=True
                                )
                                fig.update_traces(texttemplate='%{z:.2f}', textfont_size=10)

                            # Enhance layout (skip for 3D OLS Regression)
                            if selected_chart != '3D OLS Regression':
                                layout_height = 600

                                fig.update_layout(
                                    height=layout_height,
                                    plot_bgcolor='rgba(0,0,0,0)',
                                    paper_bgcolor='rgba(0,0,0,0)',
                                    showlegend=True if color_col else False
                                )

                                if selected_chart != 'Heatmap':
                                    if x_axis:
                                        fig.update_xaxes(title=x_axis.replace('_', ' '))
                                    if y_axis:
                                        fig.update_yaxes(title=y_axis.replace('_', ' '))

                                st.plotly_chart(fig, use_container_width=True)

                            # Show additional statistics if relevant
                            if selected_chart in ['Histogram', 'Box Plot', 'Violin Plot'] and y_axis:
                                col_stats1, col_stats2, col_stats3, col_stats4 = st.columns(4)
                                with col_stats1:
                                    st.metric("Mean", f"{plot_data[y_axis].mean():.3f}")
                                with col_stats2:
                                    st.metric("Median", f"{plot_data[y_axis].median():.3f}")
                                with col_stats3:
                                    st.metric("Std Dev", f"{plot_data[y_axis].std():.3f}")
                                with col_stats4:
                                    st.metric("Range", f"{plot_data[y_axis].max() - plot_data[y_axis].min():.3f}")

                        except Exception as e:
                            st.error(f"Error generating visualization: {str(e)}")
                            st.info("Please check your data and selected parameters.")

            # Download results
            csv = st.session_state.batch_final_df.to_csv(index=False)
            # Extract original filename without extension and add suffix
            original_name = uploaded_file.name.rsplit('.', 1)[0]
            download_filename = f"{original_name}_Calculated_Properties.csv"
            st.download_button(
                label="Download Results as CSV",
                data=csv,
                file_name=download_filename,
                mime="text/csv"
            )

elif input_mode == "Data Visualization":
    st.header("Data Visualization")
    st.markdown("Upload any CSV file to create interactive visualizations of your data.")

    viz_uploaded_file = st.file_uploader("Upload CSV file for visualization", type=['csv'], key="viz_upload")

    if viz_uploaded_file is not None:
        try:
            viz_df = pd.read_csv(viz_uploaded_file)
        except Exception as e:
            st.error(f"Error reading file: {str(e)}")
            st.stop()

        # Data preprocessing: Convert string numbers to numeric
        st.subheader("Data Preprocessing")

        # Store original column info before preprocessing
        original_dtypes = viz_df.dtypes.to_dict()
        converted_columns = []

        # Try to convert string columns that might contain numbers
        for col in viz_df.columns:
            if viz_df[col].dtype == 'object':  # String/object columns
                # Try to convert to numeric, keeping original if fails
                numeric_converted = pd.to_numeric(viz_df[col], errors='coerce')

                # If conversion was successful for most values (>50%), use the numeric version
                non_null_original = viz_df[col].notna().sum()
                non_null_converted = numeric_converted.notna().sum()

                if non_null_converted > 0 and (non_null_converted / non_null_original) >= 0.5:
                    viz_df[col] = numeric_converted
                    converted_columns.append(col)

        # Show preprocessing results
        if converted_columns:
            st.success(f"✅ Converted {len(converted_columns)} columns to numeric: {', '.join(converted_columns)}")

            # Show conversion details in an expander
            with st.expander("📊 View Conversion Details"):
                for col in converted_columns:
                    original_type = str(original_dtypes[col])
                    new_type = str(viz_df[col].dtype)
                    non_null_count = viz_df[col].notna().sum()
                    total_count = len(viz_df)
                    st.write(f"• **{col}**: {original_type} → {new_type} ({non_null_count}/{total_count} values converted)")
        else:
            st.info("ℹ️ No additional numeric conversions needed.")

        st.subheader("Data Preview")
        st.info(f"Loaded {len(viz_df)} rows and {len(viz_df.columns)} columns from {viz_uploaded_file.name}")
        st.dataframe(viz_df.head(), use_container_width=True)

        # Show data info
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("**Data Summary:**")
            st.write(f"• Total rows: {len(viz_df)}")
            st.write(f"• Total columns: {len(viz_df.columns)}")
            numeric_cols_viz = viz_df.select_dtypes(include=[np.number]).columns.tolist()
            st.write(f"• Numeric columns: {len(numeric_cols_viz)}")
            st.write(f"• Text columns: {len(viz_df.columns) - len(numeric_cols_viz)}")

        with col2:
            if len(numeric_cols_viz) > 0:
                st.markdown("**Numeric Columns:**")
                for col in numeric_cols_viz[:10]:  # Show first 10
                    st.write(f"• {col}")
                if len(numeric_cols_viz) > 10:
                    st.write(f"• ... and {len(numeric_cols_viz) - 10} more")

        # Enhanced Interactive Visualization System for standalone mode
        all_cols_viz = viz_df.columns.tolist()

        if len(numeric_cols_viz) >= 1:
            st.subheader("🎨 Interactive Data Visualization")

            # Chart type selection
            chart_types_viz = {
                'Scatter Plot': {'requires': ['x', 'y'], 'desc': 'Compare two numeric variables'},
                'Histogram': {'requires': ['x'], 'desc': 'Distribution of a single variable'},
                'Box Plot': {'requires': ['y'], 'desc': 'Statistical summary with quartiles'},
                'Violin Plot': {'requires': ['y'], 'desc': 'Distribution shape with density'},
                'Line Plot': {'requires': ['x', 'y'], 'desc': 'Trends over continuous data'},
                'Heatmap': {'requires': [], 'desc': 'Correlation matrix of all numeric variables'}
            }

            col1_viz, col2_viz = st.columns([1, 2])
            with col1_viz:
                selected_chart_viz = st.selectbox(
                    "📊 Chart Type:",
                    options=list(chart_types_viz.keys()),
                    key="chart_type_select_viz"
                )
            with col2_viz:
                st.info(f"💡 {chart_types_viz[selected_chart_viz]['desc']}")

            # Dynamic input controls based on chart type
            requires_viz = chart_types_viz[selected_chart_viz]['requires']

            # Create input columns
            input_cols_viz = st.columns(4)

            # Initialize variables
            x_axis_viz, y_axis_viz, color_col_viz, size_col_viz = None, None, None, None

            with input_cols_viz[0]:
                if 'x' in requires_viz:
                    x_axis_viz = st.selectbox(
                        "🔢 X-axis:",
                        options=numeric_cols_viz,
                        index=0,
                        key="viz_x_axis_standalone"
                    )

            with input_cols_viz[1]:
                if 'y' in requires_viz:
                    available_y_viz = [col for col in numeric_cols_viz if col != x_axis_viz] if x_axis_viz else numeric_cols_viz
                    y_axis_viz = st.selectbox(
                        "📈 Y-axis:",
                        options=available_y_viz,
                        index=0 if available_y_viz else 0,
                        key="viz_y_axis_standalone"
                    )
                elif selected_chart_viz in ['Histogram', 'Box Plot', 'Violin Plot']:
                    y_axis_viz = st.selectbox(
                        "📊 Variable:",
                        options=numeric_cols_viz,
                        index=0,
                        key="viz_single_var_standalone"
                    )

            with input_cols_viz[2]:
                # Color mapping options
                color_options_viz = ['None'] + all_cols_viz
                color_col_viz = st.selectbox(
                    "🎨 Color by:",
                    options=color_options_viz,
                    key="viz_color_col_standalone"
                )
                if color_col_viz == 'None':
                    color_col_viz = None

            with input_cols_viz[3]:
                # Additional options based on chart type
                if selected_chart_viz == 'Scatter Plot':
                    show_trendline_viz = st.checkbox("📈 Trend Line", value=False, key="viz_trendline_standalone")
                elif selected_chart_viz == 'Histogram':
                    bins_viz = st.slider("📊 Bins:", min_value=5, max_value=50, value=20, key="viz_bins_standalone")
                elif selected_chart_viz in ['Box Plot', 'Violin Plot']:
                    show_points_viz = st.checkbox("🔹 Show Points", value=True, key="viz_points_standalone")

            # Color scale options (if color column is selected)
            color_scale_viz = 'viridis'  # Default
            if color_col_viz and color_col_viz in numeric_cols_viz:
                st.markdown("**🌈 Color Scale Options:**")
                color_scale_options_viz = {
                    'Viridis': 'viridis',
                    'Plasma': 'plasma',
                    'Inferno': 'inferno',
                    'Rainbow': 'rainbow',
                    'RdYlBu': 'RdYlBu',
                    'Blue-Red': 'RdBu',
                    'Spectral': 'Spectral'
                }
                color_scale_col1_viz, color_scale_col2_viz = st.columns(2)
                with color_scale_col1_viz:
                    color_scale_viz = st.selectbox(
                        "Color Scale:",
                        options=list(color_scale_options_viz.values()),
                        format_func=lambda x: [k for k, v in color_scale_options_viz.items() if v == x][0],
                        index=0,
                        key="viz_color_scale_standalone"
                    )
                with color_scale_col2_viz:
                    reverse_scale_viz = st.checkbox("🔄 Reverse Scale", value=False, key="viz_reverse_scale_standalone")
                    if reverse_scale_viz:
                        color_scale_viz = color_scale_viz + '_r'

            # Size parameter for scatter plots
            size_col_viz = None
            if selected_chart_viz == 'Scatter Plot':
                st.markdown("**📏 Size Options:**")
                size_options_viz = ['None'] + numeric_cols_viz
                size_col_viz = st.selectbox(
                    "Size by:",
                    options=size_options_viz,
                    key="viz_size_col_standalone"
                )
                if size_col_viz == 'None':
                    size_col_viz = None

            # Generate visualization
            if True:  # Always generate visualization
                try:
                    # Prepare data
                    plot_data_viz = viz_df.copy()

                    # Handle color and size columns
                    color_param_viz = {}
                    if color_col_viz:
                        if color_col_viz in numeric_cols_viz:
                            # Continuous color scale
                            color_param_viz = {
                                'color': color_col_viz,
                                'color_continuous_scale': color_scale_viz
                            }
                        else:
                            # Discrete color scale
                            color_param_viz = {'color': color_col_viz}

                    # Add size parameter for scatter plots
                    if size_col_viz and selected_chart_viz == 'Scatter Plot':
                        color_param_viz['size'] = size_col_viz

                    # Create the selected chart
                    if selected_chart_viz == 'Scatter Plot' and x_axis_viz and y_axis_viz:
                        fig = px.scatter(
                            plot_data_viz,
                            x=x_axis_viz,
                            y=y_axis_viz,
                            title=f"{x_axis_viz.replace('_', ' ')} vs {y_axis_viz.replace('_', ' ')}",
                            trendline="ols" if show_trendline_viz else None,
                            **color_param_viz
                        )

                        if not color_col_viz:
                            # Add correlation coefficient and regression stats
                            correlation = plot_data_viz[x_axis_viz].corr(plot_data_viz[y_axis_viz])

                            # Calculate regression statistics if trendline is shown
                            if show_trendline_viz:

                                # Prepare data for regression
                                X = plot_data_viz[[x_axis_viz]].dropna()
                                y = plot_data_viz[y_axis_viz].dropna()

                                # Ensure same indices
                                common_idx = X.index.intersection(y.index)
                                X = X.loc[common_idx]
                                y = y.loc[common_idx]

                                if len(X) > 1:
                                    # Fit linear regression
                                    reg = LinearRegression()
                                    reg.fit(X, y)

                                    # Calculate R²
                                    y_pred = reg.predict(X)
                                    r2 = r2_score(y, y_pred)

                                    # Get slope and intercept
                                    slope = reg.coef_[0]
                                    intercept = reg.intercept_

                                    # Format equation with actual variable names
                                    sign = "+" if intercept >= 0 else "-"
                                    y_name = y_axis_viz.replace('_', ' ')
                                    x_name = x_axis_viz.replace('_', ' ')
                                    equation = f"{y_name} = {slope:.3f} × {x_name} {sign} {abs(intercept):.3f}"

                                    # Add regression info
                                    fig.add_annotation(
                                        x=0.02, y=0.98,
                                        xref="paper", yref="paper",
                                        text=f"<b>Regression Statistics:</b><br>" +
                                             f"Equation: {equation}<br>" +
                                             f"R² = {r2:.3f}<br>" +
                                             f"Correlation: {correlation:.3f}",
                                        showarrow=False,
                                        bgcolor="rgba(255,255,255,0.9)",
                                        bordercolor="black",
                                        borderwidth=1,
                                        align="left",
                                        font=dict(size=11)
                                    )
                            else:
                                # Just show correlation if no trendline
                                fig.add_annotation(
                                    x=0.02, y=0.98,
                                    xref="paper", yref="paper",
                                    text=f"Correlation: {correlation:.3f}",
                                    showarrow=False,
                                    bgcolor="rgba(255,255,255,0.8)",
                                    bordercolor="black",
                                    borderwidth=1
                                )

                    elif selected_chart_viz == 'Histogram' and y_axis_viz:
                        fig = px.histogram(
                            plot_data_viz,
                            x=y_axis_viz,
                            title=f"Distribution of {y_axis_viz.replace('_', ' ')}",
                            nbins=bins_viz,
                            **color_param_viz
                        )

                        # Add mean and median lines
                        mean_val = plot_data_viz[y_axis_viz].mean()
                        median_val = plot_data_viz[y_axis_viz].median()
                        fig.add_vline(x=mean_val, line_dash="dash", line_color="red",
                                    annotation_text=f"Mean: {mean_val:.2f}")
                        fig.add_vline(x=median_val, line_dash="dot", line_color="blue",
                                    annotation_text=f"Median: {median_val:.2f}")

                    elif selected_chart_viz == 'Box Plot' and y_axis_viz:
                        fig = px.box(
                            plot_data_viz,
                            y=y_axis_viz,
                            title=f"Box Plot of {y_axis_viz.replace('_', ' ')}",
                            points="all" if show_points_viz else False,
                            **color_param_viz
                        )

                    elif selected_chart_viz == 'Violin Plot' and y_axis_viz:
                        fig = px.violin(
                            plot_data_viz,
                            y=y_axis_viz,
                            title=f"Violin Plot of {y_axis_viz.replace('_', ' ')}",
                            points="all" if show_points_viz else False,
                            **color_param_viz
                        )

                    elif selected_chart_viz == 'Line Plot' and x_axis_viz and y_axis_viz:
                        # Sort data for line plot
                        plot_data_sorted_viz = plot_data_viz.sort_values(x_axis_viz)
                        fig = px.line(
                            plot_data_sorted_viz,
                            x=x_axis_viz,
                            y=y_axis_viz,
                            title=f"{y_axis_viz.replace('_', ' ')} over {x_axis_viz.replace('_', ' ')}",
                            **color_param_viz
                        )

                    elif selected_chart_viz == 'Heatmap':
                        # Create correlation heatmap
                        corr_matrix = plot_data_viz[numeric_cols_viz].corr()
                        fig = px.imshow(
                            corr_matrix,
                            title="Properties Correlation Heatmap",
                            color_continuous_scale='RdBu_r',
                            aspect='auto',
                            text_auto=True
                        )
                        fig.update_traces(texttemplate='%{z:.2f}', textfont_size=10)

                    # Enhance layout
                    layout_height_viz = 600

                    fig.update_layout(
                        height=layout_height_viz,
                        plot_bgcolor='rgba(0,0,0,0)',
                        paper_bgcolor='rgba(0,0,0,0)',
                        showlegend=True if color_col_viz else False
                    )

                    if selected_chart_viz != 'Heatmap':
                        if x_axis_viz:
                            fig.update_xaxes(title=x_axis_viz.replace('_', ' '))
                        if y_axis_viz:
                            fig.update_yaxes(title=y_axis_viz.replace('_', ' '))

                    st.plotly_chart(fig, use_container_width=True)

                    # Show additional statistics if relevant
                    if selected_chart_viz in ['Histogram', 'Box Plot', 'Violin Plot'] and y_axis_viz:
                        col_stats1, col_stats2, col_stats3, col_stats4 = st.columns(4)
                        with col_stats1:
                            st.metric("Mean", f"{plot_data_viz[y_axis_viz].mean():.3f}")
                        with col_stats2:
                            st.metric("Median", f"{plot_data_viz[y_axis_viz].median():.3f}")
                        with col_stats3:
                            st.metric("Std Dev", f"{plot_data_viz[y_axis_viz].std():.3f}")
                        with col_stats4:
                            st.metric("Range", f"{plot_data_viz[y_axis_viz].max() - plot_data_viz[y_axis_viz].min():.3f}")

                except Exception as e:
                    st.error(f"Error generating visualization: {str(e)}")
                    st.info("Please check your data and selected parameters.")
        else:
            st.warning("No numeric columns found in the data. Visualizations require at least one numeric column.")
            st.info("The uploaded file should contain numeric data for creating charts.")

elif input_mode == "3D Regression Analysis":
    st.header("📐 3D OLS Regression Analysis")
    st.markdown("Perform 3D Ordinary Least Squares regression: **Z = b₀ + b₁·X + b₂·Y**")
    st.markdown("Perfect for SAR (Structure-Activity Relationship) analysis and multi-variate modeling")

    from regression_3d import perform_3d_regression, RegressionSummary

    # File upload
    reg_file = st.file_uploader("Upload CSV file with your data", type=['csv'], key="reg_upload")

    if reg_file is not None:
        try:
            reg_df = pd.read_csv(reg_file)
        except Exception as e:
            st.error(f"Error reading file: {str(e)}")
            st.stop()

        # Data preprocessing
        for col in reg_df.columns:
            if reg_df[col].dtype == 'object':
                numeric_converted = pd.to_numeric(reg_df[col], errors='coerce')
                non_null_original = reg_df[col].notna().sum()
                non_null_converted = numeric_converted.notna().sum()
                if non_null_converted > 0 and (non_null_converted / non_null_original) >= 0.5:
                    reg_df[col] = numeric_converted

        st.subheader("Data Preview")
        st.info(f"Loaded {len(reg_df)} rows and {len(reg_df.columns)} columns from {reg_file.name}")
        st.dataframe(reg_df.head(10), use_container_width=True)

        # Get numeric columns
        numeric_cols_reg = reg_df.select_dtypes(include=[np.number]).columns.tolist()

        if len(numeric_cols_reg) >= 3:
            st.markdown("---")
            st.subheader("🔧 Variable Selection")
            st.markdown("Select three numeric variables for regression: **Z = b₀ + b₁·X + b₂·Y**")

            col1, col2, col3 = st.columns(3)

            with col1:
                x_var = st.selectbox(
                    "📊 Independent Variable 1 (X):",
                    options=numeric_cols_reg,
                    index=0,
                    key="reg_x_var",
                    help="First predictor variable"
                )

            with col2:
                available_y = [col for col in numeric_cols_reg if col != x_var]
                y_var = st.selectbox(
                    "📊 Independent Variable 2 (Y):",
                    options=available_y,
                    index=0 if available_y else 0,
                    key="reg_y_var",
                    help="Second predictor variable"
                )

            with col3:
                available_z = [col for col in numeric_cols_reg if col != x_var and col != y_var]
                z_var = st.selectbox(
                    "🎯 Dependent Variable (Z):",
                    options=available_z,
                    index=0 if available_z else 0,
                    key="reg_z_var",
                    help="Response variable to be predicted"
                )

            # Add example suggestion for SAR data
            st.info(f"💡 **Selected model:** {z_var} = b₀ + b₁·{x_var} + b₂·{y_var}")

            # Perform regression button
            if st.button("🚀 Perform 3D OLS Regression", type="primary"):
                try:
                    with st.spinner("Fitting 3D OLS regression model..."):
                        # Perform regression
                        model, summary = perform_3d_regression(reg_df, x_var, y_var, z_var)

                    st.success("✅ Regression analysis complete!")

                    # Display statsmodels-style summary
                    st.markdown("---")
                    st.subheader("📋 OLS Regression Results")

                    # Text summary
                    with st.expander("📊 Full Statistical Summary (statsmodels style)", expanded=True):
                        st.code(summary.get_summary_text(), language=None)

                    # Coefficient table
                    st.markdown("### 📈 Coefficient Table")
                    coef_df = summary.get_summary_dataframe()
                    st.dataframe(coef_df, use_container_width=True)

                    # Key metrics in columns
                    st.markdown("### 🎯 Key Model Statistics")
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.metric("R² (Goodness of Fit)", f"{summary.r_squared:.4f}")
                        quality = "Excellent" if summary.r_squared >= 0.9 else "Good" if summary.r_squared >= 0.7 else "Moderate" if summary.r_squared >= 0.5 else "Poor"
                        st.caption(f"Quality: {quality}")
                    with col2:
                        st.metric("Adjusted R²", f"{summary.adj_r_squared:.4f}")
                        st.caption("Adjusted for # of predictors")
                    with col3:
                        st.metric("RMSE", f"{summary.model.rmse:.4f}")
                        st.caption("Root Mean Squared Error")
                    with col4:
                        st.metric("F-statistic", f"{summary.f_statistic:.2f}")
                        st.caption(f"p-value: {summary.f_pvalue:.2e}")

                    # Equation display
                    st.markdown("### ✏️ Fitted Plane Equation")
                    equation = model.get_equation_string(decimals=3)
                    equation_display = equation.replace('Z', z_var).replace('X', x_var).replace('Y', y_var)
                    st.code(equation_display, language=None)

                    # Interpretation
                    st.markdown("### 📖 Interpretation")
                    st.markdown(f"""
                    **Model Equation:** `{equation_display}`

                    **What this means:**
                    - **Intercept (b₀ = {model.b0:.3f})**: Baseline value of {z_var} when both {x_var} and {y_var} are zero
                    - **{x_var} coefficient (b₁ = {model.b1:.3f})**: For each 1-unit increase in {x_var}, {z_var} changes by {model.b1:.3f} (holding {y_var} constant)
                    - **{y_var} coefficient (b₂ = {model.b2:.3f})**: For each 1-unit increase in {y_var}, {z_var} changes by {model.b2:.3f} (holding {x_var} constant)
                    - **R² = {summary.r_squared:.3f}**: This model explains {summary.r_squared*100:.1f}% of the variance in {z_var}
                    """)

                    # Coefficient significance
                    st.markdown("### 🔍 Statistical Significance of Coefficients")
                    sig_data = []
                    for i, (var_name, coef, se, t, p) in enumerate([
                        ('Intercept', model.b0, summary.se_b0, summary.t_b0, summary.p_b0),
                        (x_var, model.b1, summary.se_b1, summary.t_b1, summary.p_b1),
                        (y_var, model.b2, summary.se_b2, summary.t_b2, summary.p_b2)
                    ]):
                        significance = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                        sig_label = "Highly significant" if p < 0.001 else "Significant" if p < 0.01 else "Moderately significant" if p < 0.05 else "Not significant"
                        sig_data.append({
                            'Variable': var_name,
                            'Coefficient': f"{coef:.4f}",
                            'p-value': f"{p:.4e}",
                            'Significance': significance,
                            'Interpretation': sig_label
                        })

                    sig_df = pd.DataFrame(sig_data)
                    st.dataframe(sig_df, use_container_width=True)
                    st.caption("Significance codes: *** p<0.001, ** p<0.01, * p<0.05, ns: not significant")

                    # 3D Visualization
                    st.markdown("---")
                    st.markdown("### 🎨 3D Visualization")

                    # Get plane mesh
                    X_mesh, Y_mesh, Z_mesh = model.get_plane_mesh(num_points=25)

                    # Create 3D scatter plot
                    fig = go.Figure()

                    # Color by residuals
                    residuals_abs = np.abs(model.residuals)

                    # Add data points
                    fig.add_trace(go.Scatter3d(
                        x=model.x,
                        y=model.y,
                        z=model.z,
                        mode='markers',
                        marker=dict(
                            size=7,
                            color=residuals_abs,
                            colorscale='Reds',
                            showscale=True,
                            colorbar=dict(title='|Residual|', x=1.1)
                        ),
                        name='Data Points',
                        text=[f"{x_var}: {x:.3f}<br>{y_var}: {y:.3f}<br>{z_var}: {z:.3f}<br>Residual: {r:.3f}"
                              for x, y, z, r in zip(model.x, model.y, model.z, model.residuals)],
                        hovertemplate='%{text}<extra></extra>'
                    ))

                    # Add fitted plane
                    fig.add_trace(go.Surface(
                        x=X_mesh,
                        y=Y_mesh,
                        z=Z_mesh,
                        opacity=0.7,
                        colorscale='Blues',
                        showscale=False,
                        name='Fitted Plane',
                        hovertemplate='Fitted Plane<br>Predicted Value<extra></extra>'
                    ))

                    # Update layout
                    fig.update_layout(
                        title=f"3D OLS Regression: {z_var} vs {x_var} and {y_var}",
                        scene=dict(
                            xaxis_title=x_var,
                            yaxis_title=y_var,
                            zaxis_title=z_var,
                            camera=dict(
                                eye=dict(x=1.5, y=1.5, z=1.3)
                            )
                        ),
                        height=700,
                        showlegend=True
                    )

                    st.plotly_chart(fig, use_container_width=True)
                    st.caption("💡 Points are colored by their absolute residual (distance from fitted plane). Use mouse to rotate the 3D plot.")

                    # Residual Analysis
                    st.markdown("---")
                    st.markdown("### 📊 Residual Analysis")

                    col1, col2 = st.columns(2)

                    with col1:
                        # Residual distribution
                        fig_res_hist = px.histogram(
                            x=model.residuals,
                            nbins=20,
                            title="Residual Distribution",
                            labels={'x': 'Residual'},
                            color_discrete_sequence=['#4ECDC4']
                        )
                        fig_res_hist.add_vline(x=0, line_dash="dash", line_color="red", annotation_text="Zero")
                        fig_res_hist.update_layout(height=400, showlegend=False)
                        st.plotly_chart(fig_res_hist, use_container_width=True)
                        st.caption("Residuals should be normally distributed around zero")

                    with col2:
                        # Predicted vs Actual
                        predicted = model.predict(model.x, model.y)
                        fig_pred = px.scatter(
                            x=predicted,
                            y=model.z,
                            title="Predicted vs Actual Values",
                            labels={'x': f'Predicted {z_var}', 'y': f'Actual {z_var}'},
                            color_discrete_sequence=['#FF6B6B']
                        )
                        # Add perfect prediction line
                        min_val = min(predicted.min(), model.z.min())
                        max_val = max(predicted.max(), model.z.max())
                        fig_pred.add_trace(go.Scatter(
                            x=[min_val, max_val],
                            y=[min_val, max_val],
                            mode='lines',
                            line=dict(color='gray', dash='dash'),
                            name='Perfect Prediction',
                            showlegend=True
                        ))
                        fig_pred.update_layout(height=400)
                        st.plotly_chart(fig_pred, use_container_width=True)
                        st.caption("Points should lie close to the diagonal line for good fit")

                    # Diagnostic tests summary
                    st.markdown("### 🔬 Diagnostic Tests")
                    diag_col1, diag_col2, diag_col3 = st.columns(3)

                    with diag_col1:
                        st.metric("Durbin-Watson", f"{summary.durbin_watson:.3f}")
                        dw_interpretation = "No autocorrelation" if 1.5 < summary.durbin_watson < 2.5 else "Possible autocorrelation"
                        st.caption(f"{dw_interpretation} (ideal: ~2.0)")

                    with diag_col2:
                        st.metric("Jarque-Bera (JB)", f"{summary.jb_statistic:.3f}")
                        jb_interpretation = "Normal residuals" if summary.jb_pvalue > 0.05 else "Non-normal residuals"
                        st.caption(f"{jb_interpretation} (p={summary.jb_pvalue:.3f})")

                    with diag_col3:
                        st.metric("Condition Number", f"{summary.condition_number:.1f}")
                        cn_interpretation = "Low multicollinearity" if summary.condition_number < 30 else "Moderate" if summary.condition_number < 100 else "High multicollinearity"
                        st.caption(cn_interpretation)

                    # Download results
                    st.markdown("---")
                    st.markdown("### 💾 Export Results")

                    # Create summary report
                    report_lines = [
                        "=" * 80,
                        "3D OLS REGRESSION ANALYSIS REPORT",
                        "=" * 80,
                        f"\nFile: {reg_file.name}",
                        f"Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
                        f"\nModel: {z_var} = b0 + b1·{x_var} + b2·{y_var}",
                        "\n" + summary.get_summary_text(),
                        "\n" + "=" * 80,
                        "END OF REPORT",
                        "=" * 80
                    ]
                    report_text = "\n".join(report_lines)

                    col1, col2 = st.columns(2)
                    with col1:
                        st.download_button(
                            label="📄 Download Full Report (TXT)",
                            data=report_text,
                            file_name=f"3D_OLS_Regression_Report_{z_var}.txt",
                            mime="text/plain"
                        )

                    with col2:
                        # Export coefficients and stats
                        stats_dict = summary.get_statistics_dict()
                        stats_df = pd.DataFrame([stats_dict])
                        st.download_button(
                            label="📊 Download Statistics (CSV)",
                            data=stats_df.to_csv(index=False),
                            file_name=f"3D_OLS_Statistics_{z_var}.csv",
                            mime="text/csv"
                        )

                except Exception as e:
                    st.error(f"Error displaying regression results: {str(e)}")
                    import traceback
                    with st.expander("🐛 Debug Information"):
                        st.code(traceback.format_exc())

        else:
            st.warning(f"⚠️ Need at least 3 numeric columns for 3D regression. Found {len(numeric_cols_reg)} numeric columns.")
            st.info("Please upload a CSV file with at least 3 numeric columns (2 independent variables and 1 dependent variable).")

    else:
        st.info("👆 Upload a CSV file to get started with 3D regression analysis")
        st.markdown("""
        ### 📚 What is 3D OLS Regression?

        3D Ordinary Least Squares regression fits a **plane** to your data in 3-dimensional space:

        **Z = b₀ + b₁·X + b₂·Y**

        Where:
        - **Z** is the dependent variable (what you want to predict)
        - **X** and **Y** are independent variables (predictors)
        - **b₀, b₁, b₂** are coefficients calculated to minimize prediction error

        ### 🎯 Use Cases:

        **SAR Analysis (Structure-Activity Relationship)**:
        - Predict **biological activity (pKi)** from **LogP** and **TPSA**
        - Model **IC50** values from molecular descriptors

        **General Applications**:
        - Multi-variate modeling with two predictors
        - Understanding combined effects of two variables
        - Scientific data analysis requiring 3D relationships

        ### 📊 Output Includes:

        - ✅ Complete statistical summary (like statsmodels)
        - ✅ Coefficient table with p-values and confidence intervals
        - ✅ R², Adjusted R², RMSE, F-statistic
        - ✅ Interactive 3D visualization
        - ✅ Residual analysis plots
        - ✅ Diagnostic tests (Durbin-Watson, Jarque-Bera)
        - ✅ Downloadable report and statistics

        ### 🚀 Example Data Structure:

        Your CSV should have at least 3 numeric columns:

        | LogP | TPSA | pKi  |
        |------|------|------|
        | 2.3  | 45.2 | 7.1  |
        | 3.1  | 62.8 | 6.8  |
        | 1.9  | 38.5 | 7.5  |
        | ...  | ...  | ...  |

        Upload your data above to start!
        """)

# Information section - Context-aware based on mode
if input_mode == "3D Regression Analysis":
    with st.expander("ℹ️ 3D Regression Analysis - Help & Documentation"):
        st.markdown("""
        ### 📐 3D Ordinary Least Squares (OLS) Regression

        **Mathematical Model:**
        ```
        Z = b₀ + b₁·X + b₂·Y + ε
        ```

        Where:
        - **Z** = Dependent variable (response, what you want to predict)
        - **X, Y** = Independent variables (predictors, features)
        - **b₀** = Intercept (baseline value when X=0, Y=0)
        - **b₁** = Coefficient for X (effect of X on Z, holding Y constant)
        - **b₂** = Coefficient for Y (effect of Y on Z, holding X constant)
        - **ε** = Error term (residual)

        ---

        ### 📊 Statistical Output Explained

        #### **Model Quality Metrics:**
        - **R² (R-squared)**: Proportion of variance explained (0-1, higher is better)
          - 0.9+ = Excellent fit
          - 0.7-0.9 = Good fit
          - 0.5-0.7 = Moderate fit
          - <0.5 = Poor fit

        - **Adjusted R²**: R² adjusted for number of predictors (penalizes overfitting)

        - **RMSE (Root Mean Squared Error)**: Average prediction error in Z units

        - **F-statistic**: Tests if model is better than just using mean of Z
          - High F-value + low p-value = model is significant

        #### **Coefficient Statistics:**
        - **Coefficient**: Estimated effect size
        - **Std Error**: Uncertainty in coefficient estimate
        - **t-statistic**: Coefficient / Std Error (measures significance)
        - **P>|t|**: Probability that coefficient is actually zero
          - p < 0.001 (***) = Highly significant
          - p < 0.01 (**) = Significant
          - p < 0.05 (*) = Moderately significant
          - p > 0.05 (ns) = Not significant
        - **[0.025, 0.975]**: 95% confidence interval for coefficient

        #### **Diagnostic Tests:**
        - **Durbin-Watson**: Tests for autocorrelation (ideal value ≈ 2.0)
          - 1.5-2.5 = No autocorrelation
          - <1.5 or >2.5 = Possible autocorrelation

        - **Jarque-Bera (JB)**: Tests if residuals are normally distributed
          - p > 0.05 = Residuals are normal (good)
          - p < 0.05 = Non-normal residuals (potential issue)

        - **Condition Number**: Tests for multicollinearity
          - <30 = Low multicollinearity (good)
          - 30-100 = Moderate multicollinearity
          - >100 = High multicollinearity (X and Y are too correlated)

        - **AIC/BIC**: Model comparison metrics (lower is better)

        ---

        ### 🎯 Use Cases

        #### **SAR Analysis (Structure-Activity Relationship):**
        - **Predict biological activity (pKi, IC50)** from molecular descriptors
        - Example: `pKi = b₀ + b₁·LogP + b₂·TPSA`
        - Understand how lipophilicity (LogP) and polarity (TPSA) affect binding

        #### **QSAR (Quantitative Structure-Activity Relationship):**
        - Model drug efficacy from chemical properties
        - Optimize lead compounds
        - Predict activity for new molecules

        #### **General Scientific Data:**
        - Any relationship with 2 independent variables and 1 dependent
        - Physics, chemistry, biology, economics, engineering

        ---

        ### 🎨 Enhanced Visualization Features

        #### **6 Professional Color Schemes:**
        1. **Default** - Classic red/blue
        2. **Professional** - Viridis/YlGnBu
        3. **Dark Mode** - Plasma on dark background
        4. **High Contrast** - Hot/Ice colors
        5. **Publication** - Print-ready styling
        6. **Colorblind Safe** - Accessible Cividis colors

        #### **8 Camera Presets:**
        - **Default**: Standard 3D perspective
        - **Top View**: Looking down Z-axis
        - **Side Views**: X and Y axis perspectives
        - **Isometric**: 45° engineering view
        - **Bird's Eye**: High-angle overview
        - **Close-up/Far View**: Zoom control

        #### **Interactive Features:**
        - **Residual Lines**: Show vertical distance from points to plane
        - **Color by Residuals**: Highlight poorly-fitted points
        - **Adjustable Opacity**: See through plane to view all data
        - **High-Res Export**: PNG at 3200x2400 for publications
        - **Interactive HTML**: Share rotatable 3D plots

        ---

        ### 📝 OLS Assumptions (Check These!)

        1. **Linearity**: Z has a linear relationship with X and Y
        2. **Independence**: Observations are independent
        3. **Homoscedasticity**: Constant variance of residuals
        4. **Normality**: Residuals are normally distributed
        5. **No Multicollinearity**: X and Y are not too correlated

        **How to Check:**
        - View residual plots (should be random scatter)
        - Check Jarque-Bera test (p > 0.05 for normality)
        - Check Condition Number (<30 for low multicollinearity)
        - Look at residual distribution histogram (bell-shaped)

        ---

        ### 💡 Tips for Best Results

        1. **Data Preparation:**
           - Remove outliers or use robust regression
           - Check for missing values
           - Ensure sufficient sample size (n ≥ 30 recommended)

        2. **Variable Selection:**
           - Choose X and Y that make scientific sense
           - Avoid highly correlated X and Y
           - Center/scale variables if very different ranges

        3. **Interpretation:**
           - Always check R² AND residual plots
           - Look for patterns in residuals (shouldn't be any)
           - Consider coefficient significance (p-values)
           - Report confidence intervals for coefficients

        4. **Visualization:**
           - Try different camera angles to see all data
           - Enable residual lines to spot outliers
           - Use colorblind-safe scheme for presentations
           - Export high-res for publications

        ---

        ### 📚 References

        **Mathematical Foundation:**
        - Based on: https://sapiencespace.com/breaking-down-3d-linear-regression/
        - Classic OLS theory with analytical solution

        **Statistical Tests:**
        - Durbin-Watson: Tests for serial correlation
        - Jarque-Bera: Normality test (Jarque & Bera, 1980)
        - F-test: Overall model significance

        **Software:**
        - Similar output to statsmodels (Python)
        - Compatible with R lm() function output

        ---

        ### 🆘 Troubleshooting

        **"Singular matrix" or "Perfect multicollinearity" error:**
        - X and Y are perfectly correlated (one is a multiple of the other)
        - Solution: Choose different variables

        **Very low R² (<0.3):**
        - Linear model may not fit your data
        - Try transforming variables (log, sqrt)
        - Consider non-linear models

        **High Condition Number (>100):**
        - X and Y are too correlated
        - Solution: Use only one of them, or PCA

        **Non-normal residuals (low JB p-value):**
        - May indicate outliers or wrong model
        - Check for outliers in data
        - Consider robust regression

        ---

        **Developed by:** Yashwanth Reddy for ITR-UIC
        **Part of:** Chemo-Informatics Toolkit
        **Module:** 3D OLS Regression Analysis
        """)
else:
    # Default information for other modes
    with st.expander("ℹ️ Information & Property Explanations"):
        st.markdown(PropertyExplanations.get_explanations())

# Visualization Help Section
with st.expander("📊 Visualization Guide"):
    st.markdown("""
    ### 🎨 Enhanced Visualization System

    **Chart Types Available:**
    - **Scatter Plot**: Compare two numeric variables, with optional trend lines
    - **Histogram**: Show distribution of a single variable with mean/median lines
    - **Box Plot**: Statistical summary with quartiles and outliers
    - **Violin Plot**: Distribution shape with density curves
    - **Line Plot**: Trends over continuous data
    - **Heatmap**: Correlation matrix of all numeric properties

    **Color Mapping Options:**
    - Color by any column (numeric or categorical)
    - Multiple color scales: Viridis, Plasma, Inferno, Rainbow, etc.
    - Reverse color scales available
    - Continuous scales for numeric data, discrete for categorical

    **Advanced Features:**
    - **Size Mapping**: For scatter plots, map point size to any numeric column
    - **Interactive Controls**: Trend lines, bin counts, point visibility
    - **Statistical Overlays**: Mean/median lines, correlation coefficients

    **Use Cases:**
    - **Multi-dimensional Analysis**: Use X, Y, Color, and Size to explore 4 variables simultaneously
    - **Property Relationships**: Use scatter plots with trend lines to find correlations
    - **Distribution Analysis**: Use histograms/box plots to understand property ranges
    """)
