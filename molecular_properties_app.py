import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from io import StringIO
from molecular_calculator import MolecularCalculator, PropertyExplanations
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


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
input_mode = st.sidebar.radio("Select input mode:", ["Single Molecule", "Batch Processing", "Data Visualization"])

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

        st.subheader("Data Preview")
        st.info(f"Loaded {len(df)} rows from {uploaded_file.name}")
        st.dataframe(df.head(), use_container_width=True)

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

        # Initialize batch results session state
        if 'batch_results_df' not in st.session_state:
            st.session_state.batch_results_df = None
        if 'batch_final_df' not in st.session_state:
            st.session_state.batch_final_df = None

        # Clear previous results when new file is uploaded
        if 'current_batch_file' not in st.session_state:
            st.session_state.current_batch_file = None

        if st.session_state.current_batch_file != uploaded_file.name:
            st.session_state.current_batch_file = uploaded_file.name
            st.session_state.batch_results_df = None
            st.session_state.batch_final_df = None
            st.session_state.batch_selected_properties = set()  # Clear property selections too

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

                # Create results DataFrame and store in session state
                results_df = pd.DataFrame(results)
                final_df = pd.concat([df, results_df], axis=1)

                # Store results in session state
                st.session_state.batch_results_df = results_df
                st.session_state.batch_final_df = final_df

                status_text.text("Processing complete!")
            else:
                st.error("Please select a SMILES column and at least one property.")

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
                numeric_cols = st.session_state.batch_results_df.select_dtypes(include=[np.number]).columns.tolist()
                all_cols = st.session_state.batch_final_df.columns.tolist()

                if len(numeric_cols) >= 1:
                    st.subheader("🎨 Interactive Data Visualization")

                    # Chart type selection
                    chart_types = {
                        'Scatter Plot': {'requires': ['x', 'y'], 'desc': 'Compare two numeric variables'},
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
                    x_axis, y_axis, color_col, size_col = None, None, None, None

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

                            # Enhance layout
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

# Information section
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