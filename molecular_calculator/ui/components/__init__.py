"""UI components package.

This module exports reusable UI components for the Streamlit application.
"""

from .file_uploader import (
    render_file_uploader,
    validate_uploaded_file,
    read_uploaded_file,
    detect_smiles_column,
    detect_name_column,
    render_column_selector,
    render_file_upload_section,
    create_download_button,
    create_excel_download_button,
)

from .property_selector import (
    get_all_property_groups,
    get_all_properties,
    render_property_selector,
    render_compact_property_selector,
    render_property_group_selector,
    render_property_summary,
    get_property_description,
)

from .molecule_viewer import (
    render_structure_viewer_hint,
    embed_structure_viewer,
    render_smiles_input,
    render_format_selector,
    render_molecule_info,
    prepare_chart_customdata,
    get_structure_viewer_guide,
)

from .results_display import (
    render_property_card,
    render_properties_grid,
    render_properties_by_group,
    render_calculation_result,
    render_dataframe_preview,
    render_batch_results,
    render_batch_statistics,
    render_rule_compliance_summary,
    render_processing_progress,
    render_error_summary,
    render_property_explanations,
)

from .charts import (
    create_scatter_plot,
    create_histogram,
    create_box_plot,
    create_violin_plot,
    create_bar_chart,
    create_3d_scatter,
    create_correlation_heatmap,
    create_pair_plot,
    render_chart_with_viewer,
    get_available_chart_types,
    render_chart_controls,
)

from .visualization import (
    render_distribution_plots,
    render_interactive_visualization,
)

from .structure_viewer import (
    get_structure_viewer_component,
    get_structure_viewer_hint,
)

from .threejs_regression import (
    get_threejs_regression_component,
    get_threejs_colorbar,
)

from .plotly_3d_utils import (
    add_axis_arrows,
    add_origin_crosshair,
    add_residual_vectors,
    add_predicted_markers,
    create_fitted_plane_mesh,
)

from .interference_display import (
    render_interference_header,
    render_interference_metrics,
    render_batch_interference_metrics,
    render_flag_summary_table,
    render_flagged_compounds_table,
    render_interference_section,
    calculate_batch_interference_flags,
)

__all__ = [
    # File uploader
    "render_file_uploader",
    "validate_uploaded_file",
    "read_uploaded_file",
    "detect_smiles_column",
    "detect_name_column",
    "render_column_selector",
    "render_file_upload_section",
    "create_download_button",
    "create_excel_download_button",
    # Property selector
    "get_all_property_groups",
    "get_all_properties",
    "render_property_selector",
    "render_compact_property_selector",
    "render_property_group_selector",
    "render_property_summary",
    "get_property_description",
    # Molecule viewer
    "render_structure_viewer_hint",
    "embed_structure_viewer",
    "render_smiles_input",
    "render_format_selector",
    "render_molecule_info",
    "prepare_chart_customdata",
    "get_structure_viewer_guide",
    # Results display
    "render_property_card",
    "render_properties_grid",
    "render_properties_by_group",
    "render_calculation_result",
    "render_dataframe_preview",
    "render_batch_results",
    "render_batch_statistics",
    "render_rule_compliance_summary",
    "render_processing_progress",
    "render_error_summary",
    "render_property_explanations",
    # Charts
    "create_scatter_plot",
    "create_histogram",
    "create_box_plot",
    "create_violin_plot",
    "create_bar_chart",
    "create_3d_scatter",
    "create_correlation_heatmap",
    "create_pair_plot",
    "render_chart_with_viewer",
    "get_available_chart_types",
    "render_chart_controls",
    # Visualization
    "render_distribution_plots",
    "render_interactive_visualization",
    # Structure viewer
    "get_structure_viewer_component",
    "get_structure_viewer_hint",
    # ThreeJS regression
    "get_threejs_regression_component",
    "get_threejs_colorbar",
    # Plotly 3D utilities
    "add_axis_arrows",
    "add_origin_crosshair",
    "add_residual_vectors",
    "add_predicted_markers",
    "create_fitted_plane_mesh",
    # Interference display
    "render_interference_header",
    "render_interference_metrics",
    "render_batch_interference_metrics",
    "render_flag_summary_table",
    "render_flagged_compounds_table",
    "render_interference_section",
    "calculate_batch_interference_flags",
]
