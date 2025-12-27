"""UI package for Streamlit application.

This module exports UI pages and components for the Molecular Properties Calculator.
"""

from .pages import (
    render_single_molecule_page,
    render_batch_processing_page,
)

from .components import (
    # File uploader
    render_file_uploader,
    render_file_upload_section,
    create_download_button,
    create_excel_download_button,
    # Property selector
    render_property_selector,
    get_all_property_groups,
    # Results display
    render_properties_by_group,
    render_batch_results,
    render_property_explanations,
    # Charts
    create_scatter_plot,
    create_3d_scatter,
    render_chart_with_viewer,
    get_available_chart_types,
)

from .theme import (
    ChartTheme,
    ColorPalette,
    DEFAULT_THEME,
    COLOR_SCALES,
    CHART_SIZES,
    apply_theme,
    apply_3d_theme,
    get_color_scale,
    get_categorical_colors,
    get_chart_size,
    inject_custom_css,
)

__all__ = [
    # Pages
    "render_single_molecule_page",
    "render_batch_processing_page",
    # Components
    "render_file_uploader",
    "render_file_upload_section",
    "create_download_button",
    "create_excel_download_button",
    "render_property_selector",
    "get_all_property_groups",
    "render_properties_by_group",
    "render_batch_results",
    "render_property_explanations",
    "create_scatter_plot",
    "create_3d_scatter",
    "render_chart_with_viewer",
    "get_available_chart_types",
    # Theme
    "ChartTheme",
    "ColorPalette",
    "DEFAULT_THEME",
    "COLOR_SCALES",
    "CHART_SIZES",
    "apply_theme",
    "apply_3d_theme",
    "get_color_scale",
    "get_categorical_colors",
    "get_chart_size",
    "inject_custom_css",
]
