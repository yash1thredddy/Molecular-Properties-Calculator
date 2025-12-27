# ==============================================================================
# Chart Export Utilities
# ==============================================================================
"""
Utilities for exporting Plotly charts to various formats.
Requires kaleido package for static image export.
"""

import io
import logging
from typing import Optional, Literal, Dict, Any

import plotly.graph_objects as go

logger = logging.getLogger(__name__)

# Export format configurations
EXPORT_FORMATS: Dict[str, Dict[str, Any]] = {
    'png': {
        'extension': '.png',
        'mime_type': 'image/png',
        'description': 'PNG Image',
    },
    'svg': {
        'extension': '.svg',
        'mime_type': 'image/svg+xml',
        'description': 'SVG Vector Image',
    },
    'pdf': {
        'extension': '.pdf',
        'mime_type': 'application/pdf',
        'description': 'PDF Document',
    },
    'html': {
        'extension': '.html',
        'mime_type': 'text/html',
        'description': 'Interactive HTML',
    },
}

ExportFormat = Literal['png', 'svg', 'pdf', 'html']


def export_chart_to_bytes(
    fig: go.Figure,
    format: ExportFormat = 'png',
    width: int = 1200,
    height: int = 800,
    scale: float = 2.0,
) -> bytes:
    """
    Export a Plotly figure to bytes in the specified format.

    Args:
        fig: Plotly figure to export
        format: Export format ('png', 'svg', 'pdf', 'html')
        width: Image width in pixels
        height: Image height in pixels
        scale: Scale factor for image resolution

    Returns:
        Bytes of the exported image/document

    Raises:
        ValueError: If format is not supported
        RuntimeError: If export fails
    """
    if format not in EXPORT_FORMATS:
        raise ValueError(f"Unsupported format: {format}. Supported: {list(EXPORT_FORMATS.keys())}")

    try:
        if format == 'html':
            # HTML export doesn't use kaleido
            html_str = fig.to_html(include_plotlyjs='cdn', full_html=True)
            return html_str.encode('utf-8')
        else:
            # Use kaleido for static image export
            img_bytes = fig.to_image(
                format=format,
                width=width,
                height=height,
                scale=scale,
            )
            return img_bytes

    except Exception as e:
        logger.error(f"Failed to export chart to {format}: {e}")
        raise RuntimeError(f"Chart export failed: {e}")


def export_chart(
    fig: go.Figure,
    filename: str,
    format: ExportFormat = 'png',
    width: int = 1200,
    height: int = 800,
    scale: float = 2.0,
) -> str:
    """
    Export a Plotly figure to a file.

    Args:
        fig: Plotly figure to export
        filename: Output filename (extension will be added if missing)
        format: Export format ('png', 'svg', 'pdf', 'html')
        width: Image width in pixels
        height: Image height in pixels
        scale: Scale factor for image resolution

    Returns:
        Path to the exported file

    Raises:
        ValueError: If format is not supported
        RuntimeError: If export fails
    """
    format_info = EXPORT_FORMATS.get(format)
    if not format_info:
        raise ValueError(f"Unsupported format: {format}")

    # Ensure proper extension
    if not filename.endswith(format_info['extension']):
        filename = filename + format_info['extension']

    try:
        img_bytes = export_chart_to_bytes(fig, format, width, height, scale)

        mode = 'w' if format == 'html' else 'wb'
        with open(filename, mode) as f:
            if format == 'html':
                f.write(img_bytes.decode('utf-8'))
            else:
                f.write(img_bytes)

        logger.info(f"Chart exported to {filename}")
        return filename

    except Exception as e:
        logger.error(f"Failed to save chart to {filename}: {e}")
        raise


def get_download_data(
    fig: go.Figure,
    format: ExportFormat = 'png',
    width: int = 1200,
    height: int = 800,
    scale: float = 2.0,
) -> tuple:
    """
    Get chart data ready for Streamlit download button.

    Args:
        fig: Plotly figure to export
        format: Export format
        width: Image width in pixels
        height: Image height in pixels
        scale: Scale factor for image resolution

    Returns:
        Tuple of (bytes_data, mime_type, file_extension)

    Example:
        data, mime, ext = get_download_data(fig, 'png')
        st.download_button(
            label="Download Chart",
            data=data,
            file_name=f"chart{ext}",
            mime=mime
        )
    """
    format_info = EXPORT_FORMATS[format]
    img_bytes = export_chart_to_bytes(fig, format, width, height, scale)

    return (
        img_bytes,
        format_info['mime_type'],
        format_info['extension'],
    )


def create_export_buttons(
    fig: go.Figure,
    st_container,
    base_filename: str = "chart",
    formats: Optional[list] = None,
    width: int = 1200,
    height: int = 800,
) -> None:
    """
    Create Streamlit download buttons for multiple export formats.

    Args:
        fig: Plotly figure to export
        st_container: Streamlit container to add buttons to
        base_filename: Base filename without extension
        formats: List of formats to offer (default: all)
        width: Image width
        height: Image height

    Example:
        create_export_buttons(fig, st.sidebar, "regression_plot")
    """
    if formats is None:
        formats = ['png', 'svg', 'pdf', 'html']

    cols = st_container.columns(len(formats))

    for i, fmt in enumerate(formats):
        format_info = EXPORT_FORMATS.get(fmt)
        if format_info:
            try:
                data, mime, ext = get_download_data(fig, fmt, width, height)
                cols[i].download_button(
                    label=f"{fmt.upper()}",
                    data=data,
                    file_name=f"{base_filename}{ext}",
                    mime=mime,
                    key=f"export_{base_filename}_{fmt}",
                )
            except Exception as e:
                cols[i].error(f"{fmt.upper()} export failed")
                logger.error(f"Export button creation failed for {fmt}: {e}")
