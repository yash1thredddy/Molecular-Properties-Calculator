"""Custom Monaco-based formula editor Streamlit component.

v1 (spike): renders a Monaco editor with column + function autocomplete and
returns the formula string. The frontend lives in ``frontend/`` and is built
to ``frontend/build`` (committed for Streamlit Community Cloud, which does not
run ``npm build``).

Dev mode: set ``FORMULA_EDITOR_DEV=1`` and run the Vite dev server
(``npm run dev`` in ``frontend/``, default http://localhost:5173``) for hot reload.
Release mode (default): serves the committed ``frontend/build`` directory.
"""
from __future__ import annotations

import os
from typing import List, Optional

import streamlit.components.v1 as components

_DEV = os.environ.get("FORMULA_EDITOR_DEV") == "1"
_DIR = os.path.dirname(os.path.abspath(__file__))
_BUILD_DIR = os.path.join(_DIR, "frontend", "build")

if _DEV:
    _component = components.declare_component(
        "formula_editor", url="http://localhost:5173"
    )
else:
    _component = components.declare_component("formula_editor", path=_BUILD_DIR)


def formula_editor(
    value: str = "",
    columns: Optional[List[str]] = None,
    functions: Optional[List[str]] = None,
    height: int = 160,
    key: Optional[str] = None,
) -> str:
    """Render the Monaco formula editor and return the current formula string.

    Args:
        value: Initial formula text.
        columns: Column names offered as ``[Column]`` autocomplete entries.
        functions: Function names offered as ``func()`` autocomplete entries.
        height: Editor height in pixels.
        key: Streamlit widget key.

    Returns:
        The current formula text (empty string until the user types).
    """
    return _component(
        value=value,
        columns=columns or [],
        functions=functions or [],
        height=height,
        key=key,
        default=value,
    )
