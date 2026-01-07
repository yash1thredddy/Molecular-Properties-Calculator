"""Pytest configuration for Tests/ directory.

This file ensures the project root is in the Python path
so tests can import from molecular_calculator.
"""

import sys
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
