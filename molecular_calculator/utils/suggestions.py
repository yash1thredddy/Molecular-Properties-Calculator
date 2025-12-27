# ==============================================================================
# Smart Variable Suggestions
# ==============================================================================
"""
Utilities for intelligently suggesting variables based on data characteristics.
Helps users identify relevant columns and correlations.
"""

import logging
from typing import List, Tuple, Optional, Dict, Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ==============================================================================
# Column Type Detection
# ==============================================================================

# Common SMILES column names (case-insensitive matching)
SMILES_COLUMN_PATTERNS = [
    'smiles', 'smile', 'canonical_smiles', 'canonicalsmiles',
    'structure', 'mol', 'molecule', 'compound', 'ligand',
    'smiles_string', 'smi', 'input_smiles', 'input',
]

# Common identifier column names
ID_COLUMN_PATTERNS = [
    'id', 'name', 'compound_id', 'mol_id', 'molecule_id', 'molecule_name',
    'compound_name', 'ligand_name', 'cas', 'pubchem_cid', 'chembl_id',
    'index', 'idx', 'row', 'sample', 'entry',
]

# Common activity/property column names
ACTIVITY_COLUMN_PATTERNS = [
    'ic50', 'ec50', 'ki', 'kd', 'activity', 'potency',
    'pic50', 'pec50', 'pki', 'pkd', 'log',
    'binding', 'affinity', 'inhibition',
]

# Common molecular property patterns
PROPERTY_COLUMN_PATTERNS = [
    'mw', 'molecular_weight', 'molwt', 'mass',
    'logp', 'alogp', 'clogp', 'logd',
    'tpsa', 'psa', 'polar_surface_area',
    'hbd', 'hba', 'h_bond', 'hydrogen',
    'rotatable', 'rot_bonds', 'num_rot',
    'rings', 'aromatic', 'hetero',
    'heavy_atoms', 'num_atoms',
]


def detect_smiles_column(df: pd.DataFrame) -> Optional[str]:
    """
    Detect the most likely SMILES column in a DataFrame.

    Args:
        df: DataFrame to analyze

    Returns:
        Column name if found, None otherwise
    """
    columns_lower = {col.lower(): col for col in df.columns}

    # Check for exact pattern matches
    for pattern in SMILES_COLUMN_PATTERNS:
        if pattern in columns_lower:
            return columns_lower[pattern]

    # Check for partial matches
    for col_lower, col_orig in columns_lower.items():
        for pattern in SMILES_COLUMN_PATTERNS:
            if pattern in col_lower:
                return col_orig

    # Heuristic: look for string columns with typical SMILES characteristics
    for col in df.columns:
        if df[col].dtype == 'object':
            sample = df[col].dropna().head(10)
            if len(sample) > 0:
                # SMILES typically contain C, c, (, ), =, #, etc.
                smiles_chars = set('CcNnOoSsPpFfClBrI()[]=#@+-0123456789')
                is_smiles_like = all(
                    set(str(val)).issubset(smiles_chars) and len(str(val)) > 5
                    for val in sample
                )
                if is_smiles_like:
                    return col

    return None


def detect_id_column(df: pd.DataFrame) -> Optional[str]:
    """
    Detect the most likely identifier column in a DataFrame.

    Args:
        df: DataFrame to analyze

    Returns:
        Column name if found, None otherwise
    """
    columns_lower = {col.lower(): col for col in df.columns}

    for pattern in ID_COLUMN_PATTERNS:
        if pattern in columns_lower:
            return columns_lower[pattern]

    # Partial matches
    for col_lower, col_orig in columns_lower.items():
        for pattern in ID_COLUMN_PATTERNS:
            if pattern in col_lower:
                return col_orig

    return None


def detect_column_type(col_name: str) -> str:
    """
    Categorize a column based on its name.

    Args:
        col_name: Column name to categorize

    Returns:
        Category string: 'smiles', 'id', 'activity', 'property', or 'unknown'
    """
    col_lower = col_name.lower()

    for pattern in SMILES_COLUMN_PATTERNS:
        if pattern in col_lower:
            return 'smiles'

    for pattern in ID_COLUMN_PATTERNS:
        if pattern in col_lower:
            return 'id'

    for pattern in ACTIVITY_COLUMN_PATTERNS:
        if pattern in col_lower:
            return 'activity'

    for pattern in PROPERTY_COLUMN_PATTERNS:
        if pattern in col_lower:
            return 'property'

    return 'unknown'


# ==============================================================================
# Correlation Analysis
# ==============================================================================

def suggest_correlated_pairs(
    df: pd.DataFrame,
    min_correlation: float = 0.5,
    max_pairs: int = 10,
) -> List[Tuple[str, str, float]]:
    """
    Find highly correlated numeric column pairs.

    Args:
        df: DataFrame to analyze
        min_correlation: Minimum absolute correlation threshold
        max_pairs: Maximum number of pairs to return

    Returns:
        List of (col1, col2, correlation) tuples, sorted by absolute correlation
    """
    # Get numeric columns only
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

    if len(numeric_cols) < 2:
        return []

    try:
        corr_matrix = df[numeric_cols].corr()
    except Exception as e:
        logger.warning(f"Failed to compute correlation matrix: {e}")
        return []

    pairs = []
    seen = set()

    for i, col1 in enumerate(numeric_cols):
        for col2 in numeric_cols[i + 1:]:
            if (col1, col2) in seen or (col2, col1) in seen:
                continue

            corr = corr_matrix.loc[col1, col2]
            if pd.notna(corr) and abs(corr) >= min_correlation:
                pairs.append((col1, col2, corr))
                seen.add((col1, col2))

    # Sort by absolute correlation
    pairs.sort(key=lambda x: abs(x[2]), reverse=True)

    return pairs[:max_pairs]


def suggest_regression_variables(
    df: pd.DataFrame,
    target_col: Optional[str] = None,
) -> Dict[str, List[str]]:
    """
    Suggest variables for regression analysis.

    Args:
        df: DataFrame to analyze
        target_col: Optional target column (y variable)

    Returns:
        Dictionary with 'x_candidates', 'y_candidates', 'z_candidates'
    """
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()

    # Categorize columns
    activity_cols = []
    property_cols = []
    other_cols = []

    for col in numeric_cols:
        col_type = detect_column_type(col)
        if col_type == 'activity':
            activity_cols.append(col)
        elif col_type == 'property':
            property_cols.append(col)
        else:
            other_cols.append(col)

    # Y candidates: activity columns are usually targets
    y_candidates = activity_cols + other_cols
    if target_col and target_col in y_candidates:
        y_candidates.remove(target_col)
        y_candidates.insert(0, target_col)

    # X candidates: property columns are usually predictors
    x_candidates = property_cols + other_cols

    # Z candidates: for 3D plots, similar to X
    z_candidates = property_cols + other_cols

    return {
        'x_candidates': x_candidates,
        'y_candidates': y_candidates,
        'z_candidates': z_candidates,
    }


# ==============================================================================
# Variable Recommendations
# ==============================================================================

def recommend_visualization_variables(
    df: pd.DataFrame,
    chart_type: str = 'scatter',
) -> Dict[str, Optional[str]]:
    """
    Recommend variables for a specific chart type.

    Args:
        df: DataFrame to analyze
        chart_type: Type of chart ('scatter', 'histogram', 'box', 'bar', '3d')

    Returns:
        Dictionary with recommended column names
    """
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    categorical_cols = df.select_dtypes(include=['object', 'category']).columns.tolist()

    # Filter out likely SMILES columns from categorical
    smiles_col = detect_smiles_column(df)
    if smiles_col and smiles_col in categorical_cols:
        categorical_cols.remove(smiles_col)

    recommendations: Dict[str, Optional[str]] = {}

    if chart_type == 'scatter':
        # Find best correlated pair
        pairs = suggest_correlated_pairs(df, min_correlation=0.3, max_pairs=1)
        if pairs:
            recommendations['x'] = pairs[0][0]
            recommendations['y'] = pairs[0][1]
        else:
            recommendations['x'] = numeric_cols[0] if numeric_cols else None
            recommendations['y'] = numeric_cols[1] if len(numeric_cols) > 1 else None
        recommendations['color'] = categorical_cols[0] if categorical_cols else None

    elif chart_type == 'histogram':
        # Prefer activity columns for histograms
        suggestions = suggest_regression_variables(df)
        recommendations['x'] = suggestions['y_candidates'][0] if suggestions['y_candidates'] else None

    elif chart_type == 'box':
        recommendations['x'] = categorical_cols[0] if categorical_cols else None
        suggestions = suggest_regression_variables(df)
        recommendations['y'] = suggestions['y_candidates'][0] if suggestions['y_candidates'] else None

    elif chart_type == 'bar':
        recommendations['x'] = categorical_cols[0] if categorical_cols else None
        recommendations['y'] = numeric_cols[0] if numeric_cols else None

    elif chart_type == '3d':
        if len(numeric_cols) >= 3:
            # Find three related columns
            pairs = suggest_correlated_pairs(df, min_correlation=0.2, max_pairs=5)
            used_cols = set()
            for col1, col2, _ in pairs:
                used_cols.add(col1)
                used_cols.add(col2)
                if len(used_cols) >= 3:
                    break

            cols_list = list(used_cols)[:3]
            if len(cols_list) < 3:
                cols_list = numeric_cols[:3]

            recommendations['x'] = cols_list[0]
            recommendations['y'] = cols_list[1]
            recommendations['z'] = cols_list[2] if len(cols_list) > 2 else None

    return recommendations


def get_column_stats(df: pd.DataFrame, col: str) -> Dict[str, Any]:
    """
    Get summary statistics for a column.

    Args:
        df: DataFrame
        col: Column name

    Returns:
        Dictionary with column statistics
    """
    if col not in df.columns:
        return {}

    series = df[col]
    stats: Dict[str, Any] = {
        'name': col,
        'dtype': str(series.dtype),
        'non_null': int(series.count()),
        'null_count': int(series.isna().sum()),
        'null_pct': round(series.isna().mean() * 100, 1),
    }

    if pd.api.types.is_numeric_dtype(series):
        stats.update({
            'min': float(series.min()) if pd.notna(series.min()) else None,
            'max': float(series.max()) if pd.notna(series.max()) else None,
            'mean': float(series.mean()) if pd.notna(series.mean()) else None,
            'median': float(series.median()) if pd.notna(series.median()) else None,
            'std': float(series.std()) if pd.notna(series.std()) else None,
        })
    else:
        stats.update({
            'unique': int(series.nunique()),
            'top': series.mode().iloc[0] if len(series.mode()) > 0 else None,
        })

    return stats
