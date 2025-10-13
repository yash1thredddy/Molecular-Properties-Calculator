"""
Ligand Efficiency Indices (LEI) Calculator

This module provides functionality for calculating various ligand efficiency indices
based on the AtlasCBS methodology (Cele Abad-Zapatero, A. Cortes-Cabrera 2013).

Calculates 8 Ligand Efficiency Indices:
- NSEI: Normalized Surface Efficiency Index
- NBEI: Normalized Binding Efficiency Index
- BEI: Binding Efficiency Index
- SEI: Surface Efficiency Index
- nBEI: Alternative binding efficiency index
- mBEI: Molecular binding efficiency index
- LEH: Ligand Efficiency (Hopkins)
- LEP: Ligand Efficiency (Polar atoms)

Developed by: Yashwanth Reddy for ITR-UIC
Part of: Chemo-Informatics Toolkit
Based on: AtlasCBS methodology
"""

import pandas as pd
import numpy as np
import math
from typing import Dict, Any, Optional, List, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors


class DependencyChecker:
    """
    Intelligent dependency checker for LEI calculations.
    Detects existing columns and determines what needs to be calculated.
    """

    # Define column aliases for flexible detection
    COLUMN_ALIASES = {
        'pki': ['pKi', 'pki', 'PKI', 'pIC50', 'pic50', 'p_Ki', 'p_ki'],
        'mw': ['MW', 'mw', 'Molecular_Weight', 'molecular_weight', 'MolWt', 'mol_wt'],
        'tpsa': ['TPSA', 'tpsa', 'PSA', 'psa', 'Polar_Surface_Area', 'polar_surface_area'],
        'heavy_atoms': ['Heavy_Atom_Count', 'heavy_atoms', 'HeavyAtomCount', 'nHeavy', 'n_heavy', 'HeavyAtoms'],
        'smiles': ['SMILES', 'smiles', 'Smiles', 'smi', 'SMI', 'canonical_smiles', 'CANONICAL_SMILES']
    }

    # Define dependencies for each LEI
    LEI_DEPENDENCIES = {
        'NSEI': ['pki', 'polar_atoms'],
        'NBEI': ['pki', 'heavy_atoms'],
        'BEI': ['pki', 'mw'],
        'SEI': ['pki', 'tpsa'],
        'nBEI': ['pki', 'heavy_atoms'],
        'mBEI': ['pki', 'mw'],
        'LEH': ['pki', 'heavy_atoms'],
        'LEP': ['pki', 'polar_atoms']
    }

    @staticmethod
    def detect_column(df: pd.DataFrame, standard_name: str) -> Optional[str]:
        """
        Detect a specific column type in DataFrame (case-insensitive).

        Args:
            df: Pandas DataFrame
            standard_name: Standard name from COLUMN_ALIASES keys

        Returns:
            Actual column name if found, None otherwise
        """
        if standard_name not in DependencyChecker.COLUMN_ALIASES:
            return None

        aliases = DependencyChecker.COLUMN_ALIASES[standard_name]

        for col in df.columns:
            if col in aliases or col.lower() in [alias.lower() for alias in aliases]:
                return col
        return None

    @staticmethod
    def detect_all_columns(df: pd.DataFrame) -> Dict[str, Optional[str]]:
        """
        Detect all relevant columns in DataFrame.

        Args:
            df: Pandas DataFrame

        Returns:
            Dictionary mapping standard names to actual column names
            e.g., {'pki': 'pKi', 'mw': 'Molecular_Weight', 'smiles': 'SMILES'}
        """
        detected = {}
        for standard_name in DependencyChecker.COLUMN_ALIASES.keys():
            detected[standard_name] = DependencyChecker.detect_column(df, standard_name)
        return detected

    @staticmethod
    def check_lei_dependencies(df: pd.DataFrame, selected_leis: List[str], manual_mappings: Dict[str, str] = None) -> Dict[str, Any]:
        """
        Check dependencies for selected LEI calculations.

        Args:
            df: Pandas DataFrame
            selected_leis: List of LEI names to calculate
            manual_mappings: Optional manual column mappings to override auto-detection

        Returns:
            Dictionary with dependency status:
            {
                'detected_columns': {...},
                'available_deps': [...],
                'missing_deps': [...],
                'can_calculate': [...],
                'cannot_calculate': [...],
                'needs_calculation': [...],
                'status_by_lei': {...}
            }
        """
        detected_cols = DependencyChecker.detect_all_columns(df)

        # Override with manual mappings if provided
        if manual_mappings:
            for std_name, col_name in manual_mappings.items():
                detected_cols[std_name] = col_name

        result = {
            'detected_columns': detected_cols,
            'available_deps': [],
            'missing_deps': [],
            'can_calculate': [],
            'cannot_calculate': [],
            'needs_calculation': [],
            'status_by_lei': {}
        }

        # Track what dependencies are available
        available_deps = set()
        for dep_name, col_name in detected_cols.items():
            if col_name is not None:
                available_deps.add(dep_name)

        result['available_deps'] = list(available_deps)

        # Check if we have SMILES (needed to calculate missing properties)
        has_smiles = detected_cols['smiles'] is not None

        # Check each selected LEI
        for lei in selected_leis:
            if lei not in DependencyChecker.LEI_DEPENDENCIES:
                continue

            required_deps = DependencyChecker.LEI_DEPENDENCIES[lei]
            lei_status = {
                'required': required_deps,
                'available': [],
                'missing': [],
                'needs_calc': [],
                'can_proceed': False
            }

            for dep in required_deps:
                if dep in available_deps:
                    lei_status['available'].append(dep)
                elif dep in ['polar_atoms', 'heavy_atoms', 'mw', 'tpsa']:
                    # These can be calculated from SMILES
                    if has_smiles:
                        lei_status['needs_calc'].append(dep)
                        if dep not in result['needs_calculation']:
                            result['needs_calculation'].append(dep)
                    else:
                        lei_status['missing'].append(dep)
                else:
                    lei_status['missing'].append(dep)

            # Can proceed if all deps are available or can be calculated
            lei_status['can_proceed'] = len(lei_status['missing']) == 0

            if lei_status['can_proceed']:
                result['can_calculate'].append(lei)
            else:
                result['cannot_calculate'].append(lei)
                for missing in lei_status['missing']:
                    if missing not in result['missing_deps']:
                        result['missing_deps'].append(missing)

            result['status_by_lei'][lei] = lei_status

        return result

    @staticmethod
    def generate_status_message(check_result: Dict[str, Any]) -> str:
        """
        Generate a human-readable status message from dependency check results.

        Args:
            check_result: Result from check_lei_dependencies()

        Returns:
            Formatted status message
        """
        detected = check_result['detected_columns']
        can_calc = check_result['can_calculate']
        cannot_calc = check_result['cannot_calculate']
        needs_calc = check_result['needs_calculation']

        messages = []

        # Show detected columns
        messages.append("📊 **Detected Columns:**")
        for dep_name, col_name in detected.items():
            if col_name:
                messages.append(f"  ✓ {dep_name.upper()}: '{col_name}'")

        # Show what needs to be calculated
        if needs_calc:
            messages.append(f"\n🔧 **Will Calculate from SMILES:**")
            for dep in needs_calc:
                messages.append(f"  → {dep.replace('_', ' ').title()}")

        # Show what's ready
        if can_calc:
            messages.append(f"\n✅ **Ready to Calculate:** {', '.join(can_calc)}")

        # Show what's blocked
        if cannot_calc:
            messages.append(f"\n❌ **Cannot Calculate (missing dependencies):** {', '.join(cannot_calc)}")
            if check_result['missing_deps']:
                messages.append(f"   Missing: {', '.join(check_result['missing_deps'])}")

        return "\n".join(messages)


class LigandEfficiencyCalculator:
    """
    Calculator for Ligand Efficiency Indices (LEIs).
    Based on AtlasCBS methodology.
    """

    # Physical constants
    R_CONST = 0.00198  # Gas constant in kcal/(mol·K)
    T_KELVIN = 300.0   # Temperature in Kelvin

    @staticmethod
    def count_polar_atoms(smiles: str) -> int:
        """
        Count polar atoms (N, O) in a molecule.

        Args:
            smiles: SMILES string

        Returns:
            Number of polar atoms (N + O)
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return 0

            count = 0
            for atom in mol.GetAtoms():
                if atom.GetSymbol() in ['N', 'O']:
                    count += 1
            return count
        except:
            return 0

    @staticmethod
    def count_heavy_atoms(smiles: str) -> int:
        """
        Count heavy atoms (non-hydrogen) in a molecule.

        Args:
            smiles: SMILES string

        Returns:
            Number of heavy atoms
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return 0
            return Descriptors.HeavyAtomCount(mol)
        except:
            return 0

    @staticmethod
    def calculate_ki_from_pki(pki: float) -> float:
        """
        Convert pKi to Ki (Molar).
        Ki = 10^(-pKi)

        Args:
            pki: pKi value

        Returns:
            Ki in Molar
        """
        return math.pow(10, -pki)

    @staticmethod
    def calculate_delta_g(ki: float, temperature: float = T_KELVIN) -> float:
        """
        Calculate Gibbs free energy (ΔG) from Ki.
        ΔG = RT × ln(Ki)

        Args:
            ki: Ki value in Molar
            temperature: Temperature in Kelvin (default: 300K)

        Returns:
            ΔG in kcal/mol
        """
        if ki <= 0:
            return 0.0
        return LigandEfficiencyCalculator.R_CONST * temperature * math.log(ki)

    @staticmethod
    def calculate_lei_values(
        pki: float,
        mw: Optional[float] = None,
        tpsa: Optional[float] = None,
        heavy_atoms: Optional[int] = None,
        polar_atoms: Optional[int] = None,
        pki_error: float = 0.0
    ) -> Dict[str, float]:
        """
        Calculate all Ligand Efficiency Indices.

        Args:
            pki: pKi value (required)
            mw: Molecular weight in Da (optional)
            tpsa: Topological polar surface area in Ų (optional)
            heavy_atoms: Number of heavy atoms (optional)
            polar_atoms: Number of polar atoms (optional)
            pki_error: Error in pKi measurement (default: 0.0)

        Returns:
            Dictionary with calculated LEI values
        """
        results = {}

        # Calculate Ki and ΔG (needed for several LEIs)
        ki = LigandEfficiencyCalculator.calculate_ki_from_pki(pki)
        delta_g = LigandEfficiencyCalculator.calculate_delta_g(ki)
        ki_error = math.log(10) * ki * pki_error if pki_error > 0 else 0.0

        # NSEI: pKi / Polar Atom Count
        if polar_atoms is not None and polar_atoms > 0:
            results['NSEI'] = round(pki / polar_atoms, 4)
            results['NSEI_error'] = round(pki_error / polar_atoms, 4) if pki_error > 0 else 0.0

        # NBEI: pKi / Heavy Atom Count
        if heavy_atoms is not None and heavy_atoms > 0:
            results['NBEI'] = round(pki / heavy_atoms, 4)
            results['NBEI_error'] = round(pki_error / heavy_atoms, 4) if pki_error > 0 else 0.0

        # BEI: pKi / (MW/1000)
        if mw is not None and mw > 0:
            results['BEI'] = round(pki / (mw / 1000.0), 4)
            results['BEI_error'] = round(pki_error / (mw / 1000.0), 4) if pki_error > 0 else 0.0

        # SEI: pKi / (TPSA/100)
        if tpsa is not None and tpsa > 0:
            results['SEI'] = round(pki / (tpsa / 100.0), 4)
            results['SEI_error'] = round(pki_error / (tpsa / 100.0), 4) if pki_error > 0 else 0.0

        # nBEI: -log10(Ki / Heavy Atom Count)
        if heavy_atoms is not None and heavy_atoms > 0:
            ratio = ki / heavy_atoms
            results['nBEI'] = round(-math.log10(ratio), 4)
            results['nBEI_error'] = round(pki_error, 4) if pki_error > 0 else 0.0

        # mBEI: -log10(Ki / MW)
        if mw is not None and mw > 0:
            ratio = ki / mw
            results['mBEI'] = round(-math.log10(ratio), 4)
            results['mBEI_error'] = round(pki_error, 4) if pki_error > 0 else 0.0

        # LEH (Hopkins): -ΔG / Heavy Atom Count
        if heavy_atoms is not None and heavy_atoms > 0:
            results['LEH'] = round(-delta_g / heavy_atoms, 4)
            delta_g_error = LigandEfficiencyCalculator.R_CONST * LigandEfficiencyCalculator.T_KELVIN * ki_error / ki if ki > 0 else 0
            results['LEH_error'] = round(delta_g_error / heavy_atoms, 4) if pki_error > 0 else 0.0

        # LEP (Polar): -ΔG / Polar Atom Count
        if polar_atoms is not None and polar_atoms > 0:
            results['LEP'] = round(-delta_g / polar_atoms, 4)
            delta_g_error = LigandEfficiencyCalculator.R_CONST * LigandEfficiencyCalculator.T_KELVIN * ki_error / ki if ki > 0 else 0
            results['LEP_error'] = round(delta_g_error / polar_atoms, 4) if pki_error > 0 else 0.0

        return results

    @staticmethod
    def calculate_lei_from_row(
        row: pd.Series,
        detected_columns: Dict[str, Optional[str]],
        selected_leis: List[str],
        calculate_missing: bool = True
    ) -> Dict[str, float]:
        """
        Calculate LEIs for a single row in a DataFrame.

        Args:
            row: Pandas Series (DataFrame row)
            detected_columns: Column mapping from DependencyChecker
            selected_leis: List of LEI names to calculate
            calculate_missing: Whether to calculate missing dependencies from SMILES

        Returns:
            Dictionary with calculated LEI values
        """
        # Extract values from row
        pki = None
        if detected_columns['pki'] and detected_columns['pki'] in row.index:
            pki = row[detected_columns['pki']]
            if pd.isna(pki):
                return {}
        else:
            return {}  # pKi is required

        # Get available values from row
        mw = row[detected_columns['mw']] if detected_columns['mw'] and detected_columns['mw'] in row.index else None
        tpsa = row[detected_columns['tpsa']] if detected_columns['tpsa'] and detected_columns['tpsa'] in row.index else None
        heavy_atoms = row[detected_columns['heavy_atoms']] if detected_columns['heavy_atoms'] and detected_columns['heavy_atoms'] in row.index else None

        # Get SMILES if available (for calculating missing values)
        smiles = None
        if detected_columns['smiles'] and detected_columns['smiles'] in row.index:
            smiles = row[detected_columns['smiles']]

        # Calculate missing dependencies from SMILES if needed
        polar_atoms = None
        if calculate_missing and smiles and pd.notna(smiles):
            if heavy_atoms is None or pd.isna(heavy_atoms):
                heavy_atoms = LigandEfficiencyCalculator.count_heavy_atoms(smiles)
            polar_atoms = LigandEfficiencyCalculator.count_polar_atoms(smiles)

        # Handle NaN values
        if mw is not None and pd.isna(mw):
            mw = None
        if tpsa is not None and pd.isna(tpsa):
            tpsa = None
        if heavy_atoms is not None and pd.isna(heavy_atoms):
            heavy_atoms = None

        # Calculate LEIs
        all_leis = LigandEfficiencyCalculator.calculate_lei_values(
            pki=float(pki),
            mw=float(mw) if mw is not None else None,
            tpsa=float(tpsa) if tpsa is not None else None,
            heavy_atoms=int(heavy_atoms) if heavy_atoms is not None else None,
            polar_atoms=int(polar_atoms) if polar_atoms is not None else None
        )

        # Filter to only selected LEIs
        filtered_leis = {}
        for lei in selected_leis:
            if lei in all_leis:
                filtered_leis[lei] = all_leis[lei]
            # Include error columns if they exist
            error_col = f"{lei}_error"
            if error_col in all_leis:
                filtered_leis[error_col] = all_leis[error_col]

        return filtered_leis

    @staticmethod
    def process_batch(
        df: pd.DataFrame,
        selected_leis: List[str],
        show_errors: bool = False,
        manual_mappings: Dict[str, str] = None
    ) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """
        Process a batch DataFrame and calculate selected LEIs.

        Args:
            df: Input DataFrame
            selected_leis: List of LEI names to calculate
            show_errors: Whether to include error columns
            manual_mappings: Optional manual column mappings

        Returns:
            Tuple of (result DataFrame, status dictionary)
        """
        # Check dependencies with manual mappings
        dependency_check = DependencyChecker.check_lei_dependencies(df, selected_leis, manual_mappings)
        detected_columns = dependency_check['detected_columns']

        # Filter to LEIs that can be calculated
        calculable_leis = dependency_check['can_calculate']

        if not calculable_leis:
            return df, {
                'success': False,
                'message': 'No LEIs can be calculated with available data',
                'dependency_check': dependency_check
            }

        # Calculate LEIs for each row
        lei_results = []
        for idx, row in df.iterrows():
            lei_values = LigandEfficiencyCalculator.calculate_lei_from_row(
                row=row,
                detected_columns=detected_columns,
                selected_leis=calculable_leis,
                calculate_missing=True
            )

            # Remove error columns if not requested
            if not show_errors:
                lei_values = {k: v for k, v in lei_values.items() if not k.endswith('_error')}

            lei_results.append(lei_values)

        # Create results DataFrame
        lei_df = pd.DataFrame(lei_results)

        # Combine with original DataFrame
        result_df = pd.concat([df, lei_df], axis=1)

        status = {
            'success': True,
            'calculated_leis': calculable_leis,
            'skipped_leis': [lei for lei in selected_leis if lei not in calculable_leis],
            'dependency_check': dependency_check,
            'rows_processed': len(df)
        }

        return result_df, status


def get_lei_property_group() -> Dict[str, list]:
    """
    Get LEI property group for UI display.

    Returns:
        Dictionary with LEI property group
    """
    return {
        "Ligand Efficiency Indices": [
            'NSEI', 'NBEI', 'BEI', 'SEI',
            'nBEI', 'mBEI', 'LEH', 'LEP'
        ]
    }


def get_lei_descriptions() -> Dict[str, str]:
    """
    Get descriptions for each LEI.

    Returns:
        Dictionary mapping LEI names to descriptions
    """
    return {
        'NSEI': 'Normalized Surface Efficiency Index (pKi / Polar Atoms)',
        'NBEI': 'Normalized Binding Efficiency Index (pKi / Heavy Atoms)',
        'BEI': 'Binding Efficiency Index (pKi / (MW/1000))',
        'SEI': 'Surface Efficiency Index (pKi / (TPSA/100))',
        'nBEI': 'Alternative Binding Efficiency Index (-log10(Ki / Heavy Atoms))',
        'mBEI': 'Molecular Binding Efficiency Index (-log10(Ki / MW))',
        'LEH': 'Ligand Efficiency Hopkins (-ΔG / Heavy Atoms)',
        'LEP': 'Ligand Efficiency Polar (-ΔG / Polar Atoms)'
    }
