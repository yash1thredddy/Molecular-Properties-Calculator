"""
Molecular Properties Calculator Backend

This module provides core functionality for calculating molecular properties
from chemical structures (SMILES, InChI, InChI Key).

Developed by: Yashwanth Reddy for ITR-UIC
Part of: Chemo-Informatics Toolkit
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit import RDLogger
import re
import requests
import time
from typing import Dict, Any, Optional, Union

# Import Ligand Efficiency functionality
from ligand_efficiency import (
    DependencyChecker,
    LigandEfficiencyCalculator,
    get_lei_property_group,
    get_lei_descriptions
)


class MolecularCalculator:
    """Main class for molecular property calculations"""

    @staticmethod
    def suppress_rdkit_warnings(suppress: bool = True):
        """
        Suppress or enable RDKit warning messages

        Args:
            suppress: If True, suppress warnings; if False, enable warnings
        """
        if suppress:
            RDLogger.DisableLog('rdApp.*')
        else:
            RDLogger.EnableLog('rdApp.*')

    @staticmethod
    def convert_inchi_key_to_smiles(inchi_key: str, timeout: int = 10) -> Optional[str]:
        """
        Convert InChI Key to SMILES using Chemical Identifier Resolver (CIR)

        Args:
            inchi_key: InChI Key string
            timeout: Request timeout in seconds

        Returns:
            SMILES string or None if conversion fails
        """
        try:
            # Use NIH Chemical Identifier Resolver
            cir_url = f"https://cactus.nci.nih.gov/chemical/structure/{inchi_key}/smiles"

            response = requests.get(cir_url, timeout=timeout)

            if response.status_code == 200:
                smiles = response.text.strip()
                # Validate the SMILES
                if smiles and not smiles.startswith("Error") and not smiles.startswith("<"):
                    # Test if SMILES is valid
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is not None:
                        return smiles

            # Fallback: Try PubChem API
            pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchi_key}/property/IsomericSMILES/JSON"

            response = requests.get(pubchem_url, timeout=timeout)

            if response.status_code == 200:
                data = response.json()
                if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                    properties = data['PropertyTable']['Properties']
                    if properties and 'IsomericSMILES' in properties[0]:
                        smiles = properties[0]['IsomericSMILES']
                        # Validate the SMILES
                        mol = Chem.MolFromSmiles(smiles)
                        if mol is not None:
                            return smiles

        except Exception:
            pass

        return None

    @staticmethod
    def convert_to_smiles(input_text: str, input_type: str, enable_online_lookup: bool = True) -> Optional[str]:
        """
        Convert InChI or InChI Key to SMILES

        Args:
            input_text: Input molecular structure string
            input_type: Type of input ('smiles', 'inchi', 'inchi_key')
            enable_online_lookup: Allow online database lookup for InChI Keys

        Returns:
            SMILES string or None if conversion fails
        """
        try:
            if input_type.lower() == 'smiles':
                return input_text
            elif input_type.lower() == 'inchi':
                mol = Chem.MolFromInchi(input_text)
                if mol is None:
                    return None

                # Sanitize molecule to resolve stereochemical conflicts
                try:
                    Chem.SanitizeMol(mol)
                except:
                    # If sanitization fails, try without stereo
                    Chem.RemoveStereochemistry(mol)
                    Chem.SanitizeMol(mol)

                return Chem.MolToSmiles(mol)
            elif input_type.lower() == 'inchi_key':
                if enable_online_lookup:
                    return MolecularCalculator.convert_inchi_key_to_smiles(input_text)
                else:
                    return None  # Cannot convert without database lookup
        except Exception as e:
            return None

    @staticmethod
    def calculate_molecular_properties(smiles: str) -> Dict[str, Any]:
        """
        Calculate various molecular properties from SMILES

        Args:
            smiles: SMILES string

        Returns:
            Dictionary containing calculated properties
        """
        if not smiles or smiles is None:
            return {}

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {}

            # Check if molecule is valid (has atoms)
            if mol.GetNumAtoms() == 0:
                return {}

            properties = {}

            # Basic properties
            properties['Molecular_Weight'] = round(Descriptors.MolWt(mol), 3)
            properties['Heavy_Atom_Count'] = Descriptors.HeavyAtomCount(mol)
            properties['Atom_Count'] = mol.GetNumAtoms()
            properties['Bond_Count'] = mol.GetNumBonds()

            # Lipinski Rule of Five
            properties['LogP'] = round(Descriptors.MolLogP(mol), 3)
            properties['HB_Donors'] = Descriptors.NumHDonors(mol)
            properties['HB_Acceptors'] = Descriptors.NumHAcceptors(mol)
            properties['TPSA'] = round(Descriptors.TPSA(mol), 2)
            # Calculate 10×PSA/MW ratio (useful for predicting membrane permeability)
            properties['10xPSA_MW'] = round((10 * properties['TPSA']) / properties['Molecular_Weight'], 3) if properties['Molecular_Weight'] > 0 else 0
            properties['Rotatable_Bonds'] = Descriptors.NumRotatableBonds(mol)

            # Drug-likeness
            properties['QED'] = round(QED.qed(mol), 4)

            # Ring properties
            properties['Aromatic_Rings'] = Descriptors.NumAromaticRings(mol)
            properties['Aliphatic_Rings'] = Descriptors.NumAliphaticRings(mol)
            properties['Saturated_Rings'] = Descriptors.NumSaturatedRings(mol)
            properties['Heteroatoms'] = Descriptors.NumHeteroatoms(mol)

            # Crippen descriptors
            properties['CrippenLogP'] = round(Crippen.MolLogP(mol), 3)
            properties['CrippenMR'] = round(Crippen.MolMR(mol), 3)

            # Try to calculate additional descriptors with error handling
            try:
                properties['BertzCT'] = round(Descriptors.BertzCT(mol), 3)
            except:
                properties['BertzCT'] = None

            try:
                properties['LabuteASA'] = round(Descriptors.LabuteASA(mol), 3)
            except:
                properties['LabuteASA'] = None

            # Try rdMolDescriptors with fallback
            try:
                from rdkit.Chem import rdMolDescriptors
                properties['Chi0'] = round(rdMolDescriptors.Chi0(mol), 3)
                properties['Chi1'] = round(rdMolDescriptors.Chi1(mol), 3)
            except:
                properties['Chi0'] = None
                properties['Chi1'] = None

            # Additional simple descriptors
            try:
                properties['Formal_Charge'] = Chem.rdmolops.GetFormalCharge(mol)
            except:
                properties['Formal_Charge'] = 0

            # Fragment-based descriptors
            try:
                properties['Ring_Count'] = rdMolDescriptors.CalcNumRings(mol) if hasattr(rdMolDescriptors, 'CalcNumRings') else None
            except:
                properties['Ring_Count'] = None

            # Lipinski Rule compliance (0 = compliant, 1 = violates rule)
            lipinski_violations = sum([
                properties['Molecular_Weight'] > 500,
                properties['LogP'] > 5,
                properties['HB_Donors'] > 5,
                properties['HB_Acceptors'] > 10
            ])
            properties['Lipinski_Violations'] = 1 if lipinski_violations > 0 else 0

            # Veber Rule compliance (0 = compliant, 1 = violates rule)
            veber_violations = sum([
                properties['TPSA'] > 140,
                properties['Rotatable_Bonds'] > 10
            ])
            properties['Veber_Violations'] = 1 if veber_violations > 0 else 0

            # Filter out None values
            properties = {k: v for k, v in properties.items() if v is not None}

            return properties
        except Exception as e:
            return {}

    @staticmethod
    def detect_smiles_column(df: pd.DataFrame) -> Optional[str]:
        """
        Detect SMILES column with case-insensitive matching

        Args:
            df: Pandas DataFrame

        Returns:
            Column name containing SMILES or None if not found
        """
        possible_names = ['smiles', 'SMILES', 'Smiles', 'smi', 'SMI', 'canonical_smiles', 'CANONICAL_SMILES']
        for col in df.columns:
            if col in possible_names or col.lower() in [name.lower() for name in possible_names]:
                return col
        return None

    @staticmethod
    def detect_input_format(input_text: str) -> str:
        """
        Auto-detect input format

        Args:
            input_text: Input molecular structure string

        Returns:
            Detected format ('inchi', 'inchi_key', or 'smiles')
        """
        if input_text.startswith('InChI='):
            return 'inchi'
        elif re.match(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$', input_text):
            return 'inchi_key'
        else:
            return 'smiles'

    @staticmethod
    def get_property_groups(include_lei: bool = False) -> Dict[str, list]:
        """
        Get organized property groups for UI display

        Args:
            include_lei: Include Ligand Efficiency Indices (requires pKi)

        Returns:
            Dictionary with property groups and their properties
        """
        groups = {
            "Basic Properties": ['Molecular_Weight', 'Heavy_Atom_Count', 'Atom_Count', 'Bond_Count', 'Formal_Charge'],
            "Lipinski Properties": ['LogP', 'HB_Donors', 'HB_Acceptors', 'TPSA', '10xPSA_MW', 'Rotatable_Bonds'],
            "Drug-likeness": ['QED'],
            "Rule Violations": ['Lipinski_Violations', 'Veber_Violations'],
            "Ring Properties": ['Aromatic_Rings', 'Aliphatic_Rings', 'Saturated_Rings', 'Ring_Count', 'Heteroatoms'],
            "Complexity": ['BertzCT', 'Chi0', 'Chi1'],
            "Additional": ['CrippenLogP', 'CrippenMR', 'LabuteASA']
        }

        # Add LEI properties if requested
        if include_lei:
            lei_group = get_lei_property_group()
            groups.update(lei_group)

        return groups

    @classmethod
    def process_batch(cls, df: pd.DataFrame, smiles_col: str, selected_properties: set = None) -> pd.DataFrame:
        """
        Process batch of molecules and calculate properties

        Args:
            df: DataFrame containing molecules
            smiles_col: Name of column containing SMILES
            selected_properties: Set of properties to calculate (None for all)

        Returns:
            DataFrame with calculated properties
        """
        results = []

        for idx, row in df.iterrows():
            smiles = row[smiles_col]

            # Auto-detect and convert if needed
            input_format = cls.detect_input_format(str(smiles)) if pd.notna(smiles) else 'smiles'
            if input_format != 'smiles':
                smiles = cls.convert_to_smiles(str(smiles), input_format)

            properties = cls.calculate_molecular_properties(smiles) if pd.notna(smiles) else {}

            # Filter properties based on selection
            if properties and selected_properties:
                properties = {k: v for k, v in properties.items() if k in selected_properties}

            results.append(properties)

        # Create results DataFrame
        results_df = pd.DataFrame(results)
        final_df = pd.concat([df, results_df], axis=1)

        return final_df


class ThreeDOLSRegression:
    """
    3D Ordinary Least Squares (OLS) Regression Calculator

    Fits a plane to 3D data: Z = b0 + b1*X + b2*Y
    Based on analytical formulas from partial differentiation of the cost function
    """

    def __init__(self, x: np.ndarray, y: np.ndarray, z: np.ndarray, original_indices: np.ndarray = None):
        """
        Initialize 3D OLS regression with data

        Args:
            x: 1D array of X values
            y: 1D array of Y values
            z: 1D array of Z values (dependent variable)
            original_indices: Optional array of original indices (for tracking valid data points)
        """
        # Convert to numpy arrays and remove any NaN values
        self.x = np.array(x, dtype=float)
        self.y = np.array(y, dtype=float)
        self.z = np.array(z, dtype=float)
        
        # Store or create original indices
        if original_indices is not None:
            self.original_indices = np.array(original_indices)
        else:
            self.original_indices = np.arange(len(self.x))

        # Create mask for valid data (no NaN or Inf)
        valid_mask = (
            np.isfinite(self.x) &
            np.isfinite(self.y) &
            np.isfinite(self.z)
        )

        self.x = self.x[valid_mask]
        self.y = self.y[valid_mask]
        self.z = self.z[valid_mask]
        self.valid_indices = self.original_indices[valid_mask]  # Track which rows were kept

        self.n = len(self.x)

        if self.n < 3:
            raise ValueError("At least 3 valid data points are required for 3D OLS regression")

        # Calculate coefficients
        self.b0, self.b1, self.b2 = self._calculate_coefficients()

        # Calculate statistics
        self.r_squared = self._calculate_r_squared()
        self.residuals = self._calculate_residuals()
        self.rmse = self._calculate_rmse()

    def _calculate_coefficients(self) -> tuple:
        """
        Calculate OLS coefficients using analytical formulas

        Based on formulas from: https://sapiencespace.com/breaking-down-3d-linear-regression/

        Returns:
            Tuple of (b0, b1, b2) coefficients
        """
        # Calculate means
        x_mean = np.mean(self.x)
        y_mean = np.mean(self.y)
        z_mean = np.mean(self.z)

        # Calculate deviations from mean
        x_dev = self.x - x_mean
        y_dev = self.y - y_mean
        z_dev = self.z - z_mean

        # Calculate intermediate sums (Sa, Sb, Sc, Sd, Se)
        Sa = np.sum(x_dev ** 2)  # Σ(x - x_mean)²
        Sb = np.sum(y_dev ** 2)  # Σ(y - y_mean)²
        Sc = np.sum(x_dev * y_dev)  # Σ(x - x_mean)(y - y_mean)
        Sd = np.sum(x_dev * z_dev)  # Σ(x - x_mean)(z - z_mean)
        Se = np.sum(y_dev * z_dev)  # Σ(y - y_mean)(z - z_mean)

        # Calculate denominator for slope coefficients
        denominator = Sa * Sb - Sc ** 2

        if abs(denominator) < 1e-10:
            # Handle degenerate case (perfect multicollinearity)
            raise ValueError("X and Y variables are perfectly collinear")

        # Calculate slope coefficients (from normal equations inversion)
        # b1 = (Sb*Sd - Sc*Se)/D ; b2 = (Sa*Se - Sc*Sd)/D
        b1 = (Sd * Sb - Sc * Se) / denominator
        b2 = (Sa * Se - Sc * Sd) / denominator

        # Calculate intercept
        b0 = z_mean - b1 * x_mean - b2 * y_mean

        return b0, b1, b2

    def _calculate_residuals(self) -> np.ndarray:
        """Calculate residuals (observed - predicted)"""
        z_predicted = self.predict(self.x, self.y)
        return self.z - z_predicted

    def _calculate_r_squared(self) -> float:
        """
        Calculate R² (coefficient of determination)

        R² = 1 - (SS_res / SS_tot)
        where SS_res = Σ(observed - predicted)²
              SS_tot = Σ(observed - mean)²

        Returns:
            R² value between -∞ and 1 (closer to 1 is better)
        """
        z_predicted = self.predict(self.x, self.y)
        z_mean = np.mean(self.z)

        # Sum of squares of residuals
        ss_res = np.sum((self.z - z_predicted) ** 2)

        # Total sum of squares
        ss_tot = np.sum((self.z - z_mean) ** 2)

        if ss_tot < 1e-10:
            # Handle case where all Z values are the same
            return 1.0 if ss_res < 1e-10 else 0.0

        r_squared = 1 - (ss_res / ss_tot)

        return r_squared

    def _calculate_rmse(self) -> float:
        """
        Calculate Root Mean Squared Error

        RMSE = sqrt(Σ(observed - predicted)² / n)

        Returns:
            RMSE value
        """
        return np.sqrt(np.mean(self.residuals ** 2))

    def predict(self, x: Union[float, np.ndarray], y: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Predict Z values using the fitted plane equation

        Args:
            x: X value(s) for prediction
            y: Y value(s) for prediction

        Returns:
            Predicted Z value(s)
        """
        return self.b0 + self.b1 * np.array(x) + self.b2 * np.array(y)

    def get_equation_string(self, decimals: int = 4) -> str:
        """
        Get the plane equation as a formatted string

        Args:
            decimals: Number of decimal places for coefficients

        Returns:
            Equation string in format: Z = b0 + b1*X + b2*Y
        """
        b0_str = f"{self.b0:.{decimals}f}"
        b1_str = f"{abs(self.b1):.{decimals}f}"
        b2_str = f"{abs(self.b2):.{decimals}f}"

        # Format with proper signs
        b1_sign = "+" if self.b1 >= 0 else "-"
        b2_sign = "+" if self.b2 >= 0 else "-"

        equation = f"Z = {b0_str} {b1_sign} {b1_str}·X {b2_sign} {b2_str}·Y"

        return equation

    def get_statistics(self) -> Dict[str, float]:
        """
        Get regression statistics

        Returns:
            Dictionary containing regression statistics
        """
        return {
            'b0': self.b0,
            'b1': self.b1,
            'b2': self.b2,
            'R²': self.r_squared,
            'RMSE': self.rmse,
            'n': self.n
        }

    def get_plane_mesh(self, num_points: int = 20) -> tuple:
        """
        Generate mesh grid for plotting the fitted plane

        Args:
            num_points: Number of points in each dimension for the mesh

        Returns:
            Tuple of (X_mesh, Y_mesh, Z_mesh) for 3D plotting
        """
        # Create grid based on data range with some padding
        x_min, x_max = self.x.min(), self.x.max()
        y_min, y_max = self.y.min(), self.y.max()

        x_range = x_max - x_min
        y_range = y_max - y_min

        # Add 10% padding
        x_min -= 0.1 * x_range
        x_max += 0.1 * x_range
        y_min -= 0.1 * y_range
        y_max += 0.1 * y_range

        # Create mesh grid
        X_mesh = np.linspace(x_min, x_max, num_points)
        Y_mesh = np.linspace(y_min, y_max, num_points)
        X_mesh, Y_mesh = np.meshgrid(X_mesh, Y_mesh)

        # Calculate Z values on the mesh
        Z_mesh = self.predict(X_mesh, Y_mesh)

        return X_mesh, Y_mesh, Z_mesh


class PropertyExplanations:
    """Class containing property explanations and documentation"""

    @staticmethod
    def get_explanations() -> str:
        """
        Get detailed explanations of all molecular properties

        Returns:
            Markdown formatted string with property explanations
        """
        return """
        ### Supported Input Formats
        - **SMILES**: Simplified Molecular Input Line Entry System
        - **InChI**: International Chemical Identifier
        - **InChI Key**: Hashed version of InChI (database lookup required)

        ### File Formats
        - **CSV**: Comma-separated values
        - **XLSX**: Excel spreadsheet format

        ### Property Explanations

        #### Basic Properties
        - **Molecular_Weight**: Molecular weight in Daltons (Da)
        - **Heavy_Atom_Count**: Number of non-hydrogen atoms
        - **Atom_Count**: Total number of atoms (including hydrogens)
        - **Bond_Count**: Total number of bonds
        - **Formal_Charge**: Net formal charge of the molecule

        #### Lipinski Properties (Rule of Five)
        - **LogP**: Partition coefficient (lipophilicity, -2 to 6 typical range)
        - **HB_Donors**: Hydrogen bond donors (≤5 for drug-likeness)
        - **HB_Acceptors**: Hydrogen bond acceptors (≤10 for drug-likeness)
        - **TPSA**: Topological polar surface area in Ų (≤140 for oral bioavailability)
        - **10xPSA_MW**: PSA/MW ratio scaled by 10 (lower values indicate better membrane permeability)
        - **Rotatable_Bonds**: Number of rotatable bonds (≤10 for oral bioavailability)

        #### Drug-likeness
        - **QED**: Quantitative Estimate of Drug-likeness (0-1 scale, higher is better)

        #### Rule Violations (Binary: 0=Compliant, 1=Violates)
        - **Lipinski_Violations**:
          - 0 = Passes Lipinski Rule of Five (MW≤500, LogP≤5, HBD≤5, HBA≤10)
          - 1 = Violates at least one Lipinski rule
        - **Veber_Violations**:
          - 0 = Passes Veber Rule (TPSA≤140, RotBonds≤10)
          - 1 = Violates at least one Veber rule

        #### Ring Properties
        - **Aromatic_Rings**: Number of aromatic rings
        - **Aliphatic_Rings**: Number of non-aromatic rings
        - **Saturated_Rings**: Number of saturated rings
        - **Ring_Count**: Total number of rings
        - **Heteroatoms**: Number of non-carbon, non-hydrogen atoms

        #### Complexity
        - **BertzCT**: Bertz complexity index (higher = more complex)
        - **Chi0**: Chi connectivity index 0 (molecular connectivity)
        - **Chi1**: Chi connectivity index 1 (molecular connectivity)

        #### Additional Descriptors
        - **CrippenLogP**: Crippen's LogP calculation method
        - **CrippenMR**: Crippen's molar refractivity
        - **LabuteASA**: Labute's approximate surface area

        ### Usage Tips
        - No properties are selected by default - choose what you need
        - Perfect for files with existing calculations - add only missing properties
        - Common SMILES column names are automatically detected
        - Invalid molecules will result in empty property values
        - All calculations use RDKit library
        - Results are compatible with StarDrop, Pipeline Pilot, and other software

        ### InChI Key Conversion
        - **Enabled by default** for convenience
        - Uses **NIH Chemical Identifier Resolver** (primary)
        - **PubChem API** fallback for redundancy
        - Can be disabled in Settings for privacy/offline use
        - Requires internet connection for database lookup

        ### Ligand Efficiency Indices (LEI) - Requires pKi Values
        Based on AtlasCBS methodology (Cele Abad-Zapatero, A. Cortes-Cabrera 2013)

        #### Available LEIs:
        - **NSEI**: Normalized Surface Efficiency Index (pKi / Polar Atoms)
        - **NBEI**: Normalized Binding Efficiency Index (pKi / Heavy Atoms)
        - **BEI**: Binding Efficiency Index (pKi / (MW/1000))
        - **SEI**: Surface Efficiency Index (pKi / (TPSA/100))
        - **nBEI**: Alternative Binding Efficiency Index (-log10(Ki / Heavy Atoms))
        - **mBEI**: Molecular Binding Efficiency Index (-log10(Ki / MW))
        - **LEH**: Ligand Efficiency Hopkins (-ΔG / Heavy Atoms) [Hopkins et al. DDT, 2004]
        - **LEP**: Ligand Efficiency Polar (-ΔG / Polar Atoms)

        #### LEI Requirements:
        - **pKi column required** in uploaded file
        - Optional columns (will be calculated from SMILES if missing):
          - MW (Molecular Weight)
          - TPSA (Topological Polar Surface Area)
          - Heavy_Atom_Count
        - System automatically detects available columns
        - Only calculates missing values from SMILES when needed

        ### Development Information
        Developed by: **Yashwanth Reddy** for **ITR-UIC**
        Part of: **Chemo-Informatics Toolkit**
        """
