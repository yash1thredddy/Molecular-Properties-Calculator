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
    def get_property_groups() -> Dict[str, list]:
        """
        Get organized property groups for UI display

        Returns:
            Dictionary with property groups and their properties
        """
        return {
            "Basic Properties": ['Molecular_Weight', 'Heavy_Atom_Count', 'Atom_Count', 'Bond_Count', 'Formal_Charge'],
            "Lipinski Properties": ['LogP', 'HB_Donors', 'HB_Acceptors', 'TPSA', 'Rotatable_Bonds'],
            "Drug-likeness": ['QED'],
            "Rule Violations": ['Lipinski_Violations', 'Veber_Violations'],
            "Ring Properties": ['Aromatic_Rings', 'Aliphatic_Rings', 'Saturated_Rings', 'Ring_Count', 'Heteroatoms'],
            "Complexity": ['BertzCT', 'Chi0', 'Chi1'],
            "Additional": ['CrippenLogP', 'CrippenMR', 'LabuteASA']
        }

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

        ### Development Information
        Developed by: **Yashwanth Reddy** for **ITR-UIC**
        Part of: **Chemo-Informatics Toolkit**
        """