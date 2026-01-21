"""Molecular property calculation service.

This module provides the core property calculation functionality
using RDKit descriptors.
"""

import logging
from typing import Dict, Any, Optional, List, Set

from rdkit import Chem, RDLogger
from rdkit.Chem import QED, Descriptors, Crippen, rdMolDescriptors

from molecular_calculator.models import (
    MolecularProperties,
    CalculationResult,
    PROPERTY_GROUPS,
)
from molecular_calculator.utils.exceptions import (
    InvalidSMILESError,
    PropertyCalculationError,
)

logger = logging.getLogger(__name__)


class PropertyCalculator:
    """Service for calculating molecular properties from SMILES.

    Uses RDKit to compute various molecular descriptors including
    physicochemical properties, drug-likeness metrics, and rule violations.

    Example:
        >>> calculator = PropertyCalculator()
        >>> result = calculator.calculate("CCO")
        >>> if result.success:
        ...     print(f"MW: {result.properties.molecular_weight}")
    """

    def __init__(self, suppress_warnings: bool = True):
        """Initialize the property calculator.

        Args:
            suppress_warnings: Whether to suppress RDKit warning messages
        """
        if suppress_warnings:
            RDLogger.DisableLog('rdApp.*')

    @staticmethod
    def enable_rdkit_warnings():
        """Re-enable RDKit warning messages."""
        RDLogger.EnableLog('rdApp.*')

    @staticmethod
    def disable_rdkit_warnings():
        """Disable RDKit warning messages."""
        RDLogger.DisableLog('rdApp.*')

    def _parse_molecule(self, smiles: str) -> Optional[Chem.Mol]:
        """Parse SMILES string into RDKit molecule.

        Args:
            smiles: SMILES string to parse

        Returns:
            RDKit Mol object or None if parsing fails
        """
        if not smiles or not isinstance(smiles, str):
            return None

        try:
            mol = Chem.MolFromSmiles(smiles.strip())
            if mol is not None and mol.GetNumAtoms() > 0:
                return mol
        except Exception as e:
            logger.debug(f"SMILES parsing failed: {e}")

        return None

    def _calculate_basic_properties(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Calculate basic molecular properties.

        Args:
            mol: RDKit Mol object

        Returns:
            Dictionary of basic properties
        """
        return {
            'molecular_weight': round(Descriptors.MolWt(mol), 3),
            'heavy_atom_count': Descriptors.HeavyAtomCount(mol),
            'atom_count': mol.GetNumAtoms(),
            'bond_count': mol.GetNumBonds(),
        }

    def _calculate_lipinski_properties(
        self,
        mol: Chem.Mol,
        molecular_weight: float,
        heavy_atom_count: int
    ) -> Dict[str, Any]:
        """Calculate Lipinski rule-related properties.

        Args:
            mol: RDKit Mol object
            molecular_weight: Pre-calculated molecular weight
            heavy_atom_count: Pre-calculated heavy atom count

        Returns:
            Dictionary of Lipinski properties
        """
        tpsa = round(Descriptors.TPSA(mol), 2)

        # Calculate 10×PSA/MW ratio
        psa_mw_ratio = 0.0
        if molecular_weight > 0:
            psa_mw_ratio = round((10 * tpsa) / molecular_weight, 3)

        # Calculate NPOLoNHA (Polar Atoms / Heavy Atoms)
        # NPOL = Number of heteroatoms (N, O, S, P, halogens, etc.)
        npol = Descriptors.NumHeteroatoms(mol)
        npol_nha = 0.0
        if heavy_atom_count > 0:
            npol_nha = round(npol / heavy_atom_count, 3)

        return {
            'logp': round(Descriptors.MolLogP(mol), 3),
            'hb_donors': Descriptors.NumHDonors(mol),
            'hb_acceptors': Descriptors.NumHAcceptors(mol),
            'tpsa': tpsa,
            'psa_mw_ratio': psa_mw_ratio,
            'npol_nha': npol_nha,
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
        }

    def _calculate_ring_properties(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Calculate ring-related properties.

        Args:
            mol: RDKit Mol object

        Returns:
            Dictionary of ring properties
        """
        props = {
            'aromatic_rings': Descriptors.NumAromaticRings(mol),
            'aliphatic_rings': Descriptors.NumAliphaticRings(mol),
            'saturated_rings': Descriptors.NumSaturatedRings(mol),
            'heteroatoms': Descriptors.NumHeteroatoms(mol),
        }

        # Ring count with fallback
        try:
            props['ring_count'] = rdMolDescriptors.CalcNumRings(mol)
        except Exception:
            props['ring_count'] = None

        return props

    def _calculate_complexity_properties(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Calculate molecular complexity properties.

        Args:
            mol: RDKit Mol object

        Returns:
            Dictionary of complexity properties
        """
        props = {}

        # BertzCT
        try:
            props['bertz_ct'] = round(Descriptors.BertzCT(mol), 3)
        except Exception:
            props['bertz_ct'] = None

        # Chi connectivity indices
        try:
            props['chi0'] = round(rdMolDescriptors.Chi0(mol), 3)
            props['chi1'] = round(rdMolDescriptors.Chi1(mol), 3)
        except Exception:
            props['chi0'] = None
            props['chi1'] = None

        return props

    def _calculate_additional_properties(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Calculate additional molecular properties.

        Args:
            mol: RDKit Mol object

        Returns:
            Dictionary of additional properties
        """
        props = {
            'crippen_logp': round(Crippen.MolLogP(mol), 3),
            'crippen_mr': round(Crippen.MolMR(mol), 3),
        }

        # LabuteASA with fallback
        try:
            props['labute_asa'] = round(Descriptors.LabuteASA(mol), 3)
        except Exception:
            props['labute_asa'] = None

        # Formal charge
        try:
            props['formal_charge'] = Chem.rdmolops.GetFormalCharge(mol)
        except Exception:
            props['formal_charge'] = 0

        return props

    def _calculate_rule_violations(
        self,
        molecular_weight: float,
        logp: float,
        hb_donors: int,
        hb_acceptors: int,
        tpsa: float,
        rotatable_bonds: int
    ) -> Dict[str, int]:
        """Calculate drug-likeness rule violations.

        Args:
            molecular_weight: Molecular weight
            logp: LogP value
            hb_donors: Number of H-bond donors
            hb_acceptors: Number of H-bond acceptors
            tpsa: Topological polar surface area
            rotatable_bonds: Number of rotatable bonds

        Returns:
            Dictionary with violation flags (0 = compliant, 1 = violates)
        """
        # Lipinski Rule of Five
        lipinski_violations = sum([
            molecular_weight > 500,
            logp > 5,
            hb_donors > 5,
            hb_acceptors > 10
        ])

        # Veber Rule
        veber_violations = sum([
            tpsa > 140,
            rotatable_bonds > 10
        ])

        return {
            'lipinski_violations': 1 if lipinski_violations > 0 else 0,
            'veber_violations': 1 if veber_violations > 0 else 0,
        }

    def calculate(
        self,
        smiles: str,
        selected_properties: Set[str] = None  # Reserved for future optimization
    ) -> CalculationResult:
        """Calculate molecular properties from SMILES.

        Args:
            smiles: SMILES string
            selected_properties: Reserved for future use. Currently all properties
                                are calculated regardless of this parameter.
                                Use calculate_as_dict() with selected_properties
                                to filter the returned properties.

        Returns:
            CalculationResult with calculated properties

        Example:
            >>> calc = PropertyCalculator()
            >>> result = calc.calculate("CCO")
            >>> print(result.properties.molecular_weight)
            46.069

        Note:
            All properties are always calculated internally due to
            interdependencies. Use calculate_as_dict() with selected_properties
            to get a filtered subset of properties.
        """
        if not smiles:
            return CalculationResult(
                success=False,
                error="SMILES string is empty"
            )

        # Parse molecule
        mol = self._parse_molecule(smiles)
        if mol is None:
            return CalculationResult(
                success=False,
                smiles=smiles,
                error=f"Invalid SMILES: {smiles}"
            )

        try:
            # Calculate all property groups
            basic = self._calculate_basic_properties(mol)
            lipinski = self._calculate_lipinski_properties(
                mol,
                basic['molecular_weight'],
                basic['heavy_atom_count']
            )
            ring = self._calculate_ring_properties(mol)
            complexity = self._calculate_complexity_properties(mol)
            additional = self._calculate_additional_properties(mol)

            # Calculate QED
            try:
                qed_value = round(QED.qed(mol), 4)
            except Exception:
                qed_value = None

            # Calculate rule violations
            violations = self._calculate_rule_violations(
                basic['molecular_weight'],
                lipinski['logp'],
                lipinski['hb_donors'],
                lipinski['hb_acceptors'],
                lipinski['tpsa'],
                lipinski['rotatable_bonds']
            )

            # Build MolecularProperties object
            props = MolecularProperties(
                molecular_weight=basic['molecular_weight'],
                heavy_atom_count=basic['heavy_atom_count'],
                atom_count=basic['atom_count'],
                bond_count=basic['bond_count'],
                formal_charge=additional.get('formal_charge', 0),
                logp=lipinski['logp'],
                hb_donors=lipinski['hb_donors'],
                hb_acceptors=lipinski['hb_acceptors'],
                tpsa=lipinski['tpsa'],
                psa_mw_ratio=lipinski['psa_mw_ratio'],
                npol_nha=lipinski['npol_nha'],
                rotatable_bonds=lipinski['rotatable_bonds'],
                qed=qed_value,
                aromatic_rings=ring['aromatic_rings'],
                aliphatic_rings=ring['aliphatic_rings'],
                saturated_rings=ring['saturated_rings'],
                ring_count=ring.get('ring_count'),
                heteroatoms=ring['heteroatoms'],
                bertz_ct=complexity.get('bertz_ct'),
                chi0=complexity.get('chi0'),
                chi1=complexity.get('chi1'),
                crippen_logp=additional['crippen_logp'],
                crippen_mr=additional['crippen_mr'],
                labute_asa=additional.get('labute_asa'),
                lipinski_violations=violations['lipinski_violations'],
                veber_violations=violations['veber_violations'],
            )

            return CalculationResult(
                success=True,
                smiles=smiles,
                properties=props
            )

        except Exception as e:
            logger.error(f"Property calculation failed for {smiles}: {e}")
            return CalculationResult(
                success=False,
                smiles=smiles,
                error=str(e)
            )

    def calculate_as_dict(
        self,
        smiles: str,
        selected_properties: Set[str] = None
    ) -> Dict[str, Any]:
        """Calculate properties and return as dictionary.

        Convenience method that returns properties directly as a dict.

        Args:
            smiles: SMILES string
            selected_properties: Optional set of property names to filter

        Returns:
            Dictionary of calculated properties (empty if failed)
        """
        result = self.calculate(smiles, selected_properties)

        if result.success and result.properties:
            props_dict = result.properties.to_dict()

            # Filter if specific properties requested
            if selected_properties:
                props_dict = {
                    k: v for k, v in props_dict.items()
                    if k in selected_properties
                }

            return props_dict

        return {}

    @staticmethod
    def get_property_groups(include_lei: bool = False) -> Dict[str, List[str]]:
        """Get organized property groups for UI display.

        Args:
            include_lei: Include Ligand Efficiency Indices (requires pKi)

        Returns:
            Dictionary with property groups and their properties
        """
        groups = dict(PROPERTY_GROUPS)

        if include_lei:
            from molecular_calculator.models import LEI_PROPERTY_GROUP
            groups.update(LEI_PROPERTY_GROUP)

        return groups

    @staticmethod
    def get_all_property_names() -> List[str]:
        """Get list of all available property names.

        Returns:
            List of all property names across all groups
        """
        all_props = []
        for props in PROPERTY_GROUPS.values():
            all_props.extend(props)
        return all_props


# Singleton instance for convenience
_property_calculator: Optional[PropertyCalculator] = None
_property_calculator_lock = __import__('threading').Lock()


def get_property_calculator() -> PropertyCalculator:
    """Get the shared property calculator instance.

    Thread-safe singleton accessor using double-checked locking.

    Returns:
        PropertyCalculator singleton instance
    """
    global _property_calculator
    if _property_calculator is None:
        with _property_calculator_lock:
            # Double-check after acquiring lock
            if _property_calculator is None:
                _property_calculator = PropertyCalculator()
    return _property_calculator
