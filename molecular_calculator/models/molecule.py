"""Data models for molecular structures and properties.

This module defines data classes for representing molecules,
their properties, and calculation results.
"""

from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List
from enum import Enum


class InputFormat(Enum):
    """Enumeration of supported molecular input formats."""
    SMILES = "smiles"
    INCHI = "inchi"
    INCHI_KEY = "inchi_key"
    UNKNOWN = "unknown"


class PropertyGroup(Enum):
    """Enumeration of property categories."""
    BASIC = "Basic Properties"
    LIPINSKI = "Lipinski Properties"
    DRUGLIKENESS = "Drug-likeness"
    RULES = "Rule Violations"
    RINGS = "Ring Properties"
    COMPLEXITY = "Complexity"
    ADDITIONAL = "Additional"
    LEI = "Ligand Efficiency Indices"
    ASSAY_INTERFERENCE = "Assay Interference"


@dataclass
class MoleculeInput:
    """Represents a molecular structure input.

    Attributes:
        value: The input string (SMILES, InChI, or InChI Key)
        format: Detected or specified input format
        name: Optional molecule name/identifier

    Example:
        >>> mol_input = MoleculeInput(value="CCO", format=InputFormat.SMILES, name="Ethanol")
    """
    value: str
    format: InputFormat = InputFormat.UNKNOWN
    name: Optional[str] = None

    def __post_init__(self):
        """Validate and normalize the input."""
        # Always validate - empty strings should also be rejected
        self.value = self.value.strip()
        if not self.value:
            raise ValueError("Input value cannot be empty or whitespace-only")


@dataclass
class MolecularProperties:
    """Stores calculated molecular properties.

    All properties are optional as not all may be calculable for every molecule.

    Attributes:
        molecular_weight: Molecular weight in Daltons
        heavy_atom_count: Number of non-hydrogen atoms
        atom_count: Total number of atoms
        bond_count: Total number of bonds
        logp: Partition coefficient (lipophilicity)
        hb_donors: Hydrogen bond donors
        hb_acceptors: Hydrogen bond acceptors
        tpsa: Topological polar surface area
        psa_mw_ratio: 10×PSA/MW ratio
        npol_nha: Polar atoms / Heavy atoms ratio
        rotatable_bonds: Number of rotatable bonds
        qed: Quantitative Estimate of Drug-likeness
        aromatic_rings: Number of aromatic rings
        aliphatic_rings: Number of non-aromatic rings
        saturated_rings: Number of saturated rings
        ring_count: Total number of rings
        heteroatoms: Number of heteroatoms
        formal_charge: Net formal charge
        crippen_logp: Crippen's LogP calculation
        crippen_mr: Crippen's molar refractivity
        bertz_ct: Bertz complexity index
        labute_asa: Labute's approximate surface area
        chi0: Chi connectivity index 0
        chi1: Chi connectivity index 1
        lipinski_violations: 0 if compliant, 1 if violates
        veber_violations: 0 if compliant, 1 if violates
    """
    # Basic properties
    molecular_weight: Optional[float] = None
    heavy_atom_count: Optional[int] = None
    atom_count: Optional[int] = None
    bond_count: Optional[int] = None
    formal_charge: Optional[int] = None

    # Lipinski properties
    logp: Optional[float] = None
    hb_donors: Optional[int] = None
    hb_acceptors: Optional[int] = None
    tpsa: Optional[float] = None
    psa_mw_ratio: Optional[float] = None
    npol_nha: Optional[float] = None
    rotatable_bonds: Optional[int] = None

    # Drug-likeness
    qed: Optional[float] = None

    # Ring properties
    aromatic_rings: Optional[int] = None
    aliphatic_rings: Optional[int] = None
    saturated_rings: Optional[int] = None
    ring_count: Optional[int] = None
    heteroatoms: Optional[int] = None

    # Complexity
    bertz_ct: Optional[float] = None
    chi0: Optional[float] = None
    chi1: Optional[float] = None

    # Additional descriptors
    crippen_logp: Optional[float] = None
    crippen_mr: Optional[float] = None
    labute_asa: Optional[float] = None

    # Rule violations (binary)
    lipinski_violations: Optional[int] = None
    veber_violations: Optional[int] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert properties to dictionary format.

        Returns:
            Dictionary with property names as keys, filtering out None values.
        """
        # Map attribute names to display names
        name_mapping = {
            'molecular_weight': 'Molecular_Weight',
            'heavy_atom_count': 'Heavy_Atom_Count',
            'atom_count': 'Atom_Count',
            'bond_count': 'Bond_Count',
            'formal_charge': 'Formal_Charge',
            'logp': 'LogP',
            'hb_donors': 'HB_Donors',
            'hb_acceptors': 'HB_Acceptors',
            'tpsa': 'TPSA',
            'psa_mw_ratio': '10xPSA_MW',
            'npol_nha': 'NPOLoNHA',
            'rotatable_bonds': 'Rotatable_Bonds',
            'qed': 'QED',
            'aromatic_rings': 'Aromatic_Rings',
            'aliphatic_rings': 'Aliphatic_Rings',
            'saturated_rings': 'Saturated_Rings',
            'ring_count': 'Ring_Count',
            'heteroatoms': 'Heteroatoms',
            'bertz_ct': 'BertzCT',
            'chi0': 'Chi0',
            'chi1': 'Chi1',
            'crippen_logp': 'CrippenLogP',
            'crippen_mr': 'CrippenMR',
            'labute_asa': 'LabuteASA',
            'lipinski_violations': 'Lipinski_Violations',
            'veber_violations': 'Veber_Violations',
        }

        result = {}
        for attr, display_name in name_mapping.items():
            value = getattr(self, attr)
            if value is not None:
                result[display_name] = value
        return result

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'MolecularProperties':
        """Create MolecularProperties from a dictionary.

        Args:
            data: Dictionary with property names (can use display names or attribute names)

        Returns:
            MolecularProperties instance
        """
        # Reverse mapping from display names to attribute names
        reverse_mapping = {
            'Molecular_Weight': 'molecular_weight',
            'Heavy_Atom_Count': 'heavy_atom_count',
            'Atom_Count': 'atom_count',
            'Bond_Count': 'bond_count',
            'Formal_Charge': 'formal_charge',
            'LogP': 'logp',
            'HB_Donors': 'hb_donors',
            'HB_Acceptors': 'hb_acceptors',
            'TPSA': 'tpsa',
            '10xPSA_MW': 'psa_mw_ratio',
            'NPOLoNHA': 'npol_nha',
            'Rotatable_Bonds': 'rotatable_bonds',
            'QED': 'qed',
            'Aromatic_Rings': 'aromatic_rings',
            'Aliphatic_Rings': 'aliphatic_rings',
            'Saturated_Rings': 'saturated_rings',
            'Ring_Count': 'ring_count',
            'Heteroatoms': 'heteroatoms',
            'BertzCT': 'bertz_ct',
            'Chi0': 'chi0',
            'Chi1': 'chi1',
            'CrippenLogP': 'crippen_logp',
            'CrippenMR': 'crippen_mr',
            'LabuteASA': 'labute_asa',
            'Lipinski_Violations': 'lipinski_violations',
            'Veber_Violations': 'veber_violations',
        }

        kwargs = {}
        for key, value in data.items():
            # Try display name first, then attribute name
            attr_name = reverse_mapping.get(key, key)
            if hasattr(cls, '__dataclass_fields__') and attr_name in cls.__dataclass_fields__:
                kwargs[attr_name] = value

        return cls(**kwargs)


@dataclass
class CalculationResult:
    """Result of a molecular property calculation.

    Attributes:
        success: Whether the calculation succeeded
        smiles: The SMILES string used for calculation
        properties: Calculated molecular properties
        error: Error message if calculation failed
        warnings: List of warning messages

    Example:
        >>> result = CalculationResult(success=True, smiles="CCO", properties=props)
        >>> if result.success:
        ...     print(result.properties.molecular_weight)
    """
    success: bool
    smiles: Optional[str] = None
    properties: Optional[MolecularProperties] = None
    error: Optional[str] = None
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary format."""
        if self.success and self.properties:
            return self.properties.to_dict()
        return {}


@dataclass
class ConversionResult:
    """Result of a format conversion operation.

    Attributes:
        success: Whether the conversion succeeded
        smiles: The resulting SMILES string
        source_format: Original input format
        error: Error message if conversion failed
    """
    success: bool
    smiles: Optional[str] = None
    source_format: Optional[InputFormat] = None
    error: Optional[str] = None


@dataclass
class LigandEfficiencyIndices:
    """Stores calculated Ligand Efficiency Indices.

    Based on AtlasCBS methodology (Cele Abad-Zapatero, A. Cortes-Cabrera 2013).

    Attributes:
        nsei: Normalized Surface Efficiency Index (pKi / Polar Atoms)
        nbei: Normalized Binding Efficiency Index (pKi / Heavy Atoms)
        bei: Binding Efficiency Index (pKi / (MW/1000))
        sei: Surface Efficiency Index (pKi / (TPSA/100))
        nbei_alt: Alternative Binding Efficiency (-log10(Ki / Heavy Atoms))
        mbei: Molecular Binding Efficiency (-log10(Ki / MW))
        leh: Ligand Efficiency Hopkins (-ΔG / Heavy Atoms)
        lep: Ligand Efficiency Polar (-ΔG / Polar Atoms)
    """
    nsei: Optional[float] = None
    nbei: Optional[float] = None
    bei: Optional[float] = None
    sei: Optional[float] = None
    nbei_alt: Optional[float] = None  # nBEI in original
    mbei: Optional[float] = None
    leh: Optional[float] = None
    lep: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format."""
        name_mapping = {
            'nsei': 'NSEI',
            'nbei': 'NBEI',
            'bei': 'BEI',
            'sei': 'SEI',
            'nbei_alt': 'nBEI',
            'mbei': 'mBEI',
            'leh': 'LEH',
            'lep': 'LEP',
        }

        result = {}
        for attr, display_name in name_mapping.items():
            value = getattr(self, attr)
            if value is not None:
                result[display_name] = value
        return result


# Property group definitions for UI display and selection.
# NOTE: Assay interference flags (PAINS, Aggregator, etc.) are calculated separately
# via get_interference_flags_from_smiles() in assay_interference.py, not through
# MolecularProperties.to_dict(). They are listed here for UI property selection purposes.
PROPERTY_GROUPS: Dict[str, List[str]] = {
    PropertyGroup.BASIC.value: [
        'Molecular_Weight', 'Heavy_Atom_Count', 'Atom_Count',
        'Bond_Count', 'Formal_Charge'
    ],
    PropertyGroup.LIPINSKI.value: [
        'LogP', 'HB_Donors', 'HB_Acceptors', 'TPSA',
        '10xPSA_MW', 'NPOLoNHA', 'Rotatable_Bonds'
    ],
    PropertyGroup.DRUGLIKENESS.value: ['QED'],
    PropertyGroup.RULES.value: ['Lipinski_Violations', 'Veber_Violations'],
    PropertyGroup.RINGS.value: [
        'Aromatic_Rings', 'Aliphatic_Rings', 'Saturated_Rings',
        'Ring_Count', 'Heteroatoms'
    ],
    PropertyGroup.COMPLEXITY.value: ['BertzCT', 'Chi0', 'Chi1'],
    PropertyGroup.ADDITIONAL.value: ['CrippenLogP', 'CrippenMR', 'LabuteASA'],
    PropertyGroup.ASSAY_INTERFERENCE.value: [
        'PAINS', 'Aggregator', 'Redox', 'Fluorescence', 'Thiol'
    ],
}

LEI_PROPERTY_GROUP: Dict[str, List[str]] = {
    PropertyGroup.LEI.value: [
        'NSEI', 'NBEI', 'BEI', 'SEI',
        'nBEI', 'mBEI', 'LEH', 'LEP'
    ]
}
