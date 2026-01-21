"""
Assay Interference Detection Module

This module detects compounds with known assay interference mechanisms identified
in the 2016 Bisson et al. "Invalid Metabolic Panaceas (IMPs)" paper.

Five core interference mechanisms are detected:
1. PAINS (Pan-Assay Interference Substructures)
2. Aggregation risk (colloidal aggregators)
3. Redox reactivity (redox-active functional groups)
4. Fluorescence interference (autofluorescent compounds)
5. Thiol reactivity (cysteine-reactive electrophiles)

RDKit Dependencies:
- FilterCatalog: Built-in PAINS filter catalog for substructure matching
- rdMolDescriptors: Aromatic ring count, rotatable bonds, exact MW
- Descriptors: MolLogP, NumHeteroatoms
- Chem.MolFromSmarts: SMARTS pattern matching for custom filters

References:
- Bisson et al. (2016) J. Med. Chem. 59, 1671-1690 (Invalid Metabolic Panaceas)
- Baell & Holloway (2010) J. Med. Chem. 53, 2719-2740 (PAINS filters)
- Shoichet laboratory aggregator research (http://www.bkslab.org)
"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import FilterCatalog, rdMolDescriptors, Descriptors

logger = logging.getLogger(__name__)


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class InterferenceFlags:
    """Container for all assay interference flags.

    Attributes:
        pains: PAINS (Pan-Assay Interference Substructures) flag
        aggregator: Aggregation risk flag
        redox: Redox reactivity flag
        fluorescence: Fluorescence interference flag
        thiol: Thiol reactivity flag
        pains_details: List of specific PAINS violations
        aggregator_reason: Reason for aggregator flag
        redox_groups: List of detected redox-active groups
        fluorescence_scaffolds: List of detected fluorescent scaffolds
        thiol_electrophiles: List of detected thiol-reactive groups
    """
    pains: bool = False
    aggregator: bool = False
    redox: bool = False
    fluorescence: bool = False
    thiol: bool = False
    pains_details: List[str] = None
    aggregator_reason: str = ""
    redox_groups: List[str] = None
    fluorescence_scaffolds: List[str] = None
    thiol_electrophiles: List[str] = None

    def __post_init__(self):
        """Initialize empty lists if None."""
        if self.pains_details is None:
            self.pains_details = []
        if self.redox_groups is None:
            self.redox_groups = []
        if self.fluorescence_scaffolds is None:
            self.fluorescence_scaffolds = []
        if self.thiol_electrophiles is None:
            self.thiol_electrophiles = []

    @property
    def total_flags(self) -> int:
        """Count total number of flags raised."""
        return sum([
            self.pains,
            self.aggregator,
            self.redox,
            self.fluorescence,
            self.thiol
        ])

    @property
    def is_clean(self) -> bool:
        """Check if compound has no interference flags."""
        return self.total_flags == 0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for DataFrame/display."""
        return {
            'PAINS': 1 if self.pains else 0,
            'Aggregator': 1 if self.aggregator else 0,
            'Redox': 1 if self.redox else 0,
            'Fluorescence': 1 if self.fluorescence else 0,
            'Thiol': 1 if self.thiol else 0,
        }

    def to_detailed_dict(self) -> Dict[str, Any]:
        """Convert to detailed dictionary including reasons."""
        return {
            'PAINS': 1 if self.pains else 0,
            'PAINS_Details': ', '.join(self.pains_details) if self.pains_details else '',
            'Aggregator': 1 if self.aggregator else 0,
            'Aggregator_Reason': self.aggregator_reason,
            'Redox': 1 if self.redox else 0,
            'Redox_Groups': ', '.join(self.redox_groups) if self.redox_groups else '',
            'Fluorescence': 1 if self.fluorescence else 0,
            'Fluorescence_Scaffolds': ', '.join(self.fluorescence_scaffolds) if self.fluorescence_scaffolds else '',
            'Thiol': 1 if self.thiol else 0,
            'Thiol_Electrophiles': ', '.join(self.thiol_electrophiles) if self.thiol_electrophiles else '',
            'Total_Flags': self.total_flags,
        }


# =============================================================================
# CONSOLIDATED SMARTS PATTERNS
# =============================================================================

# Redox-active functional groups that can cause assay interference
REDOX_PATTERNS = {
    'catechol': 'c1c(O)c(O)ccc1',  # ortho-diphenol
    'quinone': 'C1(=O)C=CC(=O)C=C1',  # benzoquinone
    'disulfide': '[S;D2]-[S;D2]',  # S-S bond
    'thiol': '[SH1]',  # free thiol
    'hydroquinone': 'c1c(O)ccc(O)c1',  # para-diphenol
    'anthraquinone': 'c1ccc2c(c1)C(=O)c1ccccc1C2=O',
    'naphthoquinone': 'C1=CC2=C(C=C1)C(=O)C=CC2=O',
}

# Fluorescent scaffolds that can interfere with fluorescence-based assays
FLUORESCENT_PATTERNS = {
    'flavonoid': 'O=c1cc(-c2ccccc2)oc2ccccc12',  # flavone core
    'coumarin': 'O=C1Oc2ccccc2C=C1',  # 2H-chromen-2-one (coumarin)
    'xanthene': 'c1ccc2c(c1)Cc1ccccc1O2',
    'naphthalene': 'c1ccc2ccccc2c1',
    'anthracene': 'c1ccc2cc3ccccc3cc2c1',
    'stilbene': 'c1ccccc1C=Cc1ccccc1',  # extended conjugation
}

# Thiol-reactive electrophiles that can modify cysteine residues
THIOL_REACTIVE_PATTERNS = {
    'michael_acceptor': '[C;$(C=C)]-[C;$(C=O)]',  # alpha,beta-unsaturated carbonyl
    'acrylamide': 'C=CC(=O)N',
    'maleimide': 'O=C1C=CC(=O)N1',
    'aldehyde': '[CH1](=O)',  # not part of carboxylic acid
    'activated_ester': '[C;$(C(=O)O)][F,Cl,Br,I]',
    'epoxide': 'C1OC1',
    'isothiocyanate': 'N=C=S',
    'vinyl_sulfone': 'C=CS(=O)(=O)',
}


# ============================================================================
# PAINS Detection
# ============================================================================

def check_pains_violations(mol: Chem.Mol) -> Tuple[bool, List[str]]:
    """
    Check for PAINS (Pan-Assay Interference Substructures) using RDKit FilterCatalog.

    Args:
        mol: RDKit Mol object

    Returns:
        Tuple[bool, List[str]]: (has_pains, list_of_pains_names)
    """
    if mol is None:
        return False, []

    try:
        # Initialize PAINS filter catalog (RDKit built-in)
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        fc = FilterCatalog.FilterCatalog(params)

        # Check for matches
        pains_names = []
        entry = fc.GetFirstMatch(mol)

        if entry is not None:
            pains_names.append(entry.GetDescription())

            # Check for additional matches
            matches = fc.GetMatches(mol)
            for match in matches:
                desc = match.GetDescription()
                if desc not in pains_names:
                    pains_names.append(desc)

        has_pains = len(pains_names) > 0

        if has_pains:
            logger.debug(f"PAINS violations detected: {', '.join(pains_names)}")

        return has_pains, pains_names

    except Exception as e:
        logger.warning(f"Error in PAINS detection: {e}")
        return False, []


# ============================================================================
# Aggregation Risk Detection
# ============================================================================

def check_aggregator_risk(mol: Chem.Mol) -> Tuple[bool, str]:
    """
    Detect aggregation risk using Shoichet laboratory heuristics.

    Risk factors (Shoichet lab criteria):
    - Multiple aromatic rings (>=3)
    - Moderate molecular weight (>300 Da)
    - Low rotatable bonds (<=2, rigid structure)
    - High lipophilicity (LogP > 3)

    Args:
        mol: RDKit Mol object

    Returns:
        Tuple[bool, str]: (is_aggregator_risk, reason)
    """
    if mol is None:
        return False, ""

    try:
        # Calculate molecular descriptors
        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        mw = rdMolDescriptors.CalcExactMolWt(mol)
        logp = Descriptors.MolLogP(mol)

        # Aggregator risk heuristic
        risk_factors = []

        if num_aromatic_rings >= 3:
            risk_factors.append(f"{num_aromatic_rings} aromatic rings")

        if mw > 300:
            risk_factors.append(f"MW={mw:.1f}")

        if num_rotatable_bonds <= 2:
            risk_factors.append(f"{num_rotatable_bonds} rotatable bonds")

        if logp > 3:
            risk_factors.append(f"LogP={logp:.2f}")

        # Risk if meets ALL four criteria (conservative)
        is_risk = len(risk_factors) >= 4
        reason = "; ".join(risk_factors) if is_risk else ""

        if is_risk:
            logger.debug(f"Aggregator risk detected: {reason}")

        return is_risk, reason

    except Exception as e:
        logger.warning(f"Error in aggregator detection: {e}")
        return False, ""


# ============================================================================
# Redox Reactivity Detection
# ============================================================================

def check_redox_reactive(mol: Chem.Mol) -> Tuple[bool, List[str]]:
    """
    Detect redox-active functional groups that can cause assay interference.

    Args:
        mol: RDKit Mol object

    Returns:
        Tuple[bool, List[str]]: (is_redox_reactive, list_of_groups)
    """
    if mol is None:
        return False, []

    try:
        detected_groups = []

        for group_name, smarts in REDOX_PATTERNS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                detected_groups.append(group_name)

        is_redox = len(detected_groups) > 0

        if is_redox:
            logger.debug(f"Redox-reactive groups detected: {', '.join(detected_groups)}")

        return is_redox, detected_groups

    except Exception as e:
        logger.warning(f"Error in redox detection: {e}")
        return False, []


# ============================================================================
# Fluorescence Interference Detection
# ============================================================================

def check_fluorescence_interference(mol: Chem.Mol) -> Tuple[bool, List[str]]:
    """
    Detect compounds likely to cause fluorescence interference in assays.

    Args:
        mol: RDKit Mol object

    Returns:
        Tuple[bool, List[str]]: (is_fluorescent, list_of_scaffold_types)
    """
    if mol is None:
        return False, []

    try:
        detected_scaffolds = []

        for scaffold_name, smarts in FLUORESCENT_PATTERNS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                detected_scaffolds.append(scaffold_name)

        # Check for extended conjugation (>4 conjugated double bonds)
        num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        if num_aromatic_rings >= 3:
            if 'extended_conjugation' not in detected_scaffolds:
                detected_scaffolds.append('extended_conjugation')

        is_fluorescent = len(detected_scaffolds) > 0

        if is_fluorescent:
            logger.debug(f"Fluorescent scaffolds detected: {', '.join(detected_scaffolds)}")

        return is_fluorescent, detected_scaffolds

    except Exception as e:
        logger.warning(f"Error in fluorescence detection: {e}")
        return False, []


# ============================================================================
# Thiol Reactivity Detection
# ============================================================================

def check_thiol_reactive(mol: Chem.Mol) -> Tuple[bool, List[str]]:
    """
    Detect electrophilic groups that react with cysteine thiols in proteins.

    Args:
        mol: RDKit Mol object

    Returns:
        Tuple[bool, List[str]]: (is_thiol_reactive, list_of_electrophiles)
    """
    if mol is None:
        return False, []

    try:
        detected_electrophiles = []

        for group_name, smarts in THIOL_REACTIVE_PATTERNS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                detected_electrophiles.append(group_name)

        is_reactive = len(detected_electrophiles) > 0

        if is_reactive:
            logger.debug(f"Thiol-reactive groups detected: {', '.join(detected_electrophiles)}")

        return is_reactive, detected_electrophiles

    except Exception as e:
        logger.warning(f"Error in thiol reactivity detection: {e}")
        return False, []


# ============================================================================
# Main Interface Functions
# ============================================================================

def calculate_interference_flags(mol: Chem.Mol) -> InterferenceFlags:
    """
    Calculate all assay interference flags for a molecule.

    Args:
        mol: RDKit Mol object

    Returns:
        InterferenceFlags dataclass with all detection results
    """
    flags = InterferenceFlags()

    if mol is None:
        return flags

    try:
        # PAINS
        flags.pains, flags.pains_details = check_pains_violations(mol)

        # Aggregator
        flags.aggregator, flags.aggregator_reason = check_aggregator_risk(mol)

        # Redox
        flags.redox, flags.redox_groups = check_redox_reactive(mol)

        # Fluorescence
        flags.fluorescence, flags.fluorescence_scaffolds = check_fluorescence_interference(mol)

        # Thiol
        flags.thiol, flags.thiol_electrophiles = check_thiol_reactive(mol)

    except Exception as e:
        logger.error(f"Error calculating interference flags: {e}")

    return flags


def get_interference_flags_from_smiles(smiles: str) -> InterferenceFlags:
    """
    Calculate all assay interference flags from a SMILES string.

    Args:
        smiles: SMILES string

    Returns:
        InterferenceFlags dataclass with all detection results
    """
    if not smiles or smiles == 'N/A':
        return InterferenceFlags()

    try:
        mol = Chem.MolFromSmiles(smiles)
        return calculate_interference_flags(mol)
    except Exception as e:
        logger.error(f"Error processing SMILES '{smiles}': {e}")
        return InterferenceFlags()


def get_interference_summary(flags: InterferenceFlags) -> Dict[str, Any]:
    """
    Get a summary of interference flags for display.

    Args:
        flags: InterferenceFlags object

    Returns:
        Dictionary with summary statistics
    """
    return {
        'total_flags': flags.total_flags,
        'is_clean': flags.is_clean,
        'flags': flags.to_dict(),
        'details': flags.to_detailed_dict(),
    }


# Flag descriptions for UI display
FLAG_DESCRIPTIONS = {
    'PAINS': {
        'name': 'PAINS',
        'full_name': 'Pan-Assay Interference',
        'description': 'Pan-Assay Interference Substructures - compounds that show activity in many assays due to interference rather than genuine binding.',
        'color': '#ff6b6b',  # Red
    },
    'Aggregator': {
        'name': 'Aggregator',
        'full_name': 'Colloidal Aggregation',
        'description': 'Compounds that form colloidal aggregates in aqueous solution, leading to non-specific protein inhibition.',
        'color': '#ffa94d',  # Orange
    },
    'Redox': {
        'name': 'Redox',
        'full_name': 'Redox Cycling',
        'description': 'Redox-active compounds that can generate reactive oxygen species or interfere with redox-based assays.',
        'color': '#ffd43b',  # Yellow
    },
    'Fluorescence': {
        'name': 'Fluorescence',
        'full_name': 'Fluorescence Interference',
        'description': 'Autofluorescent compounds that can interfere with fluorescence-based assays.',
        'color': '#74c0fc',  # Blue
    },
    'Thiol': {
        'name': 'Thiol',
        'full_name': 'Thiol Reactivity',
        'description': 'Electrophilic compounds that can covalently modify cysteine residues, causing non-specific inhibition.',
        'color': '#b197fc',  # Purple
    },
}
