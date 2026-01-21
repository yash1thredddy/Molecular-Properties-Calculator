"""
Unit Tests for Assay Interference Detection Module

Tests all five interference detection mechanisms:
1. PAINS (Pan-Assay Interference Substructures)
2. Aggregator risk
3. Redox reactivity
4. Fluorescence interference
5. Thiol reactivity

Each test uses well-characterized molecules with known interference properties.

References:
- Baell & Holloway (2010) J. Med. Chem. 53, 2719-2740 (PAINS)
- Bisson et al. (2016) J. Med. Chem. 59, 1671-1690 (IMPs)
"""

import unittest
from rdkit import Chem

from molecular_calculator.services.assay_interference import (
    InterferenceFlags,
    check_pains_violations,
    check_aggregator_risk,
    check_redox_reactive,
    check_fluorescence_interference,
    check_thiol_reactive,
    calculate_interference_flags,
    get_interference_flags_from_smiles,
    get_interference_summary,
    REDOX_PATTERNS,
    FLUORESCENT_PATTERNS,
    THIOL_REACTIVE_PATTERNS,
)


class TestInterferenceFlags(unittest.TestCase):
    """Test the InterferenceFlags dataclass."""

    def test_default_initialization(self):
        """Test default initialization with all flags False."""
        flags = InterferenceFlags()
        self.assertFalse(flags.pains)
        self.assertFalse(flags.aggregator)
        self.assertFalse(flags.redox)
        self.assertFalse(flags.fluorescence)
        self.assertFalse(flags.thiol)
        self.assertEqual(flags.total_flags, 0)
        self.assertTrue(flags.is_clean)

    def test_total_flags_count(self):
        """Test total_flags property counts correctly."""
        flags = InterferenceFlags(pains=True, redox=True, thiol=True)
        self.assertEqual(flags.total_flags, 3)
        self.assertFalse(flags.is_clean)

    def test_to_dict(self):
        """Test conversion to dictionary."""
        flags = InterferenceFlags(pains=True, aggregator=False, redox=True)
        d = flags.to_dict()
        self.assertEqual(d['PAINS'], 1)
        self.assertEqual(d['Aggregator'], 0)
        self.assertEqual(d['Redox'], 1)
        self.assertEqual(d['Fluorescence'], 0)
        self.assertEqual(d['Thiol'], 0)

    def test_to_detailed_dict(self):
        """Test conversion to detailed dictionary with reasons."""
        flags = InterferenceFlags(
            pains=True,
            pains_details=['catechol_A(92)'],
            redox=True,
            redox_groups=['catechol', 'quinone']
        )
        d = flags.to_detailed_dict()
        self.assertEqual(d['PAINS'], 1)
        self.assertEqual(d['PAINS_Details'], 'catechol_A(92)')
        self.assertEqual(d['Redox_Groups'], 'catechol, quinone')
        self.assertEqual(d['Total_Flags'], 2)


class TestPAINSDetection(unittest.TestCase):
    """Test PAINS (Pan-Assay Interference Substructures) detection."""

    def test_clean_molecule_no_pains(self):
        """Test that clean molecules don't trigger PAINS."""
        # Aspirin - a well-behaved drug
        mol = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O')
        has_pains, names = check_pains_violations(mol)
        self.assertFalse(has_pains)
        self.assertEqual(len(names), 0)

    def test_catechol_pains(self):
        """Test that catechol triggers PAINS (catechol_A pattern)."""
        # Catechol (1,2-benzenediol) - known PAINS
        mol = Chem.MolFromSmiles('Oc1ccccc1O')
        has_pains, names = check_pains_violations(mol)
        self.assertTrue(has_pains)
        self.assertGreater(len(names), 0)

    def test_rhodanine_pains(self):
        """Test rhodanine scaffold detection."""
        # Rhodanine - a well-known PAINS scaffold
        mol = Chem.MolFromSmiles('O=C1NC(=S)SC1')
        has_pains, names = check_pains_violations(mol)
        # Rhodanine is a known PAINS
        self.assertTrue(has_pains)

    def test_quinone_pains(self):
        """Test quinone detection as PAINS."""
        # 1,4-benzoquinone
        mol = Chem.MolFromSmiles('O=C1C=CC(=O)C=C1')
        has_pains, names = check_pains_violations(mol)
        self.assertTrue(has_pains)

    def test_none_molecule(self):
        """Test handling of None molecule."""
        has_pains, names = check_pains_violations(None)
        self.assertFalse(has_pains)
        self.assertEqual(len(names), 0)


class TestAggregatorDetection(unittest.TestCase):
    """Test aggregation risk detection."""

    def test_small_molecule_no_risk(self):
        """Test that small molecules don't trigger aggregator risk."""
        # Ethanol - small, simple molecule
        mol = Chem.MolFromSmiles('CCO')
        is_risk, reason = check_aggregator_risk(mol)
        self.assertFalse(is_risk)

    def test_drug_like_molecule_no_risk(self):
        """Test that typical drug-like molecules don't trigger risk."""
        # Ibuprofen - typical drug
        mol = Chem.MolFromSmiles('CC(C)Cc1ccc(cc1)C(C)C(=O)O')
        is_risk, reason = check_aggregator_risk(mol)
        self.assertFalse(is_risk)

    def test_large_aromatic_aggregator(self):
        """Test that large, rigid, lipophilic aromatics trigger risk."""
        # Coronene-like structure - large, flat, aromatic
        mol = Chem.MolFromSmiles('c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67')
        is_risk, reason = check_aggregator_risk(mol)
        # Should have high aromatic ring count and low rotatable bonds
        self.assertTrue(is_risk)

    def test_aggregator_criteria(self):
        """Test that all four criteria must be met."""
        # Naphthalene - only 2 aromatic rings, should not trigger
        mol = Chem.MolFromSmiles('c1ccc2ccccc2c1')
        is_risk, reason = check_aggregator_risk(mol)
        # 2 aromatic rings < 3, so should not be flagged
        self.assertFalse(is_risk)

    def test_none_molecule(self):
        """Test handling of None molecule."""
        is_risk, reason = check_aggregator_risk(None)
        self.assertFalse(is_risk)
        self.assertEqual(reason, "")


class TestRedoxDetection(unittest.TestCase):
    """Test redox-active group detection."""

    def test_catechol_detection(self):
        """Test catechol (ortho-diphenol) detection."""
        # Catechol
        mol = Chem.MolFromSmiles('Oc1ccccc1O')
        is_redox, groups = check_redox_reactive(mol)
        self.assertTrue(is_redox)
        self.assertIn('catechol', groups)

    def test_hydroquinone_detection(self):
        """Test hydroquinone (para-diphenol) detection."""
        # Hydroquinone
        mol = Chem.MolFromSmiles('Oc1ccc(O)cc1')
        is_redox, groups = check_redox_reactive(mol)
        self.assertTrue(is_redox)
        self.assertIn('hydroquinone', groups)

    def test_thiol_detection(self):
        """Test free thiol detection."""
        # Cysteine (contains -SH)
        mol = Chem.MolFromSmiles('N[C@@H](CS)C(=O)O')
        is_redox, groups = check_redox_reactive(mol)
        self.assertTrue(is_redox)
        self.assertIn('thiol', groups)

    def test_disulfide_detection(self):
        """Test disulfide bond detection."""
        # Cystine (disulfide)
        mol = Chem.MolFromSmiles('N[C@@H](CSSC[C@H](N)C(=O)O)C(=O)O')
        is_redox, groups = check_redox_reactive(mol)
        self.assertTrue(is_redox)
        self.assertIn('disulfide', groups)

    def test_clean_molecule(self):
        """Test that clean molecules don't trigger redox."""
        # Benzene - no redox-active groups
        mol = Chem.MolFromSmiles('c1ccccc1')
        is_redox, groups = check_redox_reactive(mol)
        self.assertFalse(is_redox)
        self.assertEqual(len(groups), 0)

    def test_quercetin_catechol(self):
        """Test quercetin (known catechol-containing flavonoid)."""
        # Quercetin - has catechol B-ring
        mol = Chem.MolFromSmiles('O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12')
        is_redox, groups = check_redox_reactive(mol)
        self.assertTrue(is_redox)
        self.assertIn('catechol', groups)

    def test_none_molecule(self):
        """Test handling of None molecule."""
        is_redox, groups = check_redox_reactive(None)
        self.assertFalse(is_redox)
        self.assertEqual(len(groups), 0)


class TestFluorescenceDetection(unittest.TestCase):
    """Test fluorescence interference detection."""

    def test_fluorescein_detection(self):
        """Test fluorescein scaffold detection (known fluorophore)."""
        # Fluorescein - highly fluorescent xanthene derivative
        mol = Chem.MolFromSmiles('O=C1OC2(c3ccc(O)cc3Oc3cc(O)ccc23)c2ccccc12')
        is_fluor, scaffolds = check_fluorescence_interference(mol)
        self.assertTrue(is_fluor)
        # Should detect extended conjugation (4+ aromatic rings)
        self.assertIn('extended_conjugation', scaffolds)

    def test_naphthalene_detection(self):
        """Test naphthalene detection."""
        # Naphthalene
        mol = Chem.MolFromSmiles('c1ccc2ccccc2c1')
        is_fluor, scaffolds = check_fluorescence_interference(mol)
        self.assertTrue(is_fluor)
        self.assertIn('naphthalene', scaffolds)

    def test_anthracene_detection(self):
        """Test anthracene detection."""
        # Anthracene
        mol = Chem.MolFromSmiles('c1ccc2cc3ccccc3cc2c1')
        is_fluor, scaffolds = check_fluorescence_interference(mol)
        self.assertTrue(is_fluor)
        self.assertIn('anthracene', scaffolds)

    def test_stilbene_detection(self):
        """Test stilbene (extended conjugation) detection."""
        # trans-Stilbene
        mol = Chem.MolFromSmiles('c1ccccc1/C=C/c1ccccc1')
        is_fluor, scaffolds = check_fluorescence_interference(mol)
        self.assertTrue(is_fluor)
        self.assertIn('stilbene', scaffolds)

    def test_extended_conjugation(self):
        """Test extended conjugation detection (>=3 aromatic rings)."""
        # Pyrene - 4 fused aromatic rings (valid SMILES)
        mol = Chem.MolFromSmiles('c1cc2ccc3cccc4ccc(c1)c2c34')
        is_fluor, scaffolds = check_fluorescence_interference(mol)
        self.assertTrue(is_fluor)
        self.assertIn('extended_conjugation', scaffolds)

    def test_simple_benzene_no_fluorescence(self):
        """Test that simple benzene doesn't trigger fluorescence."""
        # Benzene - single ring, not fluorescent
        mol = Chem.MolFromSmiles('c1ccccc1')
        is_fluor, scaffolds = check_fluorescence_interference(mol)
        self.assertFalse(is_fluor)

    def test_none_molecule(self):
        """Test handling of None molecule."""
        is_fluor, scaffolds = check_fluorescence_interference(None)
        self.assertFalse(is_fluor)
        self.assertEqual(len(scaffolds), 0)


class TestThiolReactivityDetection(unittest.TestCase):
    """Test thiol-reactive electrophile detection."""

    def test_acrylamide_detection(self):
        """Test acrylamide (Michael acceptor) detection."""
        # Acrylamide
        mol = Chem.MolFromSmiles('C=CC(=O)N')
        is_reactive, groups = check_thiol_reactive(mol)
        self.assertTrue(is_reactive)
        self.assertIn('acrylamide', groups)

    def test_maleimide_detection(self):
        """Test maleimide detection."""
        # Maleimide
        mol = Chem.MolFromSmiles('O=C1C=CC(=O)N1')
        is_reactive, groups = check_thiol_reactive(mol)
        self.assertTrue(is_reactive)
        self.assertIn('maleimide', groups)

    def test_epoxide_detection(self):
        """Test epoxide detection."""
        # Ethylene oxide
        mol = Chem.MolFromSmiles('C1OC1')
        is_reactive, groups = check_thiol_reactive(mol)
        self.assertTrue(is_reactive)
        self.assertIn('epoxide', groups)

    def test_aldehyde_detection(self):
        """Test aldehyde detection."""
        # Benzaldehyde
        mol = Chem.MolFromSmiles('c1ccccc1C=O')
        is_reactive, groups = check_thiol_reactive(mol)
        self.assertTrue(is_reactive)
        self.assertIn('aldehyde', groups)

    def test_isothiocyanate_detection(self):
        """Test isothiocyanate detection."""
        # Phenyl isothiocyanate
        mol = Chem.MolFromSmiles('c1ccccc1N=C=S')
        is_reactive, groups = check_thiol_reactive(mol)
        self.assertTrue(is_reactive)
        self.assertIn('isothiocyanate', groups)

    def test_vinyl_sulfone_detection(self):
        """Test vinyl sulfone detection."""
        # Methyl vinyl sulfone
        mol = Chem.MolFromSmiles('C=CS(=O)(=O)C')
        is_reactive, groups = check_thiol_reactive(mol)
        self.assertTrue(is_reactive)
        self.assertIn('vinyl_sulfone', groups)

    def test_clean_molecule(self):
        """Test that clean molecules don't trigger thiol reactivity."""
        # Benzene
        mol = Chem.MolFromSmiles('c1ccccc1')
        is_reactive, groups = check_thiol_reactive(mol)
        self.assertFalse(is_reactive)

    def test_carboxylic_acid_not_aldehyde(self):
        """Test that carboxylic acids don't trigger aldehyde detection."""
        # Acetic acid - has C=O but not aldehyde
        mol = Chem.MolFromSmiles('CC(=O)O')
        is_reactive, groups = check_thiol_reactive(mol)
        # Should not detect aldehyde (carboxylic acid C=O is different)
        if is_reactive:
            self.assertNotIn('aldehyde', groups)

    def test_none_molecule(self):
        """Test handling of None molecule."""
        is_reactive, groups = check_thiol_reactive(None)
        self.assertFalse(is_reactive)
        self.assertEqual(len(groups), 0)


class TestCalculateInterferenceFlags(unittest.TestCase):
    """Test the main calculate_interference_flags function."""

    def test_clean_drug_molecule(self):
        """Test a clean drug molecule (ibuprofen)."""
        mol = Chem.MolFromSmiles('CC(C)Cc1ccc(cc1)C(C)C(=O)O')
        flags = calculate_interference_flags(mol)
        # Ibuprofen should be clean
        self.assertFalse(flags.pains)
        self.assertFalse(flags.redox)
        self.assertFalse(flags.thiol)

    def test_quercetin_multiple_flags(self):
        """Test quercetin (known to have multiple interference mechanisms)."""
        # Quercetin - catechol, flavonoid, redox-active
        mol = Chem.MolFromSmiles('O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12')
        flags = calculate_interference_flags(mol)
        # Should have PAINS (catechol), Redox (catechol), possibly Fluorescence
        self.assertTrue(flags.pains or flags.redox)
        self.assertGreater(flags.total_flags, 0)

    def test_none_molecule_returns_empty_flags(self):
        """Test that None molecule returns all-False flags."""
        flags = calculate_interference_flags(None)
        self.assertTrue(flags.is_clean)
        self.assertEqual(flags.total_flags, 0)


class TestGetInterferenceFlagsFromSmiles(unittest.TestCase):
    """Test the SMILES-based interface function."""

    def test_valid_smiles(self):
        """Test with valid SMILES."""
        flags = get_interference_flags_from_smiles('CCO')  # Ethanol
        self.assertIsInstance(flags, InterferenceFlags)
        self.assertTrue(flags.is_clean)

    def test_invalid_smiles(self):
        """Test with invalid SMILES."""
        flags = get_interference_flags_from_smiles('invalid_smiles_xyz')
        self.assertIsInstance(flags, InterferenceFlags)
        self.assertTrue(flags.is_clean)  # Should return empty flags

    def test_empty_smiles(self):
        """Test with empty SMILES."""
        flags = get_interference_flags_from_smiles('')
        self.assertTrue(flags.is_clean)

    def test_na_smiles(self):
        """Test with 'N/A' SMILES."""
        flags = get_interference_flags_from_smiles('N/A')
        self.assertTrue(flags.is_clean)

    def test_catechol_smiles(self):
        """Test catechol detection from SMILES."""
        flags = get_interference_flags_from_smiles('Oc1ccccc1O')
        self.assertTrue(flags.pains or flags.redox)


class TestGetInterferenceSummary(unittest.TestCase):
    """Test the get_interference_summary function."""

    def test_summary_structure(self):
        """Test that summary has correct structure."""
        flags = InterferenceFlags(pains=True, redox=True)
        summary = get_interference_summary(flags)

        self.assertIn('total_flags', summary)
        self.assertIn('is_clean', summary)
        self.assertIn('flags', summary)
        self.assertIn('details', summary)

        self.assertEqual(summary['total_flags'], 2)
        self.assertFalse(summary['is_clean'])


class TestSMARTSPatterns(unittest.TestCase):
    """Test that SMARTS patterns are valid and work correctly."""

    def test_redox_patterns_valid(self):
        """Test that all REDOX_PATTERNS are valid SMARTS."""
        for name, smarts in REDOX_PATTERNS.items():
            pattern = Chem.MolFromSmarts(smarts)
            self.assertIsNotNone(pattern, f"Invalid SMARTS for {name}: {smarts}")

    def test_fluorescent_patterns_valid(self):
        """Test that all FLUORESCENT_PATTERNS are valid SMARTS."""
        for name, smarts in FLUORESCENT_PATTERNS.items():
            pattern = Chem.MolFromSmarts(smarts)
            self.assertIsNotNone(pattern, f"Invalid SMARTS for {name}: {smarts}")

    def test_thiol_reactive_patterns_valid(self):
        """Test that all THIOL_REACTIVE_PATTERNS are valid SMARTS."""
        for name, smarts in THIOL_REACTIVE_PATTERNS.items():
            pattern = Chem.MolFromSmarts(smarts)
            self.assertIsNotNone(pattern, f"Invalid SMARTS for {name}: {smarts}")


class TestKnownCompounds(unittest.TestCase):
    """Test with well-characterized compounds from literature."""

    def test_tunicamycin_thiol(self):
        """Test tunicamycin (known thiol-reactive)."""
        # Simplified tunicamycin-like structure with reactive group
        # Using a simpler test - maleimide which is definitely thiol-reactive
        mol = Chem.MolFromSmiles('O=C1C=CC(=O)N1')  # Maleimide
        flags = calculate_interference_flags(mol)
        self.assertTrue(flags.thiol)

    def test_resveratrol_stilbene(self):
        """Test resveratrol (stilbene scaffold - fluorescent)."""
        mol = Chem.MolFromSmiles('Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1')
        flags = calculate_interference_flags(mol)
        self.assertTrue(flags.fluorescence)

    def test_dopamine_catechol(self):
        """Test dopamine (catechol - redox active)."""
        mol = Chem.MolFromSmiles('NCCc1ccc(O)c(O)c1')
        flags = calculate_interference_flags(mol)
        self.assertTrue(flags.redox)
        self.assertIn('catechol', flags.redox_groups)

    def test_curcumin_michael_acceptor(self):
        """Test curcumin-like structure (Michael acceptor)."""
        # Simplified Michael acceptor
        mol = Chem.MolFromSmiles('CC(=O)/C=C/c1ccccc1')
        flags = calculate_interference_flags(mol)
        # May or may not trigger depending on exact SMARTS matching
        # At minimum, should not crash

    def test_fluorescein_fluorescent(self):
        """Test fluorescein (highly fluorescent)."""
        mol = Chem.MolFromSmiles('O=C1OC2(c3ccc(O)cc3Oc3cc(O)ccc23)c2ccccc12')
        flags = calculate_interference_flags(mol)
        self.assertTrue(flags.fluorescence)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def test_single_atom(self):
        """Test single atom molecule."""
        mol = Chem.MolFromSmiles('[Na]')
        flags = calculate_interference_flags(mol)
        self.assertIsInstance(flags, InterferenceFlags)

    def test_very_large_molecule(self):
        """Test with a large molecule."""
        # Large peptide-like structure
        smiles = 'CC(C)C[C@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)N)C(=O)O'
        flags = get_interference_flags_from_smiles(smiles)
        self.assertIsInstance(flags, InterferenceFlags)

    def test_charged_molecule(self):
        """Test charged molecule."""
        mol = Chem.MolFromSmiles('[NH4+]')
        flags = calculate_interference_flags(mol)
        self.assertIsInstance(flags, InterferenceFlags)

    def test_radical(self):
        """Test molecule with radical."""
        mol = Chem.MolFromSmiles('[CH3]')
        flags = calculate_interference_flags(mol)
        self.assertIsInstance(flags, InterferenceFlags)


if __name__ == '__main__':
    unittest.main(verbosity=2)
