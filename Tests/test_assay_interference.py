"""
Unit Tests for Assay Interference Detection Module

Tests the interference detection mechanisms using established, peer-reviewed sources:
1. PAINS - RDKit FilterCatalog.PAINS (Baell & Holloway 2010)
2. Aggregator - Shoichet lab heuristics (Irwin et al. 2015)
3. Thiol - Dahlin et al. (2015) HTS electrophile SMARTS
4. Redox - Quinone/catechol SMARTS (Baell 2010, Proj et al. 2022)
5. Fluorescence - Fluorophore scaffold SMARTS (Su et al. 2015)

Each test uses well-characterized molecules with known interference properties.

References:
- Baell & Holloway (2010) J. Med. Chem. 53, 2719-2740 (PAINS)
- Irwin et al. (2015) J. Med. Chem. 58, 7076-7087 (Aggregator)
- Dahlin et al. (2015) J. Med. Chem. 58, 2091-2113 (Thiol-reactive)
- Proj et al. (2022) Drug Discov. Today 27, 1733-1742 (Redox)
- Su et al. (2015) J. Chem. Inf. Model. 55, 434-445 (Fluorescence)
"""

import unittest
from rdkit import Chem

from molecular_calculator.services.assay_interference import (
    InterferenceFlags,
    check_pains_violations,
    check_aggregator_risk,
    check_brenk_alerts,
    check_nih_alerts,
    check_thiol_reactive,
    check_redox_reactive,
    check_fluorescence_interference,
    calculate_interference_flags,
    get_interference_flags_from_smiles,
    get_interference_summary,
    get_all_filter_matches,
    REDOX_ACTIVE_SMARTS as REDOX_PATTERNS,
    FLUORESCENT_SMARTS as FLUORESCENT_PATTERNS,
    THIOL_REACTIVE_SMARTS as THIOL_REACTIVE_PATTERNS,
)


class TestInterferenceFlags(unittest.TestCase):
    """Test the InterferenceFlags dataclass."""

    def test_default_initialization(self):
        """Test default initialization with all flags False."""
        flags = InterferenceFlags()
        self.assertFalse(flags.pains)
        self.assertFalse(flags.aggregator)
        self.assertFalse(flags.thiol)
        self.assertFalse(flags.redox)
        self.assertFalse(flags.fluorescence)
        self.assertEqual(flags.total_flags, 0)
        self.assertTrue(flags.is_clean)

    def test_total_flags_count(self):
        """Test total_flags property counts correctly."""
        flags = InterferenceFlags(pains=True, thiol=True, redox=True)
        self.assertEqual(flags.total_flags, 3)
        self.assertFalse(flags.is_clean)

    def test_to_dict(self):
        """Test conversion to dictionary."""
        flags = InterferenceFlags(pains=True, aggregator=False, thiol=True)
        d = flags.to_dict()
        self.assertEqual(d['PAINS'], 1)
        self.assertEqual(d['Aggregator'], 0)
        self.assertEqual(d['Thiol'], 1)
        self.assertEqual(d['Redox'], 0)
        self.assertEqual(d['Fluorescence'], 0)

    def test_to_detailed_dict(self):
        """Test conversion to detailed dictionary with reasons."""
        flags = InterferenceFlags(
            pains=True,
            pains_details=['catechol_A(92)'],
            thiol=True,
            thiol_details=['aldehyde', 'michael_acceptor']
        )
        d = flags.to_detailed_dict()
        self.assertEqual(d['PAINS'], 1)
        self.assertEqual(d['PAINS_Details'], 'catechol_A(92)')
        self.assertEqual(d['Thiol'], 1)
        # Thiol_Details is a comma-separated string
        self.assertIn('aldehyde', d['Thiol_Details'])
        self.assertIn('michael_acceptor', d['Thiol_Details'])
        self.assertEqual(d['Total_Flags'], 2)

    def test_flag_aliases(self):
        """Test that property aliases work correctly."""
        flags = InterferenceFlags(
            thiol=True,
            thiol_details=['epoxide'],
            redox=True,
            redox_details=['quinone'],
            fluorescence=True,
            fluorescence_details=['coumarin']
        )
        # Test aliases
        self.assertEqual(flags.thiol_electrophiles, ['epoxide'])
        self.assertEqual(flags.redox_groups, ['quinone'])
        self.assertEqual(flags.fluorescence_scaffolds, ['coumarin'])


class TestPAINSDetection(unittest.TestCase):
    """Test PAINS (Pan-Assay Interference Substructures) detection.

    Uses RDKit FilterCatalog.PAINS - peer-reviewed PAINS filter set.
    Reference: Baell & Holloway (2010) J. Med. Chem. 53, 2719-2740
    """

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
        # RDKit returns names like 'catechol_A(92)'
        self.assertTrue(any('catechol' in name.lower() for name in names))

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


class TestBRENKDetection(unittest.TestCase):
    """Test BRENK filter detection (legacy function, still available).

    Uses RDKit FilterCatalog.BRENK - 104 unwanted substructures.
    Reference: Brenk et al. (2008) ChemMedChem 3, 435-444
    """

    def test_clean_molecule_no_brenk(self):
        """Test that clean molecules don't trigger BRENK."""
        # Simple ethanol
        mol = Chem.MolFromSmiles('CCO')
        has_alerts, names = check_brenk_alerts(mol)
        self.assertFalse(has_alerts)

    def test_aldehyde_detection(self):
        """Test aldehyde detection (known BRENK alert)."""
        # Benzaldehyde - reactive aldehyde
        mol = Chem.MolFromSmiles('c1ccccc1C=O')
        has_alerts, names = check_brenk_alerts(mol)
        self.assertTrue(has_alerts)
        self.assertGreater(len(names), 0)

    def test_epoxide_detection(self):
        """Test epoxide detection (reactive group in BRENK)."""
        # Ethylene oxide - reactive epoxide
        mol = Chem.MolFromSmiles('C1OC1')
        has_alerts, names = check_brenk_alerts(mol)
        self.assertTrue(has_alerts)

    def test_michael_acceptor_detection(self):
        """Test Michael acceptor detection."""
        # Acrylamide - Michael acceptor
        mol = Chem.MolFromSmiles('C=CC(=O)N')
        has_alerts, names = check_brenk_alerts(mol)
        self.assertTrue(has_alerts)

    def test_thiol_detection(self):
        """Test free thiol detection (BRENK includes thiols)."""
        # Cysteine (contains -SH)
        mol = Chem.MolFromSmiles('N[C@@H](CS)C(=O)O')
        has_alerts, names = check_brenk_alerts(mol)
        self.assertTrue(has_alerts)
        # BRENK filter may use different naming (e.g., 'thiol_2')
        self.assertGreater(len(names), 0)

    def test_disulfide_detection(self):
        """Test disulfide bond detection."""
        # Cystine (disulfide)
        mol = Chem.MolFromSmiles('N[C@@H](CSSC[C@H](N)C(=O)O)C(=O)O')
        has_alerts, names = check_brenk_alerts(mol)
        self.assertTrue(has_alerts)
        # BRENK uses 'disulphide' spelling
        self.assertGreater(len(names), 0)

    def test_maleimide_detection(self):
        """Test maleimide detection."""
        # Maleimide - thiol-reactive
        mol = Chem.MolFromSmiles('O=C1C=CC(=O)N1')
        has_alerts, names = check_brenk_alerts(mol)
        self.assertTrue(has_alerts)

    def test_isothiocyanate_detection(self):
        """Test isothiocyanate detection."""
        # Phenyl isothiocyanate
        mol = Chem.MolFromSmiles('c1ccccc1N=C=S')
        has_alerts, names = check_brenk_alerts(mol)
        self.assertTrue(has_alerts)

    def test_none_molecule(self):
        """Test handling of None molecule."""
        has_alerts, names = check_brenk_alerts(None)
        self.assertFalse(has_alerts)
        self.assertEqual(len(names), 0)


class TestNIHDetection(unittest.TestCase):
    """Test NIH filter detection (legacy function, still available).

    Uses RDKit FilterCatalog.NIH - problematic functional groups.
    Reference: Jadhav et al. (2009) J. Med. Chem. 53, 37-51
    """

    def test_clean_molecule_no_nih(self):
        """Test that clean molecules don't trigger NIH alerts."""
        # Simple benzene (correct SMILES: 6-membered aromatic ring)
        mol = Chem.MolFromSmiles('c1ccccc1')
        has_alerts, names = check_nih_alerts(mol)
        self.assertFalse(has_alerts)

    def test_none_molecule(self):
        """Test handling of None molecule."""
        has_alerts, names = check_nih_alerts(None)
        self.assertFalse(has_alerts)
        self.assertEqual(len(names), 0)


class TestThiolReactiveDetection(unittest.TestCase):
    """Test Thiol-reactive detection using published SMARTS.

    Reference: Dahlin et al. (2015) J. Med. Chem. 58, 2091-2113
    """

    def test_clean_molecule_no_thiol(self):
        """Test that clean molecules don't trigger thiol alerts."""
        mol = Chem.MolFromSmiles('CCO')  # Ethanol
        has_alerts, patterns = check_thiol_reactive(mol)
        self.assertFalse(has_alerts)

    def test_epoxide_detection(self):
        """Test epoxide detection (SN2 electrophile)."""
        mol = Chem.MolFromSmiles('C1OC1')  # Ethylene oxide
        has_alerts, patterns = check_thiol_reactive(mol)
        self.assertTrue(has_alerts)
        self.assertTrue(any('epoxide' in p.lower() for p in patterns))

    def test_aldehyde_detection(self):
        """Test aldehyde detection (Schiff base former)."""
        mol = Chem.MolFromSmiles('c1ccccc1C=O')  # Benzaldehyde
        has_alerts, patterns = check_thiol_reactive(mol)
        self.assertTrue(has_alerts)
        self.assertTrue(any('aldehyde' in p.lower() for p in patterns))

    def test_isothiocyanate_detection(self):
        """Test isothiocyanate detection."""
        mol = Chem.MolFromSmiles('c1ccccc1N=C=S')  # Phenyl isothiocyanate
        has_alerts, patterns = check_thiol_reactive(mol)
        self.assertTrue(has_alerts)

    def test_none_molecule(self):
        """Test handling of None molecule."""
        has_alerts, patterns = check_thiol_reactive(None)
        self.assertFalse(has_alerts)
        self.assertEqual(len(patterns), 0)


class TestRedoxReactiveDetection(unittest.TestCase):
    """Test Redox-active detection using published SMARTS.

    Reference: Proj et al. (2022) Drug Discov. Today 27, 1733-1742
    """

    def test_clean_molecule_no_redox(self):
        """Test that clean molecules don't trigger redox alerts."""
        mol = Chem.MolFromSmiles('CCO')  # Ethanol
        has_alerts, patterns = check_redox_reactive(mol)
        self.assertFalse(has_alerts)

    def test_quinone_detection(self):
        """Test quinone detection (redox-cycling)."""
        mol = Chem.MolFromSmiles('O=C1C=CC(=O)C=C1')  # p-Benzoquinone
        has_alerts, patterns = check_redox_reactive(mol)
        self.assertTrue(has_alerts)

    def test_catechol_detection(self):
        """Test catechol detection (oxidizes to quinone)."""
        mol = Chem.MolFromSmiles('Oc1ccccc1O')  # Catechol
        has_alerts, patterns = check_redox_reactive(mol)
        self.assertTrue(has_alerts)
        self.assertTrue(any('catechol' in p.lower() for p in patterns))

    def test_none_molecule(self):
        """Test handling of None molecule."""
        has_alerts, patterns = check_redox_reactive(None)
        self.assertFalse(has_alerts)
        self.assertEqual(len(patterns), 0)


class TestFluorescenceDetection(unittest.TestCase):
    """Test Autofluorescence detection using published SMARTS.

    Reference: Su et al. (2015) J. Chem. Inf. Model. 55, 434-445
    """

    def test_clean_molecule_no_fluorescence(self):
        """Test that clean molecules don't trigger fluorescence alerts."""
        mol = Chem.MolFromSmiles('CCO')  # Ethanol
        has_alerts, patterns = check_fluorescence_interference(mol)
        self.assertFalse(has_alerts)

    def test_coumarin_detection(self):
        """Test coumarin detection (known fluorophore)."""
        mol = Chem.MolFromSmiles('O=c1ccc2ccccc2o1')  # Coumarin
        has_alerts, patterns = check_fluorescence_interference(mol)
        self.assertTrue(has_alerts)

    def test_none_molecule(self):
        """Test handling of None molecule."""
        has_alerts, patterns = check_fluorescence_interference(None)
        self.assertFalse(has_alerts)
        self.assertEqual(len(patterns), 0)


class TestAggregatorDetection(unittest.TestCase):
    """Test aggregation risk detection.

    Uses Shoichet lab published heuristics.
    Reference: Irwin et al. (2015) J. Med. Chem. 58, 7076-7087
    """

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


class TestCalculateInterferenceFlags(unittest.TestCase):
    """Test the main calculate_interference_flags function."""

    def test_clean_drug_molecule(self):
        """Test a clean drug molecule (ibuprofen)."""
        mol = Chem.MolFromSmiles('CC(C)Cc1ccc(cc1)C(C)C(=O)O')
        flags = calculate_interference_flags(mol)
        # Ibuprofen should be clean
        self.assertFalse(flags.pains)
        self.assertFalse(flags.aggregator)

    def test_quercetin_multiple_flags(self):
        """Test quercetin (known to have multiple interference mechanisms)."""
        # Quercetin - catechol, flavonoid, redox-active
        mol = Chem.MolFromSmiles('O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12')
        flags = calculate_interference_flags(mol)
        # Should have PAINS (catechol)
        self.assertTrue(flags.pains)
        self.assertGreater(flags.total_flags, 0)

    def test_none_molecule_returns_empty_flags(self):
        """Test that None molecule returns all-False flags."""
        flags = calculate_interference_flags(None)
        self.assertTrue(flags.is_clean)
        self.assertEqual(flags.total_flags, 0)

    def test_aldehyde_triggers_thiol(self):
        """Test that aldehyde triggers thiol-reactive flag."""
        mol = Chem.MolFromSmiles('c1ccccc1C=O')  # Benzaldehyde
        flags = calculate_interference_flags(mol)
        self.assertTrue(flags.thiol)

    def test_epoxide_triggers_thiol(self):
        """Test that epoxide triggers thiol-reactive flag."""
        mol = Chem.MolFromSmiles('C1OC1')  # Ethylene oxide
        flags = calculate_interference_flags(mol)
        self.assertTrue(flags.thiol)

    def test_quinone_triggers_redox(self):
        """Test that quinone triggers redox flag."""
        mol = Chem.MolFromSmiles('O=C1C=CC(=O)C=C1')  # p-Benzoquinone
        flags = calculate_interference_flags(mol)
        self.assertTrue(flags.redox)

    def test_catechol_triggers_redox(self):
        """Test that catechol triggers redox flag."""
        mol = Chem.MolFromSmiles('Oc1ccccc1O')  # Catechol
        flags = calculate_interference_flags(mol)
        self.assertTrue(flags.redox)


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
        self.assertTrue(flags.pains)  # Catechol is PAINS


class TestGetInterferenceSummary(unittest.TestCase):
    """Test the get_interference_summary function."""

    def test_summary_structure(self):
        """Test that summary has correct structure."""
        flags = InterferenceFlags(pains=True, thiol=True)
        summary = get_interference_summary(flags)

        self.assertIn('total_flags', summary)
        self.assertIn('is_clean', summary)
        self.assertIn('flags', summary)
        self.assertIn('details', summary)

        self.assertEqual(summary['total_flags'], 2)
        self.assertFalse(summary['is_clean'])


class TestGetAllFilterMatches(unittest.TestCase):
    """Test the get_all_filter_matches function."""

    def test_returns_all_catalogs(self):
        """Test that all catalogs are returned."""
        mol = Chem.MolFromSmiles('CCO')
        results = get_all_filter_matches(mol)

        self.assertIn('PAINS', results)
        self.assertIn('BRENK', results)
        self.assertIn('NIH', results)
        self.assertIn('ZINC', results)

    def test_none_molecule(self):
        """Test handling of None molecule."""
        results = get_all_filter_matches(None)
        self.assertEqual(results, {})

    def test_catechol_matches_pains(self):
        """Test that catechol triggers PAINS in all filter matches."""
        mol = Chem.MolFromSmiles('Oc1ccccc1O')
        results = get_all_filter_matches(mol)
        self.assertGreater(len(results['PAINS']), 0)


class TestSMARTSPatterns(unittest.TestCase):
    """Test that SMARTS patterns are valid."""

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


class TestSMARTSPatternMatching(unittest.TestCase):
    """Test that SMARTS patterns correctly match molecules."""

    def test_catechol_pattern_matches_catechol(self):
        """Test catechol SMARTS matches actual catechol."""
        pattern = Chem.MolFromSmarts(REDOX_PATTERNS['catechol'])
        self.assertIsNotNone(pattern, "Catechol pattern should compile")

        catechol = Chem.MolFromSmiles('Oc1ccccc1O')
        self.assertTrue(
            catechol.HasSubstructMatch(pattern),
            "Catechol pattern should match catechol"
        )

    def test_epoxide_pattern_matches_epoxide(self):
        """Test epoxide SMARTS matches three-membered ring with oxygen."""
        pattern = Chem.MolFromSmarts(THIOL_REACTIVE_PATTERNS['epoxide'])
        self.assertIsNotNone(pattern, "Epoxide pattern should compile")

        ethylene_oxide = Chem.MolFromSmiles('C1OC1')
        self.assertTrue(
            ethylene_oxide.HasSubstructMatch(pattern),
            "Epoxide pattern should match ethylene oxide"
        )

    def test_aldehyde_pattern_matches_aldehyde(self):
        """Test aldehyde SMARTS matches aldehydes."""
        pattern = Chem.MolFromSmarts(THIOL_REACTIVE_PATTERNS['aldehyde'])
        self.assertIsNotNone(pattern, "Aldehyde pattern should compile")

        benzaldehyde = Chem.MolFromSmiles('c1ccccc1C=O')
        self.assertTrue(
            benzaldehyde.HasSubstructMatch(pattern),
            "Aldehyde pattern should match benzaldehyde"
        )


class TestKnownCompounds(unittest.TestCase):
    """Test with well-characterized compounds from literature."""

    def test_maleimide_thiol_reactive(self):
        """Test maleimide triggers thiol-reactive flag (known thiol-reactive)."""
        mol = Chem.MolFromSmiles('O=C1C=CC(=O)N1')  # Maleimide
        flags = calculate_interference_flags(mol)
        # Maleimide should trigger thiol or pains
        self.assertTrue(flags.thiol or flags.pains)

    def test_dopamine_catechol(self):
        """Test dopamine (catechol - PAINS)."""
        mol = Chem.MolFromSmiles('NCCc1ccc(O)c(O)c1')
        flags = calculate_interference_flags(mol)
        # Dopamine has catechol which triggers PAINS
        self.assertTrue(flags.pains)

    def test_acrylamide_michael_acceptor(self):
        """Test acrylamide (Michael acceptor - thiol-reactive)."""
        mol = Chem.MolFromSmiles('C=CC(=O)N')
        flags = calculate_interference_flags(mol)
        self.assertTrue(flags.thiol)


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
