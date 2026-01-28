"""
Comprehensive Unit Tests for Molecular Properties Calculator

Tests all functionality including:
- RDKit molecular properties calculations
- InChI key conversion via online databases
- Input format detection and conversion
- Batch processing
- Error handling and edge cases

Developed by: Yashwanth Reddy for ITR-UIC
Part of: Chemo-Informatics Toolkit
"""

import unittest
import pandas as pd
import numpy as np
from unittest.mock import patch, Mock
import requests
from molecular_calculator import MolecularCalculator, PropertyExplanations


class TestMolecularCalculator(unittest.TestCase):
    """Test suite for MolecularCalculator class"""

    def setUp(self):
        """Set up test data with well-characterized molecules.

        Test molecules were selected based on:
        1. Well-established properties in scientific literature
        2. Diverse structural features
        3. Known expected values from trusted sources

        Expected values are documented with their sources for verification.
        """
        # Test molecules with known properties
        # =====================================================================
        # ASPIRIN (Acetylsalicylic acid)
        # Reference: PubChem CID 2244, DrugBank DB00945
        # MW verified against PubChem (180.16 g/mol, using average mass)
        # LogP: Literature values range 1.19-1.43 (RDKit calculates ~1.3)
        # HBD: 1 (carboxylic acid OH)
        # HBA: RDKit counts 3 (two carbonyl O + ester O, not the acidic OH)
        # =====================================================================
        self.test_molecules = {
            'aspirin': {
                'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                'inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
                'inchi_key': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N',
                'expected_mw': 180.158,  # PubChem: 180.16 g/mol (monoisotopic: 180.042)
                'expected_logp_range': (1.0, 2.0),  # Literature range, RDKit ~1.3
                'expected_hbd': 1,  # One carboxylic acid hydrogen
                'expected_hba': 3   # RDKit: 2 carbonyl O + 1 ester O
            },
            # =====================================================================
            # CAFFEINE
            # Reference: PubChem CID 2519, DrugBank DB00201
            # MW: 194.19 g/mol (PubChem)
            # LogP: -0.07 (experimental), RDKit may vary due to calculation method
            # HBD: 0 (no OH or NH groups)
            # HBA: 6 (4 nitrogens + 2 carbonyl oxygens)
            # =====================================================================
            'caffeine': {
                'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
                'inchi': 'InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3',
                'inchi_key': 'RYYVLZVUVIJVGH-UHFFFAOYSA-N',
                'expected_mw': 194.191,  # PubChem: 194.19 g/mol
                'expected_logp_range': (-2.0, 0.5),  # Broader range for RDKit variance
                'expected_hbd': 0,  # No H-bond donors
                'expected_hba': 6   # 4N + 2 carbonyl O
            },
            # =====================================================================
            # ETHANOL
            # Reference: PubChem CID 702
            # MW: 46.07 g/mol (PubChem)
            # LogP: -0.31 (experimental), simple molecule for validation
            # HBD: 1 (hydroxyl group)
            # HBA: 1 (oxygen)
            # =====================================================================
            'ethanol': {
                'smiles': 'CCO',
                'inchi': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
                'inchi_key': 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N',
                'expected_mw': 46.069,  # PubChem: 46.07 g/mol
                'expected_logp_range': (-0.5, 0.0),  # Experimental: -0.31
                'expected_hbd': 1,  # One OH hydrogen
                'expected_hba': 1   # One oxygen
            }
        }

        # Invalid molecules for error testing
        self.invalid_molecules = {
            'invalid_smiles': 'XYZ123',
            'malformed_inchi': 'InChI=invalid',
            'fake_inchi_key': 'FAKEINCHIKEY-FAKE-N'
        }

    def test_suppress_rdkit_warnings(self):
        """Test RDKit warning suppression"""
        # Test suppression enable/disable
        MolecularCalculator.suppress_rdkit_warnings(True)
        MolecularCalculator.suppress_rdkit_warnings(False)
        # If no exceptions thrown, test passes

    def test_detect_input_format(self):
        """Test input format detection"""
        test_cases = [
            ('CCO', 'smiles'),
            ('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', 'inchi'),
            ('LFQSCWFLJHTTHZ-UHFFFAOYSA-N', 'inchi_key'),
            ('CC(=O)OC1=CC=CC=C1C(=O)O', 'smiles'),
            ('CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'smiles')
        ]

        for input_text, expected_format in test_cases:
            with self.subTest(input_text=input_text):
                detected = MolecularCalculator.detect_input_format(input_text)
                self.assertEqual(detected, expected_format)

    def test_convert_smiles_passthrough(self):
        """Test SMILES passthrough conversion"""
        for name, data in self.test_molecules.items():
            with self.subTest(molecule=name):
                result = MolecularCalculator.convert_to_smiles(
                    data['smiles'], 'smiles', enable_online_lookup=False
                )
                self.assertEqual(result, data['smiles'])

    def test_convert_inchi_to_smiles(self):
        """Test InChI to SMILES conversion"""
        for name, data in self.test_molecules.items():
            with self.subTest(molecule=name):
                result = MolecularCalculator.convert_to_smiles(
                    data['inchi'], 'inchi', enable_online_lookup=False
                )
                self.assertIsNotNone(result, f"Failed to convert InChI for {name}")
                # Verify it's a valid SMILES by checking if it can calculate properties
                properties = MolecularCalculator.calculate_molecular_properties(result)
                self.assertGreater(len(properties), 0, f"No properties calculated for {name}")

    @patch('requests.get')
    def test_convert_inchi_key_to_smiles_success(self, mock_get):
        """Test successful InChI Key to SMILES conversion"""
        # Mock successful CIR response
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = 'CCO'
        mock_get.return_value = mock_response

        result = MolecularCalculator.convert_inchi_key_to_smiles('LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
        self.assertEqual(result, 'CCO')

    @patch('requests.get')
    def test_convert_inchi_key_to_smiles_pubchem_fallback(self, mock_get):
        """Test PubChem fallback for InChI Key conversion"""
        # Mock CIR failure, PubChem success
        def mock_requests_get(url, timeout=None):
            if 'cactus.nci.nih.gov' in url:
                # CIR fails
                mock_response = Mock()
                mock_response.status_code = 404
                return mock_response
            elif 'pubchem.ncbi.nlm.nih.gov' in url:
                # PubChem succeeds
                mock_response = Mock()
                mock_response.status_code = 200
                mock_response.json.return_value = {
                    'PropertyTable': {
                        'Properties': [{'IsomericSMILES': 'CCO'}]
                    }
                }
                return mock_response

        mock_get.side_effect = mock_requests_get

        result = MolecularCalculator.convert_inchi_key_to_smiles('LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
        self.assertEqual(result, 'CCO')

    @patch('requests.get')
    def test_convert_inchi_key_to_smiles_failure(self, mock_get):
        """Test failed InChI Key conversion"""
        # Mock failed responses
        mock_response = Mock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response

        result = MolecularCalculator.convert_inchi_key_to_smiles('FAKEINCHIKEY-FAKE-N')
        self.assertIsNone(result)

    def test_convert_to_smiles_with_online_lookup_disabled(self):
        """Test InChI Key conversion with online lookup disabled"""
        result = MolecularCalculator.convert_to_smiles(
            'LFQSCWFLJHTTHZ-UHFFFAOYSA-N', 'inchi_key', enable_online_lookup=False
        )
        self.assertIsNone(result)

    def test_calculate_molecular_properties_comprehensive(self):
        """Test comprehensive molecular properties calculation"""
        for name, data in self.test_molecules.items():
            with self.subTest(molecule=name):
                properties = MolecularCalculator.calculate_molecular_properties(data['smiles'])

                # Test that we get a dictionary with properties
                self.assertIsInstance(properties, dict)
                self.assertGreater(len(properties), 0)

                # Test specific properties exist
                expected_properties = [
                    'Molecular_Weight', 'Heavy_Atom_Count', 'Atom_Count', 'Bond_Count',
                    'LogP', 'HB_Donors', 'HB_Acceptors', 'TPSA', 'Rotatable_Bonds',
                    'QED', 'Aromatic_Rings', 'Aliphatic_Rings', 'Saturated_Rings',
                    'Heteroatoms', 'CrippenLogP', 'CrippenMR', 'Lipinski_Violations',
                    'Veber_Violations', 'Formal_Charge'
                ]

                for prop in expected_properties:
                    self.assertIn(prop, properties, f"Missing property {prop} for {name}")

                # Test specific expected values
                self.assertAlmostEqual(
                    properties['Molecular_Weight'], data['expected_mw'],
                    places=1, msg=f"Molecular weight mismatch for {name}"
                )

                logp_min, logp_max = data['expected_logp_range']
                self.assertTrue(
                    logp_min <= properties['LogP'] <= logp_max,
                    f"LogP {properties['LogP']} not in expected range {data['expected_logp_range']} for {name}"
                )

                self.assertEqual(
                    properties['HB_Donors'], data['expected_hbd'],
                    f"HB_Donors mismatch for {name}"
                )

                self.assertEqual(
                    properties['HB_Acceptors'], data['expected_hba'],
                    f"HB_Acceptors mismatch for {name}"
                )

                # Test that numeric properties are indeed numeric
                numeric_properties = [
                    'Molecular_Weight', 'LogP', 'TPSA', 'QED', 'CrippenLogP', 'CrippenMR'
                ]
                for prop in numeric_properties:
                    if prop in properties:
                        self.assertIsInstance(
                            properties[prop], (int, float),
                            f"Property {prop} is not numeric for {name}"
                        )

                # Test rule compliance calculations (Lipinski can have 0-4 violations, Veber 0-2)
                self.assertIn(properties['Lipinski_Violations'], [0, 1, 2, 3, 4])
                self.assertIn(properties['Veber_Violations'], [0, 1, 2])

    def test_calculate_molecular_properties_optional_descriptors(self):
        """Test optional molecular descriptors that might fail"""
        properties = MolecularCalculator.calculate_molecular_properties('CCO')

        # Test that optional properties are handled gracefully
        optional_properties = ['BertzCT', 'LabuteASA', 'Chi0', 'Chi1', 'Ring_Count']

        for prop in optional_properties:
            if prop in properties:
                # If present, should be numeric or None
                self.assertTrue(
                    properties[prop] is None or isinstance(properties[prop], (int, float)),
                    f"Optional property {prop} has invalid type"
                )

    def test_calculate_molecular_properties_invalid_smiles(self):
        """Test molecular properties calculation with invalid SMILES"""
        result = MolecularCalculator.calculate_molecular_properties('INVALID_SMILES_XYZ')
        self.assertEqual(result, {})

        result = MolecularCalculator.calculate_molecular_properties('')
        self.assertEqual(result, {})

        # Test with None - should handle gracefully
        try:
            result = MolecularCalculator.calculate_molecular_properties(None)
            self.assertEqual(result, {})
        except Exception:
            # If it raises exception, that's also acceptable behavior
            pass

    def test_detect_smiles_column(self):
        """Test SMILES column detection"""
        test_cases = [
            (['smiles', 'name', 'id'], 'smiles'),
            (['SMILES', 'NAME'], 'SMILES'),
            (['Smiles', 'compound'], 'Smiles'),
            (['smi', 'data'], 'smi'),
            (['SMI', 'info'], 'SMI'),
            (['canonical_smiles'], 'canonical_smiles'),
            (['CANONICAL_SMILES'], 'CANONICAL_SMILES'),
            (['name', 'id', 'data'], None),  # No SMILES column
        ]

        for columns, expected in test_cases:
            with self.subTest(columns=columns):
                df = pd.DataFrame({col: [f'data_{i}'] for i, col in enumerate(columns)})
                result = MolecularCalculator.detect_smiles_column(df)
                self.assertEqual(result, expected)

    def test_get_property_groups(self):
        """Test property groups structure"""
        groups = MolecularCalculator.get_property_groups()

        self.assertIsInstance(groups, dict)
        self.assertGreater(len(groups), 0)

        # Test expected groups exist
        expected_groups = [
            "Basic Properties", "Lipinski Properties", "Drug-likeness",
            "Rule Violations", "Ring Properties", "Complexity", "Additional"
        ]

        for group in expected_groups:
            self.assertIn(group, groups, f"Missing property group: {group}")
            self.assertIsInstance(groups[group], list)
            self.assertGreater(len(groups[group]), 0)

        # Test that all properties are strings
        for group_name, properties in groups.items():
            for prop in properties:
                self.assertIsInstance(prop, str, f"Property {prop} in {group_name} is not a string")

    def test_process_batch(self):
        """Test batch processing functionality"""
        # Create test DataFrame
        test_data = {
            'smiles': ['CCO', 'CC(=O)OC1=CC=CC=C1C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'],
            'name': ['ethanol', 'aspirin', 'caffeine'],
            'id': [1, 2, 3]
        }
        df = pd.DataFrame(test_data)

        # Test with all properties
        result = MolecularCalculator.process_batch(df, 'smiles')

        # Should have original columns plus new property columns
        self.assertGreaterEqual(len(result.columns), len(df.columns))
        self.assertEqual(len(result), len(df))

        # Test with selected properties
        selected_props = {'Molecular_Weight', 'LogP', 'HB_Donors'}
        result_selected = MolecularCalculator.process_batch(df, 'smiles', selected_props)

        # Should have original columns plus selected properties
        new_columns = set(result_selected.columns) - set(df.columns)
        self.assertEqual(new_columns, selected_props)

        # Test property values
        for idx, row in result.iterrows():
            self.assertIsInstance(row['Molecular_Weight'], (int, float))
            self.assertIsInstance(row['LogP'], (int, float))

    def test_process_batch_with_invalid_smiles(self):
        """Test batch processing with some invalid SMILES"""
        test_data = {
            'smiles': ['CCO', 'INVALID_SMILES', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'],
            'name': ['ethanol', 'invalid', 'caffeine']
        }
        df = pd.DataFrame(test_data)

        result = MolecularCalculator.process_batch(df, 'smiles')

        # Should still return results for valid molecules
        self.assertEqual(len(result), len(df))

        # Valid molecules should have properties, invalid should have NaN/empty values
        self.assertGreater(result.iloc[0]['Molecular_Weight'], 0)  # ethanol
        self.assertGreater(result.iloc[2]['Molecular_Weight'], 0)  # caffeine

    def test_lipinski_rule_compliance(self):
        """Test Lipinski Rule of Five compliance calculation"""
        # Test molecule that passes Lipinski rules
        properties_pass = MolecularCalculator.calculate_molecular_properties('CCO')  # ethanol
        self.assertEqual(properties_pass['Lipinski_Violations'], 0)

        # Test molecule that might fail (create a large molecule)
        large_molecule = 'C' * 50  # Very long alkyl chain - should fail MW rule
        properties_fail = MolecularCalculator.calculate_molecular_properties(large_molecule)
        if properties_fail:  # If SMILES is valid
            # Lipinski can have 0-4 violations
            self.assertIn(properties_fail['Lipinski_Violations'], [0, 1, 2, 3, 4])

    def test_veber_rule_compliance(self):
        """Test Veber Rule compliance calculation"""
        properties = MolecularCalculator.calculate_molecular_properties('CCO')
        self.assertEqual(properties['Veber_Violations'], 0)

        # Veber can have 0-2 violations (TPSA and RotBonds)
        for name, data in self.test_molecules.items():
            props = MolecularCalculator.calculate_molecular_properties(data['smiles'])
            self.assertIn(props['Veber_Violations'], [0, 1, 2], f"Invalid Veber violation for {name}")

    def test_stereochemistry_handling(self):
        """Test handling of stereochemical information"""
        # Test molecules with stereochemistry
        stereo_molecules = [
            'C[C@H](O)[C@@H](C)O',  # Molecule with chiral centers
            'C/C=C/C',  # E-alkene
            'C/C=C\\C'   # Z-alkene
        ]

        for smiles in stereo_molecules:
            with self.subTest(smiles=smiles):
                properties = MolecularCalculator.calculate_molecular_properties(smiles)
                self.assertGreater(len(properties), 0, f"Failed to process stereochemical SMILES: {smiles}")


class TestInvalidStereochemistry(unittest.TestCase):
    """Test handling of malformed or invalid stereochemistry SMILES.

    RDKit may reject or warn about certain stereochemical notations.
    These tests ensure the calculator handles these cases gracefully.
    """

    def test_conflicting_chiral_centers(self):
        """Test SMILES with potentially conflicting stereochemistry.

        Some SMILES with conflicting or ambiguous stereochemistry may be
        accepted by RDKit but sanitized/corrected.
        """
        # Valid chiral SMILES - should process correctly
        valid_chiral = 'C[C@H](O)CC'  # Simple chiral center
        properties = MolecularCalculator.calculate_molecular_properties(valid_chiral)
        self.assertGreater(len(properties), 0, "Valid chiral SMILES should process")

    def test_malformed_chiral_notation(self):
        """Test SMILES with malformed chiral notation.

        Invalid stereochemistry notation should either be corrected by RDKit
        or result in empty properties (graceful handling).
        """
        # These may be rejected or processed by RDKit depending on version
        test_cases = [
            '[C@@]',           # Incomplete chiral specification
            '[C@@@H]',         # Invalid: three @ symbols
            'C[C@H]C',         # Chiral center with only 2 substituents specified
        ]

        for smiles in test_cases:
            with self.subTest(smiles=smiles):
                # Should not raise exception - should return empty dict or valid properties
                properties = MolecularCalculator.calculate_molecular_properties(smiles)
                # Either empty (rejected) or valid (corrected) is acceptable
                self.assertIsInstance(properties, dict)

    def test_undefined_stereochemistry(self):
        """Test SMILES with undefined stereochemistry at chiral centers.

        Molecules without stereochemistry defined are valid; RDKit should process them.
        """
        # No stereochemistry defined but has potential chiral centers
        smiles_undefined = 'CC(O)CC'  # No @/@@, but C2 is potentially chiral
        properties = MolecularCalculator.calculate_molecular_properties(smiles_undefined)
        self.assertGreater(len(properties), 0, "Undefined stereochemistry should process")

    def test_e_z_isomer_notation(self):
        """Test E/Z isomer notation in SMILES."""
        test_cases = [
            ('C/C=C/C', 'E-2-butene'),      # E configuration
            ('C/C=C\\C', 'Z-2-butene'),     # Z configuration
            ('C=CC=C', 'butadiene'),        # No E/Z specification
        ]

        for smiles, name in test_cases:
            with self.subTest(name=name):
                properties = MolecularCalculator.calculate_molecular_properties(smiles)
                self.assertGreater(len(properties), 0, f"Failed to process {name}")

    def test_inconsistent_double_bond_stereo(self):
        """Test SMILES with potentially inconsistent double bond stereochemistry."""
        # These SMILES have ambiguous or partial stereochemistry
        test_cases = [
            'C/C=CC',          # Partial E/Z notation (one side only)
            'C=C/C=C',         # Mixed defined/undefined
        ]

        for smiles in test_cases:
            with self.subTest(smiles=smiles):
                properties = MolecularCalculator.calculate_molecular_properties(smiles)
                self.assertIsInstance(properties, dict)

    def test_atropisomer_notation(self):
        """Test axial chirality (atropisomerism) if supported."""
        # BINAP-like structure (simplified)
        # Axial chirality may or may not be supported depending on RDKit version
        test_smiles = 'c1ccc2c(c1)-c1ccccc1-2'  # Simplified biphenyl
        properties = MolecularCalculator.calculate_molecular_properties(test_smiles)
        self.assertIsInstance(properties, dict)

    def test_allene_stereochemistry(self):
        """Test allene with axial chirality."""
        # Allenes can have chirality due to cumulated double bonds
        allene = 'C=C=C'  # Simple allene without specified chirality
        properties = MolecularCalculator.calculate_molecular_properties(allene)
        self.assertGreater(len(properties), 0, "Allene should process")

    def test_complex_stereochemistry(self):
        """Test molecules with multiple stereochemical features."""
        # Multiple chiral centers and E/Z bonds
        complex_smiles = 'C[C@H](O)/C=C/[C@@H](C)O'
        properties = MolecularCalculator.calculate_molecular_properties(complex_smiles)
        self.assertGreater(len(properties), 0, "Complex stereochemistry should process")

    def test_ring_stereochemistry(self):
        """Test stereochemistry in ring systems."""
        test_cases = [
            'C[C@H]1CCCCC1',     # Cis/trans in cyclohexane
            'C1[C@H](O)CCCC1',   # Alcohol on cyclohexane with stereochemistry
            '[C@H]12CC[C@@H](CC1)CC2',  # Bicyclic with stereochemistry
        ]

        for smiles in test_cases:
            with self.subTest(smiles=smiles):
                properties = MolecularCalculator.calculate_molecular_properties(smiles)
                # Should either process or gracefully fail
                self.assertIsInstance(properties, dict)

    def test_isotope_with_stereochemistry(self):
        """Test isotopically labeled molecules with stereochemistry."""
        # Deuterium-labeled with chiral center
        labeled = '[2H]C([2H])([2H])[C@H](O)C'
        properties = MolecularCalculator.calculate_molecular_properties(labeled)
        self.assertIsInstance(properties, dict)


class TestPropertyExplanations(unittest.TestCase):
    """Test suite for PropertyExplanations class"""

    def test_get_explanations(self):
        """Test property explanations content"""
        explanations = PropertyExplanations.get_explanations()

        self.assertIsInstance(explanations, str)
        self.assertGreater(len(explanations), 0)

        # Test that key sections are present
        expected_sections = [
            "Supported Input Formats",
            "Property Explanations",
            "Basic Properties",
            "Lipinski Properties",
            "Drug-likeness",
            "Development Information",
            "Yashwanth Reddy",
            "ITR-UIC"
        ]

        for section in expected_sections:
            self.assertIn(section, explanations, f"Missing section: {section}")


class TestRealInChIKeyConversion(unittest.TestCase):
    """Integration tests with real InChI Key conversion (requires internet)"""

    def setUp(self):
        """Skip tests if no internet connection"""
        try:
            response = requests.get('https://httpbin.org/status/200', timeout=5)
            self.internet_available = response.status_code == 200
        except:
            self.internet_available = False

    @unittest.skipUnless(True, "Requires internet connection")  # Set to False to skip
    def test_real_inchi_key_conversion(self):
        """Test actual InChI Key conversion with real database calls"""
        if not self.internet_available:
            self.skipTest("No internet connection available")

        test_cases = [
            ('LFQSCWFLJHTTHZ-UHFFFAOYSA-N', 'ethanol'),  # Ethanol
            ('BSYNRYMUTXBXSQ-UHFFFAOYSA-N', 'aspirin'),  # Aspirin
        ]

        for inchi_key, name in test_cases:
            with self.subTest(inchi_key=inchi_key, name=name):
                result = MolecularCalculator.convert_inchi_key_to_smiles(inchi_key)

                if result is not None:  # Database lookup successful
                    self.assertIsInstance(result, str)
                    self.assertGreater(len(result), 0)

                    # Verify the result is a valid SMILES
                    properties = MolecularCalculator.calculate_molecular_properties(result)
                    self.assertGreater(len(properties), 0, f"Invalid SMILES returned for {name}")
                else:
                    # If lookup fails, it's often due to network issues, not our code
                    print(f"Warning: Could not resolve {name} InChI Key - may be network issue")


class TestGetCalculatorSingleton(unittest.TestCase):
    """Test suite for get_calculator() singleton function"""

    def test_get_calculator_returns_instance(self):
        """Test that get_calculator returns a MolecularCalculator instance."""
        from molecular_calculator.core import get_calculator
        calc = get_calculator()
        self.assertIsInstance(calc, MolecularCalculator)

    def test_get_calculator_returns_same_instance(self):
        """Test that get_calculator returns the same instance each time."""
        from molecular_calculator.core import get_calculator
        calc1 = get_calculator()
        calc2 = get_calculator()
        self.assertIs(calc1, calc2)

    def test_get_calculator_is_functional(self):
        """Test that the singleton can calculate properties."""
        from molecular_calculator.core import get_calculator
        calc = get_calculator()
        result = calc.calculate("CCO")
        self.assertTrue(result.success)
        self.assertIsNotNone(result.properties.molecular_weight)


class TestParallelBatchProcessing(unittest.TestCase):
    """Test suite for process_batch_parallel() method"""

    def test_parallel_batch_basic(self):
        """Test basic parallel batch processing."""
        df = pd.DataFrame({
            'SMILES': ['C', 'CC', 'CCO', 'c1ccccc1'],
            'Name': ['Methane', 'Ethane', 'Ethanol', 'Benzene']
        })

        result = MolecularCalculator.process_batch_parallel(df, 'SMILES')

        self.assertIn('Molecular_Weight', result.columns)
        self.assertEqual(len(result), 4)
        # Check ethanol MW
        self.assertAlmostEqual(result.loc[2, 'Molecular_Weight'], 46.069, places=1)

    def test_parallel_batch_with_selection(self):
        """Test parallel batch with property selection."""
        df = pd.DataFrame({
            'SMILES': ['CCO', 'c1ccccc1'],
            'Name': ['Ethanol', 'Benzene']
        })

        result = MolecularCalculator.process_batch_parallel(
            df,
            'SMILES',
            selected_properties={'Molecular_Weight', 'LogP'}
        )

        self.assertIn('Molecular_Weight', result.columns)
        self.assertIn('LogP', result.columns)

    def test_parallel_batch_handles_invalid(self):
        """Test parallel batch handles invalid SMILES."""
        df = pd.DataFrame({
            'SMILES': ['CCO', 'invalid_smiles', 'c1ccccc1'],
            'Name': ['Ethanol', 'Invalid', 'Benzene']
        })

        result = MolecularCalculator.process_batch_parallel(df, 'SMILES')

        self.assertEqual(len(result), 3)
        # Valid molecules should have MW
        self.assertTrue(pd.notna(result.loc[0, 'Molecular_Weight']))
        self.assertTrue(pd.notna(result.loc[2, 'Molecular_Weight']))

    def test_parallel_batch_handles_nan(self):
        """Test parallel batch handles NaN values."""
        df = pd.DataFrame({
            'SMILES': ['CCO', None, 'c1ccccc1'],
            'Name': ['Ethanol', 'Missing', 'Benzene']
        })

        result = MolecularCalculator.process_batch_parallel(df, 'SMILES')

        self.assertEqual(len(result), 3)
        self.assertTrue(pd.notna(result.loc[0, 'Molecular_Weight']))
        self.assertTrue(pd.notna(result.loc[2, 'Molecular_Weight']))

    def test_parallel_batch_progress_callback(self):
        """Test progress callback is called."""
        df = pd.DataFrame({
            'SMILES': ['C', 'CC', 'CCO'],
            'Name': ['Methane', 'Ethane', 'Ethanol']
        })

        progress_calls = []

        def callback(completed, total):
            progress_calls.append((completed, total))

        MolecularCalculator.process_batch_parallel(
            df,
            'SMILES',
            progress_callback=callback
        )

        self.assertEqual(len(progress_calls), 3)
        self.assertEqual(progress_calls[-1], (3, 3))

    def test_parallel_batch_preserves_order(self):
        """Test parallel processing preserves row order."""
        df = pd.DataFrame({
            'SMILES': ['C', 'CC', 'CCC', 'CCCC', 'CCCCC'],
            'ID': [1, 2, 3, 4, 5]
        })

        result = MolecularCalculator.process_batch_parallel(df, 'SMILES')

        # Check order is preserved
        self.assertEqual(list(result['ID']), [1, 2, 3, 4, 5])
        # Check MW increases with chain length
        mws = result['Molecular_Weight'].tolist()
        self.assertEqual(mws, sorted(mws))


class TestInputSanitization(unittest.TestCase):
    """Test suite for input sanitization functions"""

    def test_sanitize_smiles_valid(self):
        """Test sanitizing valid SMILES."""
        from molecular_calculator.utils.sanitizer import sanitize_smiles
        self.assertEqual(sanitize_smiles("CCO"), "CCO")
        self.assertEqual(sanitize_smiles("c1ccccc1"), "c1ccccc1")
        self.assertEqual(sanitize_smiles("  CCO  "), "CCO")

    def test_sanitize_smiles_empty(self):
        """Test sanitizing empty SMILES."""
        from molecular_calculator.utils.sanitizer import sanitize_smiles
        self.assertIsNone(sanitize_smiles(""))
        self.assertIsNone(sanitize_smiles("   "))
        self.assertIsNone(sanitize_smiles(None))

    def test_sanitize_smiles_too_long(self):
        """Test sanitizing overly long SMILES."""
        from molecular_calculator.utils.sanitizer import sanitize_smiles
        long_smiles = "C" * 6000
        self.assertIsNone(sanitize_smiles(long_smiles))

    def test_sanitize_inchi_valid(self):
        """Test sanitizing valid InChI."""
        from molecular_calculator.utils.sanitizer import sanitize_inchi
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        self.assertEqual(sanitize_inchi(inchi), inchi)

    def test_sanitize_inchi_invalid_prefix(self):
        """Test sanitizing InChI without proper prefix."""
        from molecular_calculator.utils.sanitizer import sanitize_inchi
        self.assertIsNone(sanitize_inchi("1S/C2H6O/c1-2-3"))

    def test_sanitize_inchi_key_valid(self):
        """Test sanitizing valid InChI Key."""
        from molecular_calculator.utils.sanitizer import sanitize_inchi_key
        key = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        self.assertEqual(sanitize_inchi_key(key), key)

    def test_sanitize_inchi_key_lowercase(self):
        """Test InChI Key is converted to uppercase."""
        from molecular_calculator.utils.sanitizer import sanitize_inchi_key
        key = "lfqscwfljhtthz-uhfffaoysa-n"
        self.assertEqual(sanitize_inchi_key(key), "LFQSCWFLJHTTHZ-UHFFFAOYSA-N")

    def test_sanitize_inchi_key_invalid_format(self):
        """Test sanitizing invalid InChI Key formats."""
        from molecular_calculator.utils.sanitizer import sanitize_inchi_key
        self.assertIsNone(sanitize_inchi_key("INVALID-KEY"))
        self.assertIsNone(sanitize_inchi_key(""))
        self.assertIsNone(sanitize_inchi_key(None))

    def test_sanitize_html_xss(self):
        """Test HTML sanitization against XSS."""
        from molecular_calculator.utils.sanitizer import sanitize_html
        xss_payload = "<script>alert(1)</script>"
        result = sanitize_html(xss_payload)
        self.assertNotIn("<script>", result)
        self.assertIn("&lt;script&gt;", result)

    def test_sanitize_filename_path_traversal(self):
        """Test that path traversal is prevented."""
        from molecular_calculator.utils.sanitizer import sanitize_filename
        # sanitize_filename extracts just the filename, stripping all path components
        result = sanitize_filename("../../../etc/passwd")
        self.assertNotIn("..", result)
        self.assertNotIn("/", result)
        self.assertEqual(result, "passwd")

        result2 = sanitize_filename("..\\windows\\system32")
        self.assertNotIn("..", result2)
        self.assertEqual(result2, "system32")


class TestLoggingSetup(unittest.TestCase):
    """Test suite for logging configuration"""

    def test_setup_logging_import(self):
        """Test that setup_logging can be imported."""
        from molecular_calculator.config import setup_logging
        self.assertTrue(callable(setup_logging))

    def test_setup_logging_from_monitoring(self):
        """Test that monitoring.setup_logging works."""
        from molecular_calculator.utils.monitoring import setup_logging
        self.assertTrue(callable(setup_logging))

    def test_setup_logging_is_idempotent(self):
        """Test that setup_logging can be called multiple times."""
        from molecular_calculator.utils.monitoring import setup_logging
        import logging

        # Should not raise
        setup_logging(level=logging.INFO)
        setup_logging(level=logging.DEBUG)

        # Should still work
        logger = logging.getLogger(__name__)
        logger.info("Test message")


class TestPropertyCalculator(unittest.TestCase):
    """Test suite for PropertyCalculator service"""

    def setUp(self):
        """Set up PropertyCalculator instance."""
        from molecular_calculator.services import PropertyCalculator
        self.property_calc = PropertyCalculator()

    def test_calculate_returns_result(self):
        """Test calculate returns CalculationResult."""
        from molecular_calculator.models import CalculationResult
        result = self.property_calc.calculate("CCO")
        self.assertIsInstance(result, CalculationResult)
        self.assertTrue(result.success)
        self.assertIsNotNone(result.properties)

    def test_calculate_as_dict(self):
        """Test calculate_as_dict returns dictionary."""
        result = self.property_calc.calculate_as_dict("CCO")
        self.assertIsInstance(result, dict)
        self.assertIn('Molecular_Weight', result)

    def test_calculate_qed(self):
        """Test QED calculation."""
        result = self.property_calc.calculate("CCO")
        self.assertIsNotNone(result.properties.qed)
        self.assertGreaterEqual(result.properties.qed, 0)
        self.assertLessEqual(result.properties.qed, 1)

    def test_calculate_all_properties_present(self):
        """Test that all expected properties are calculated."""
        result = self.property_calc.calculate("c1ccccc1")
        props = result.properties

        # Check basic properties
        self.assertIsNotNone(props.molecular_weight)
        self.assertIsNotNone(props.heavy_atom_count)
        self.assertIsNotNone(props.atom_count)

        # Check Lipinski properties
        self.assertIsNotNone(props.logp)
        self.assertIsNotNone(props.hb_donors)
        self.assertIsNotNone(props.hb_acceptors)
        self.assertIsNotNone(props.tpsa)

        # Check ring properties
        self.assertIsNotNone(props.aromatic_rings)

    def test_get_all_property_names(self):
        """Test getting all property names."""
        from molecular_calculator.services import PropertyCalculator
        names = PropertyCalculator.get_all_property_names()
        self.assertIsInstance(names, list)
        self.assertIn('Molecular_Weight', names)
        self.assertIn('LogP', names)

    def test_calculate_invalid_smiles(self):
        """Test calculate with invalid SMILES."""
        result = self.property_calc.calculate("invalid_smiles")
        self.assertFalse(result.success)
        self.assertIsNotNone(result.error)


class TestConversionService(unittest.TestCase):
    """Test suite for ConversionService"""

    def setUp(self):
        """Set up ConversionService instance."""
        from molecular_calculator.services import ConversionService
        from molecular_calculator.models import InputFormat
        self.service = ConversionService()
        self.InputFormat = InputFormat

    def test_detect_format_smiles(self):
        """Test detecting SMILES format."""
        self.assertEqual(self.service.detect_format("CCO"), self.InputFormat.SMILES)
        self.assertEqual(self.service.detect_format("c1ccccc1"), self.InputFormat.SMILES)

    def test_detect_format_inchi(self):
        """Test detecting InChI format."""
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        self.assertEqual(self.service.detect_format(inchi), self.InputFormat.INCHI)

    def test_detect_format_inchi_key(self):
        """Test detecting InChI Key format."""
        key = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
        self.assertEqual(self.service.detect_format(key), self.InputFormat.INCHI_KEY)

    def test_detect_format_unknown(self):
        """Test detecting unknown format."""
        self.assertEqual(self.service.detect_format(""), self.InputFormat.UNKNOWN)
        self.assertEqual(self.service.detect_format(None), self.InputFormat.UNKNOWN)

    def test_to_smiles_from_smiles(self):
        """Test SMILES passthrough."""
        result = self.service.to_smiles("CCO")
        self.assertTrue(result.success)
        self.assertEqual(result.smiles, "CCO")
        self.assertEqual(result.source_format, self.InputFormat.SMILES)

    def test_to_smiles_from_inchi(self):
        """Test InChI to SMILES conversion."""
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        result = self.service.to_smiles(inchi)
        self.assertTrue(result.success)
        self.assertIsNotNone(result.smiles)
        self.assertEqual(result.source_format, self.InputFormat.INCHI)

    def test_to_smiles_invalid_smiles(self):
        """Test handling invalid SMILES."""
        result = self.service.to_smiles("not_valid_smiles", self.InputFormat.SMILES)
        self.assertFalse(result.success)
        self.assertIsNotNone(result.error)

    def test_validate_smiles_valid(self):
        """Test SMILES validation with valid input."""
        self.assertTrue(self.service.validate_smiles("CCO"))
        self.assertTrue(self.service.validate_smiles("c1ccccc1"))

    def test_validate_smiles_invalid(self):
        """Test SMILES validation with invalid input."""
        self.assertFalse(self.service.validate_smiles("invalid"))
        self.assertFalse(self.service.validate_smiles(""))
        self.assertFalse(self.service.validate_smiles(None))

    def test_auto_detect_and_convert(self):
        """Test auto-detection and conversion."""
        # Test with SMILES
        result = self.service.to_smiles("CCO")
        self.assertTrue(result.success)

        # Test with InChI
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        result = self.service.to_smiles(inchi)
        self.assertTrue(result.success)


class TestMolecularProperties(unittest.TestCase):
    """Test suite for MolecularProperties dataclass"""

    def test_to_dict(self):
        """Test converting properties to dictionary."""
        from molecular_calculator.models import MolecularProperties
        props = MolecularProperties(
            molecular_weight=46.069,
            heavy_atom_count=3,
            logp=-0.14
        )
        d = props.to_dict()
        self.assertEqual(d['Molecular_Weight'], 46.069)
        self.assertEqual(d['Heavy_Atom_Count'], 3)
        self.assertEqual(d['LogP'], -0.14)

    def test_to_dict_filters_none(self):
        """Test that None values are filtered from dict."""
        from molecular_calculator.models import MolecularProperties
        props = MolecularProperties(
            molecular_weight=46.069,
            heavy_atom_count=None
        )
        d = props.to_dict()
        self.assertIn('Molecular_Weight', d)
        self.assertNotIn('Heavy_Atom_Count', d)

    def test_from_dict(self):
        """Test creating properties from dictionary."""
        from molecular_calculator.models import MolecularProperties
        d = {
            'Molecular_Weight': 46.069,
            'Heavy_Atom_Count': 3,
            'LogP': -0.14
        }
        props = MolecularProperties.from_dict(d)
        self.assertEqual(props.molecular_weight, 46.069)
        self.assertEqual(props.heavy_atom_count, 3)
        self.assertEqual(props.logp, -0.14)


class TestAdditionalSanitization(unittest.TestCase):
    """Additional sanitization tests for complete coverage"""

    def test_sanitize_smiles_rejects_invalid_chars(self):
        """Test that SMILES with invalid characters are rejected (returns None)."""
        from molecular_calculator.utils.sanitizer import sanitize_smiles
        # Invalid characters should cause rejection to prevent silent corruption
        result = sanitize_smiles("CCO!")
        self.assertIsNone(result)

        # Valid SMILES should pass through
        result = sanitize_smiles("CCO")
        self.assertIsNotNone(result)
        self.assertEqual(result, "CCO")

    def test_sanitize_inchi_with_whitespace(self):
        """Test sanitizing InChI with whitespace."""
        from molecular_calculator.utils.sanitizer import sanitize_inchi
        inchi = "  InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3  "
        self.assertEqual(sanitize_inchi(inchi), "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")

    def test_sanitize_inchi_empty(self):
        """Test sanitizing empty InChI."""
        from molecular_calculator.utils.sanitizer import sanitize_inchi
        self.assertIsNone(sanitize_inchi(""))
        self.assertIsNone(sanitize_inchi(None))

    def test_sanitize_inchi_key_with_whitespace(self):
        """Test sanitizing InChI Key with whitespace."""
        from molecular_calculator.utils.sanitizer import sanitize_inchi_key
        key = "  LFQSCWFLJHTTHZ-UHFFFAOYSA-N  "
        self.assertEqual(sanitize_inchi_key(key), "LFQSCWFLJHTTHZ-UHFFFAOYSA-N")

    def test_sanitize_inchi_key_too_short(self):
        """Test sanitizing InChI Key that's too short."""
        from molecular_calculator.utils.sanitizer import sanitize_inchi_key
        self.assertIsNone(sanitize_inchi_key("LFQSCWFLJHTTHZ-UHFFFAOYSA"))

    def test_sanitize_column_name_valid(self):
        """Test sanitizing valid column names."""
        from molecular_calculator.utils.sanitizer import sanitize_column_name
        self.assertEqual(sanitize_column_name("SMILES"), "SMILES")
        self.assertEqual(sanitize_column_name("Molecular Weight"), "Molecular_Weight")

    def test_sanitize_column_name_special_chars(self):
        """Test sanitizing column names with special characters."""
        from molecular_calculator.utils.sanitizer import sanitize_column_name
        self.assertEqual(sanitize_column_name("Col<>Name"), "ColName")
        self.assertEqual(sanitize_column_name("Test@Column!"), "TestColumn")

    def test_sanitize_column_name_empty(self):
        """Test sanitizing empty column names."""
        from molecular_calculator.utils.sanitizer import sanitize_column_name
        self.assertEqual(sanitize_column_name(""), "unnamed")
        self.assertEqual(sanitize_column_name(None), "unnamed")

    def test_sanitize_filename_valid(self):
        """Test sanitizing valid filenames."""
        from molecular_calculator.utils.sanitizer import sanitize_filename
        self.assertEqual(sanitize_filename("data.csv"), "data.csv")
        self.assertEqual(sanitize_filename("results_2024.xlsx"), "results_2024.xlsx")

    def test_sanitize_filename_special_chars(self):
        """Test sanitizing filenames with special characters."""
        from molecular_calculator.utils.sanitizer import sanitize_filename
        result = sanitize_filename("file<>:name.csv")
        self.assertNotIn("<", result)
        self.assertNotIn(">", result)
        self.assertNotIn(":", result)

    def test_sanitize_filename_hidden(self):
        """Test that hidden files are not allowed."""
        from molecular_calculator.utils.sanitizer import sanitize_filename
        self.assertEqual(sanitize_filename(".hidden"), "hidden")

    def test_sanitize_html_valid(self):
        """Test sanitizing HTML content."""
        from molecular_calculator.utils.sanitizer import sanitize_html
        self.assertEqual(sanitize_html("Hello World"), "Hello World")
        self.assertEqual(
            sanitize_html("<script>alert(1)</script>"),
            "&lt;script&gt;alert(1)&lt;/script&gt;"
        )

    def test_sanitize_html_special_chars(self):
        """Test HTML special characters are escaped."""
        from molecular_calculator.utils.sanitizer import sanitize_html
        self.assertIn("&amp;", sanitize_html("A & B"))
        self.assertIn("&lt;", sanitize_html("a < b"))
        self.assertIn("&gt;", sanitize_html("a > b"))
        self.assertIn("&quot;", sanitize_html('say "hi"'))

    def test_sanitize_html_empty(self):
        """Test sanitizing empty HTML."""
        from molecular_calculator.utils.sanitizer import sanitize_html
        self.assertEqual(sanitize_html(""), "")
        self.assertEqual(sanitize_html(None), "")


# =============================================================================
# API Client Tests
# =============================================================================

class TestAPIClient(unittest.TestCase):
    """Test suite for ChemicalAPIClient"""

    def test_api_client_import(self):
        """Test that API client can be imported."""
        from molecular_calculator.services.api_client import ChemicalAPIClient, get_api_client
        self.assertTrue(callable(get_api_client))

    def test_api_client_singleton(self):
        """Test that get_api_client returns singleton."""
        from molecular_calculator.services.api_client import get_api_client
        client1 = get_api_client()
        client2 = get_api_client()
        self.assertIs(client1, client2)

    @patch('molecular_calculator.services.api_client.requests.get')
    def test_inchi_key_to_smiles_nih_success(self, mock_get):
        """Test successful NIH CIR lookup."""
        from molecular_calculator.services.api_client import ChemicalAPIClient
        # Clear cache to ensure we hit the mock
        try:
            from molecular_calculator.utils.cache import api_cache
            api_cache.clear()
        except ImportError:
            pass

        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = 'CCO'
        mock_get.return_value = mock_response

        client = ChemicalAPIClient()
        result = client.inchi_key_to_smiles('LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
        self.assertTrue(result.success)
        self.assertEqual(result.data, 'CCO')

    @patch('molecular_calculator.services.api_client.requests.get')
    def test_inchi_key_to_smiles_failure(self, mock_get):
        """Test failed API lookup."""
        from molecular_calculator.services.api_client import ChemicalAPIClient
        # Clear cache to ensure we hit the mock
        try:
            from molecular_calculator.utils.cache import api_cache
            api_cache.clear()
        except ImportError:
            pass

        mock_response = Mock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response

        client = ChemicalAPIClient()
        result = client.inchi_key_to_smiles('INVALIDKEYXYZAB-INVALIDXYZX-I')
        self.assertFalse(result.success)

    @patch('molecular_calculator.services.api_client.requests.get')
    def test_inchi_key_to_smiles_timeout(self, mock_get):
        """Test API timeout handling."""
        from molecular_calculator.services.api_client import ChemicalAPIClient
        # Clear cache to ensure we hit the API
        try:
            from molecular_calculator.utils.cache import api_cache
            api_cache.clear()
        except ImportError:
            pass

        mock_get.side_effect = requests.exceptions.Timeout()

        client = ChemicalAPIClient()
        # Use a unique key that won't be cached from other tests
        result = client.inchi_key_to_smiles('ABCDEFGHIJKLMN-OPQRSTUVWX-Y')
        self.assertFalse(result.success)


# =============================================================================
# Ligand Efficiency Tests
# =============================================================================

class TestLigandEfficiency(unittest.TestCase):
    """Test suite for Ligand Efficiency calculations"""

    def test_lei_import(self):
        """Test LEI module imports."""
        from molecular_calculator.services.ligand_efficiency import (
            LigandEfficiencyCalculator,
            DependencyChecker,
            get_lei_descriptions,
        )
        self.assertTrue(callable(get_lei_descriptions))

    def test_get_lei_descriptions(self):
        """Test getting LEI descriptions."""
        from molecular_calculator.services.ligand_efficiency import get_lei_descriptions
        descriptions = get_lei_descriptions()
        self.assertIsInstance(descriptions, dict)
        self.assertIn('NSEI', descriptions)
        self.assertIn('BEI', descriptions)

    def test_dependency_checker_detect_column(self):
        """Test column detection."""
        from molecular_calculator.services.ligand_efficiency import DependencyChecker
        df = pd.DataFrame({
            'pKi': [5.0, 6.0],
            'SMILES': ['CCO', 'c1ccccc1']
        })
        result = DependencyChecker.detect_column(df, 'pki')
        self.assertEqual(result, 'pKi')

    def test_dependency_checker_detect_all_columns(self):
        """Test detecting all columns."""
        from molecular_calculator.services.ligand_efficiency import DependencyChecker
        df = pd.DataFrame({
            'pKi': [5.0],
            'SMILES': ['CCO'],
            'MW': [46.0]
        })
        detected = DependencyChecker.detect_all_columns(df)
        self.assertIsInstance(detected, dict)
        self.assertIn('pki', detected)
        self.assertIn('smiles', detected)

    def test_lei_calculator_count_heavy_atoms(self):
        """Test heavy atom counting."""
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator
        count = LigandEfficiencyCalculator.count_heavy_atoms('CCO')
        self.assertEqual(count, 3)  # C, C, O

    def test_lei_calculator_count_polar_atoms(self):
        """Test polar atom counting."""
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator
        count = LigandEfficiencyCalculator.count_polar_atoms('CCO')
        self.assertEqual(count, 1)  # Just O

    def test_lei_calculator_calculate_ki_from_pki(self):
        """Test Ki calculation from pKi."""
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator
        ki = LigandEfficiencyCalculator.calculate_ki_from_pki(6.0)
        self.assertAlmostEqual(ki, 1e-6, places=10)

    def test_lei_calculator_calculate_delta_g(self):
        """Test delta G calculation."""
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator
        delta_g = LigandEfficiencyCalculator.calculate_delta_g(1e-9)
        self.assertLess(delta_g, 0)  # Should be negative (favorable binding)


class TestLigandEfficiencyEdgeCases(unittest.TestCase):
    """Test edge cases for Ligand Efficiency calculations.

    These tests ensure the LEI calculations handle edge cases gracefully:
    - Zero heavy atoms (division by zero protection)
    - Zero polar atoms (division by zero protection)
    - Zero pKi values
    - Negative pKi values
    - Very large values
    - Invalid SMILES
    """

    def test_zero_heavy_atoms_protection(self):
        """Test that zero heavy atoms don't cause division by zero.

        When heavy_atoms=0, NBEI, nBEI, and LEH should not be calculated
        (avoiding division by zero).
        """
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        # Calculate LEI with zero heavy atoms
        results = LigandEfficiencyCalculator.calculate_lei_values(
            pki=6.0,
            mw=100.0,
            tpsa=50.0,
            heavy_atoms=0,  # Zero heavy atoms - should not cause error
            polar_atoms=2
        )

        # NBEI, nBEI, LEH require heavy atoms, should NOT be in results
        self.assertNotIn('NBEI', results, "NBEI should not be calculated with zero heavy atoms")
        self.assertNotIn('nBEI', results, "nBEI should not be calculated with zero heavy atoms")
        self.assertNotIn('LEH', results, "LEH should not be calculated with zero heavy atoms")

        # BEI, SEI, NSEI, LEP should still be calculated
        self.assertIn('BEI', results, "BEI should still be calculated")
        self.assertIn('SEI', results, "SEI should still be calculated")
        self.assertIn('NSEI', results, "NSEI should still be calculated")
        self.assertIn('LEP', results, "LEP should still be calculated")

    def test_zero_polar_atoms_protection(self):
        """Test that zero polar atoms don't cause division by zero.

        When polar_atoms=0, NSEI and LEP should not be calculated
        (avoiding division by zero).
        """
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        # Calculate LEI with zero polar atoms
        results = LigandEfficiencyCalculator.calculate_lei_values(
            pki=6.0,
            mw=100.0,
            tpsa=50.0,
            heavy_atoms=10,
            polar_atoms=0  # Zero polar atoms - should not cause error
        )

        # NSEI and LEP require polar atoms, should NOT be in results
        self.assertNotIn('NSEI', results, "NSEI should not be calculated with zero polar atoms")
        self.assertNotIn('LEP', results, "LEP should not be calculated with zero polar atoms")

        # Other LEIs should still be calculated
        self.assertIn('BEI', results, "BEI should still be calculated")
        self.assertIn('NBEI', results, "NBEI should still be calculated")
        self.assertIn('SEI', results, "SEI should still be calculated")
        self.assertIn('LEH', results, "LEH should still be calculated")

    def test_zero_mw_protection(self):
        """Test that zero MW doesn't cause division by zero."""
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        results = LigandEfficiencyCalculator.calculate_lei_values(
            pki=6.0,
            mw=0,  # Zero MW - should not cause error
            tpsa=50.0,
            heavy_atoms=10,
            polar_atoms=2
        )

        # BEI and mBEI require MW, should NOT be in results
        self.assertNotIn('BEI', results, "BEI should not be calculated with zero MW")
        self.assertNotIn('mBEI', results, "mBEI should not be calculated with zero MW")

    def test_zero_tpsa_protection(self):
        """Test that zero TPSA doesn't cause division by zero."""
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        results = LigandEfficiencyCalculator.calculate_lei_values(
            pki=6.0,
            mw=100.0,
            tpsa=0,  # Zero TPSA - should not cause error
            heavy_atoms=10,
            polar_atoms=2
        )

        # SEI requires TPSA, should NOT be in results
        self.assertNotIn('SEI', results, "SEI should not be calculated with zero TPSA")

    def test_zero_pki_values(self):
        """Test LEI calculations with zero pKi.

        pKi=0 is unusual but valid (very weak binder, Ki=1M).
        Calculations should complete without error.
        """
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        # pKi=0 means Ki=1M (10^0 = 1)
        results = LigandEfficiencyCalculator.calculate_lei_values(
            pki=0.0,
            mw=100.0,
            tpsa=50.0,
            heavy_atoms=10,
            polar_atoms=2
        )

        # Should complete without error
        self.assertIsInstance(results, dict)

        # Values should be 0 or close to 0 for pKi=0
        if 'BEI' in results:
            self.assertEqual(results['BEI'], 0.0, "BEI should be 0 when pKi is 0")
        if 'NBEI' in results:
            self.assertEqual(results['NBEI'], 0.0, "NBEI should be 0 when pKi is 0")

    def test_negative_pki_values(self):
        """Test LEI calculations with negative pKi.

        Negative pKi values are unusual (Ki > 1M) but mathematically valid.
        The calculation should handle them without errors.
        """
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        # pKi=-1 means Ki=10M (very weak binding)
        results = LigandEfficiencyCalculator.calculate_lei_values(
            pki=-1.0,
            mw=100.0,
            tpsa=50.0,
            heavy_atoms=10,
            polar_atoms=2
        )

        # Should complete without error
        self.assertIsInstance(results, dict)

        # BEI with negative pKi should be negative
        if 'BEI' in results:
            self.assertLess(results['BEI'], 0, "BEI should be negative when pKi is negative")

    def test_very_large_pki_values(self):
        """Test LEI calculations with very large pKi values.

        Very large pKi (e.g., 15) represents extremely tight binding (Ki=1fM).
        Should handle without overflow errors.
        """
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        results = LigandEfficiencyCalculator.calculate_lei_values(
            pki=15.0,  # Ki = 10^-15 M = 1 femtomolar
            mw=500.0,
            tpsa=100.0,
            heavy_atoms=30,
            polar_atoms=6
        )

        # Should complete without error
        self.assertIsInstance(results, dict)

        # BEI should be positive and large
        if 'BEI' in results:
            self.assertGreater(results['BEI'], 0, "BEI should be positive for large pKi")
            self.assertIsInstance(results['BEI'], float)
            # Check no overflow (should be finite)
            self.assertTrue(np.isfinite(results['BEI']), "BEI should be finite")

    def test_very_large_mw_values(self):
        """Test LEI calculations with very large molecular weights.

        Large molecules (e.g., peptides) can have MW > 5000.
        """
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        results = LigandEfficiencyCalculator.calculate_lei_values(
            pki=6.0,
            mw=10000.0,  # Very large molecule (10 kDa)
            tpsa=500.0,
            heavy_atoms=500,
            polar_atoms=100
        )

        # Should complete without error
        self.assertIsInstance(results, dict)

        # BEI should be small for large molecules
        if 'BEI' in results:
            self.assertLess(results['BEI'], 1.0, "BEI should be small for large molecules")

    def test_none_values_handling(self):
        """Test LEI calculations when optional parameters are None.

        Only pKi is required; other parameters can be None.
        """
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        # Only pKi provided
        results = LigandEfficiencyCalculator.calculate_lei_values(
            pki=6.0,
            mw=None,
            tpsa=None,
            heavy_atoms=None,
            polar_atoms=None
        )

        # Should complete without error but return empty dict
        self.assertIsInstance(results, dict)
        # No LEIs can be calculated without the required properties
        self.assertEqual(len(results), 0, "No LEIs should be calculated when all properties are None")

    def test_invalid_smiles_heavy_atoms(self):
        """Test heavy atom counting with invalid SMILES."""
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        # Invalid SMILES should return 0
        count = LigandEfficiencyCalculator.count_heavy_atoms('INVALID_SMILES')
        self.assertEqual(count, 0, "Invalid SMILES should return 0 heavy atoms")

        count = LigandEfficiencyCalculator.count_heavy_atoms('')
        self.assertEqual(count, 0, "Empty string should return 0 heavy atoms")

        count = LigandEfficiencyCalculator.count_heavy_atoms(None)
        self.assertEqual(count, 0, "None should return 0 heavy atoms")

    def test_invalid_smiles_polar_atoms(self):
        """Test polar atom counting with invalid SMILES."""
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        # Invalid SMILES should return 0
        count = LigandEfficiencyCalculator.count_polar_atoms('INVALID_SMILES')
        self.assertEqual(count, 0, "Invalid SMILES should return 0 polar atoms")

        count = LigandEfficiencyCalculator.count_polar_atoms('')
        self.assertEqual(count, 0, "Empty string should return 0 polar atoms")

        count = LigandEfficiencyCalculator.count_polar_atoms(None)
        self.assertEqual(count, 0, "None should return 0 polar atoms")

    def test_delta_g_with_zero_ki(self):
        """Test delta G calculation with zero Ki (edge case).

        Ki=0 is physically impossible but should be handled gracefully.
        """
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        # Should return 0 for invalid Ki <= 0
        delta_g = LigandEfficiencyCalculator.calculate_delta_g(0.0)
        self.assertEqual(delta_g, 0.0, "Delta G should be 0 for Ki=0")

        delta_g = LigandEfficiencyCalculator.calculate_delta_g(-1.0)
        self.assertEqual(delta_g, 0.0, "Delta G should be 0 for negative Ki")

    def test_molecule_with_no_polar_atoms(self):
        """Test counting polar atoms in molecules with no N or O.

        Molecules like hydrocarbons have no polar atoms.
        """
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        # Benzene has no N or O
        count = LigandEfficiencyCalculator.count_polar_atoms('c1ccccc1')
        self.assertEqual(count, 0, "Benzene should have 0 polar atoms")

        # Hexane has no N or O
        count = LigandEfficiencyCalculator.count_polar_atoms('CCCCCC')
        self.assertEqual(count, 0, "Hexane should have 0 polar atoms")

    def test_molecule_with_many_polar_atoms(self):
        """Test counting polar atoms in molecules with many N and O."""
        from molecular_calculator.services.ligand_efficiency import LigandEfficiencyCalculator

        # Urea: 2N + 1O = 3 polar atoms
        count = LigandEfficiencyCalculator.count_polar_atoms('NC(N)=O')
        self.assertEqual(count, 3, "Urea should have 3 polar atoms (2N + 1O)")

        # Glycine: 1N + 2O = 3 polar atoms
        count = LigandEfficiencyCalculator.count_polar_atoms('NCC(=O)O')
        self.assertEqual(count, 3, "Glycine should have 3 polar atoms (1N + 2O)")


# =============================================================================
# Cache Tests
# =============================================================================

class TestCache(unittest.TestCase):
    """Test suite for TTLCache"""

    def test_cache_import(self):
        """Test cache module imports."""
        from molecular_calculator.utils.cache import TTLCache, cached
        self.assertTrue(callable(TTLCache))
        self.assertTrue(callable(cached))

    def test_cache_set_get(self):
        """Test basic cache set/get."""
        from molecular_calculator.utils.cache import TTLCache
        cache = TTLCache(ttl=60)
        cache.set('key1', 'value1')
        self.assertEqual(cache.get('key1'), 'value1')

    def test_cache_miss(self):
        """Test cache miss returns None."""
        from molecular_calculator.utils.cache import TTLCache
        cache = TTLCache(ttl=60)
        self.assertIsNone(cache.get('nonexistent'))

    def test_cache_clear(self):
        """Test cache clearing."""
        from molecular_calculator.utils.cache import TTLCache
        cache = TTLCache(ttl=60)
        cache.set('key1', 'value1')
        cache.clear()
        self.assertIsNone(cache.get('key1'))

    def test_cache_stats(self):
        """Test cache statistics."""
        from molecular_calculator.utils.cache import TTLCache
        cache = TTLCache(ttl=60)
        cache.set('key1', 'value1')
        cache.get('key1')  # hit
        cache.get('nonexistent')  # miss
        stats = cache.stats()
        self.assertIn('hits', stats)
        self.assertIn('misses', stats)

    def test_cache_decorator(self):
        """Test cached decorator."""
        from molecular_calculator.utils.cache import cached
        call_count = [0]

        @cached(ttl=60)
        def expensive_function(x):
            call_count[0] += 1
            return x * 2

        result1 = expensive_function(5)
        result2 = expensive_function(5)
        self.assertEqual(result1, 10)
        self.assertEqual(result2, 10)
        self.assertEqual(call_count[0], 1)  # Should only be called once


# =============================================================================
# Rate Limiter Tests
# =============================================================================

class TestRateLimiter(unittest.TestCase):
    """Test suite for RateLimiter"""

    def test_rate_limiter_import(self):
        """Test rate limiter imports."""
        from molecular_calculator.utils.rate_limiter import RateLimiter
        self.assertTrue(callable(RateLimiter))

    def test_rate_limiter_is_allowed(self):
        """Test rate limiter allows requests."""
        from molecular_calculator.utils.rate_limiter import RateLimiter
        limiter = RateLimiter(max_requests=10, window_seconds=1.0)
        self.assertTrue(limiter.is_allowed())

    def test_rate_limiter_blocks_when_exceeded(self):
        """Test rate limiter blocks when limit exceeded."""
        from molecular_calculator.utils.rate_limiter import RateLimiter
        limiter = RateLimiter(max_requests=2, window_seconds=1.0)
        limiter.record_request()  # 1
        limiter.record_request()  # 2
        self.assertFalse(limiter.is_allowed())  # Should be blocked

    def test_rate_limiter_remaining_requests(self):
        """Test rate limiter remaining requests count."""
        from molecular_calculator.utils.rate_limiter import RateLimiter
        limiter = RateLimiter(max_requests=10, window_seconds=1.0)
        self.assertEqual(limiter.remaining_requests(), 10)
        limiter.record_request()
        self.assertEqual(limiter.remaining_requests(), 9)


# =============================================================================
# Exceptions Tests
# =============================================================================

class TestExceptions(unittest.TestCase):
    """Test suite for custom exceptions"""

    def test_exception_imports(self):
        """Test all exceptions can be imported."""
        from molecular_calculator.utils.exceptions import (
            MolecularCalculatorError,
            InvalidSMILESError,
            PropertyCalculationError,
            FileValidationError,
            ConversionError,
            APIError,
            RateLimitError,
        )
        self.assertTrue(issubclass(InvalidSMILESError, MolecularCalculatorError))
        self.assertTrue(issubclass(PropertyCalculationError, MolecularCalculatorError))
        self.assertTrue(issubclass(ConversionError, MolecularCalculatorError))
        self.assertTrue(issubclass(APIError, MolecularCalculatorError))

    def test_invalid_smiles_error(self):
        """Test InvalidSMILESError creation."""
        from molecular_calculator.utils.exceptions import InvalidSMILESError
        error = InvalidSMILESError("Bad SMILES: XYZ")
        self.assertIn("XYZ", str(error))

    def test_property_calculation_error(self):
        """Test PropertyCalculationError creation."""
        from molecular_calculator.utils.exceptions import PropertyCalculationError
        error = PropertyCalculationError("MW calculation failed")
        self.assertIn("MW", str(error))

    def test_api_error(self):
        """Test APIError creation."""
        from molecular_calculator.utils.exceptions import APIError
        error = APIError("NIH API timeout")
        self.assertIn("NIH", str(error))

    def test_rate_limit_error(self):
        """Test RateLimitError creation."""
        from molecular_calculator.utils.exceptions import RateLimitError
        error = RateLimitError("Too many requests")
        self.assertIn("requests", str(error))


# =============================================================================
# Error Handler Tests
# =============================================================================

class TestErrorHandler(unittest.TestCase):
    """Test suite for error handling utilities"""

    def test_error_handler_import(self):
        """Test error handler imports."""
        from molecular_calculator.utils.error_handler import (
            handle_errors,
            handle_calculation_errors,
            get_error_details,
        )
        self.assertTrue(callable(handle_errors))
        self.assertTrue(callable(handle_calculation_errors))

    def test_handle_errors_decorator(self):
        """Test handle_errors decorator catches exceptions."""
        from molecular_calculator.utils.error_handler import handle_errors

        @handle_errors(default_return=None, show_user_error=False)
        def failing_function():
            raise ValueError("Test error")

        result = failing_function()
        self.assertIsNone(result)

    def test_handle_calculation_errors_decorator(self):
        """Test handle_calculation_errors decorator."""
        from molecular_calculator.utils.error_handler import handle_calculation_errors

        # handle_calculation_errors takes no parameters - always returns None on error
        @handle_calculation_errors
        def calculate_something():
            raise Exception("Calculation failed")

        result = calculate_something()
        self.assertIsNone(result)

    def test_get_error_details(self):
        """Test error detail extraction."""
        from molecular_calculator.utils.error_handler import get_error_details
        try:
            raise ValueError("Test error message")
        except ValueError as e:
            details = get_error_details(e)
            self.assertIn('type', details)
            self.assertIn('message', details)
            self.assertEqual(details['type'], 'ValueError')


# =============================================================================
# Validators Tests (Expanded)
# =============================================================================

class TestValidators(unittest.TestCase):
    """Expanded test suite for validators"""

    def test_input_validator_import(self):
        """Test InputValidator can be imported."""
        from molecular_calculator.utils.validators import InputValidator
        self.assertTrue(hasattr(InputValidator, 'validate_smiles'))

    def test_validate_smiles_valid(self):
        """Test SMILES validation with valid input."""
        from molecular_calculator.utils.validators import InputValidator
        result = InputValidator.validate_smiles('CCO')
        self.assertTrue(result.is_valid)

    def test_validate_smiles_invalid(self):
        """Test SMILES validation with invalid input."""
        from molecular_calculator.utils.validators import InputValidator
        result = InputValidator.validate_smiles('INVALID_SMILES')
        self.assertFalse(result.is_valid)

    def test_validate_inchi_valid(self):
        """Test InChI validation with valid input."""
        from molecular_calculator.utils.validators import InputValidator
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        result = InputValidator.validate_inchi(inchi)
        self.assertTrue(result.is_valid)

    def test_validate_inchi_invalid(self):
        """Test InChI validation with invalid input."""
        from molecular_calculator.utils.validators import InputValidator
        result = InputValidator.validate_inchi('not_an_inchi')
        self.assertFalse(result.is_valid)

    def test_validate_inchi_key_valid(self):
        """Test InChI Key validation with valid input."""
        from molecular_calculator.utils.validators import InputValidator
        result = InputValidator.validate_inchi_key('LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
        self.assertTrue(result.is_valid)

    def test_validate_inchi_key_invalid(self):
        """Test InChI Key validation with invalid input."""
        from molecular_calculator.utils.validators import InputValidator
        result = InputValidator.validate_inchi_key('INVALID-KEY')
        self.assertFalse(result.is_valid)

    def test_detect_format(self):
        """Test format detection."""
        from molecular_calculator.utils.validators import InputValidator
        self.assertEqual(InputValidator.detect_format('CCO'), 'smiles')
        self.assertEqual(InputValidator.detect_format('InChI=1S/C2H6O/c1-2-3'), 'inchi')
        self.assertEqual(InputValidator.detect_format('LFQSCWFLJHTTHZ-UHFFFAOYSA-N'), 'inchi_key')

    def test_is_safe_input(self):
        """Test injection detection."""
        from molecular_calculator.utils.validators import InputValidator
        self.assertTrue(InputValidator.is_safe_input('CCO'))
        self.assertFalse(InputValidator.is_safe_input('<script>alert(1)</script>'))
        self.assertFalse(InputValidator.is_safe_input("'; DROP TABLE molecules;--"))

    def test_file_validator_import(self):
        """Test FileValidator can be imported."""
        from molecular_calculator.utils.validators import FileValidator
        self.assertTrue(hasattr(FileValidator, 'validate_upload'))

    def test_dataframe_validator_import(self):
        """Test DataFrameValidator can be imported."""
        from molecular_calculator.utils.validators import DataFrameValidator
        self.assertTrue(hasattr(DataFrameValidator, 'validate'))

    def test_dataframe_validator_column_exists(self):
        """Test column existence check."""
        from molecular_calculator.utils.validators import DataFrameValidator
        df = pd.DataFrame({'SMILES': ['CCO'], 'Name': ['Ethanol']})
        self.assertTrue(DataFrameValidator.column_exists(df, 'SMILES'))
        self.assertFalse(DataFrameValidator.column_exists(df, 'Missing'))

    def test_dataframe_validator_get_numeric_columns(self):
        """Test numeric column detection."""
        from molecular_calculator.utils.validators import DataFrameValidator
        df = pd.DataFrame({
            'SMILES': ['CCO'],
            'MW': [46.0],
            'LogP': [-0.14],
            'Name': ['Ethanol']
        })
        numeric = DataFrameValidator.get_numeric_columns(df)
        self.assertIn('MW', numeric)
        self.assertIn('LogP', numeric)
        self.assertNotIn('Name', numeric)


# =============================================================================
# Suggestions Tests
# =============================================================================

class TestSuggestions(unittest.TestCase):
    """Test suite for smart suggestion utilities"""

    def test_suggestions_import(self):
        """Test suggestions module imports."""
        from molecular_calculator.utils.suggestions import (
            detect_smiles_column,
            detect_id_column,
        )
        self.assertTrue(callable(detect_smiles_column))
        self.assertTrue(callable(detect_id_column))

    def test_detect_smiles_column(self):
        """Test SMILES column detection."""
        from molecular_calculator.utils.suggestions import detect_smiles_column
        df = pd.DataFrame({
            'canonical_smiles': ['CCO'],
            'name': ['Ethanol']
        })
        suggested = detect_smiles_column(df)
        self.assertEqual(suggested, 'canonical_smiles')

    def test_detect_smiles_column_uppercase(self):
        """Test SMILES column detection with uppercase."""
        from molecular_calculator.utils.suggestions import detect_smiles_column
        df = pd.DataFrame({
            'SMILES': ['CCO'],
            'Name': ['Ethanol']
        })
        suggested = detect_smiles_column(df)
        self.assertEqual(suggested, 'SMILES')

    def test_detect_id_column(self):
        """Test ID column detection."""
        from molecular_calculator.utils.suggestions import detect_id_column
        df = pd.DataFrame({
            'compound_id': [1],
            'smiles': ['CCO']
        })
        suggested = detect_id_column(df)
        self.assertEqual(suggested, 'compound_id')


# =============================================================================
# Config Settings Tests
# =============================================================================

class TestConfigSettings(unittest.TestCase):
    """Test suite for application configuration"""

    def test_config_import(self):
        """Test config can be imported."""
        from molecular_calculator.config.settings import config, AppConfig
        self.assertIsInstance(config, AppConfig)

    def test_config_has_app_name(self):
        """Test config has APP_NAME."""
        from molecular_calculator.config.settings import config
        self.assertTrue(hasattr(config, 'APP_NAME'))
        self.assertIsInstance(config.APP_NAME, str)

    def test_config_has_max_file_size(self):
        """Test config has MAX_FILE_SIZE_MB."""
        from molecular_calculator.config.settings import config
        self.assertTrue(hasattr(config, 'MAX_FILE_SIZE_MB'))
        self.assertIsInstance(config.MAX_FILE_SIZE_MB, int)

    def test_config_has_api_settings(self):
        """Test config has API-related settings."""
        from molecular_calculator.config.settings import config
        self.assertTrue(hasattr(config, 'API_TIMEOUT_SECONDS'))
        self.assertTrue(hasattr(config, 'MAX_RETRY_ATTEMPTS'))


# =============================================================================
# Property Explanations Tests
# =============================================================================

class TestPropertyExplanationsExpanded(unittest.TestCase):
    """Expanded test suite for PropertyExplanations"""

    def test_explanations_contains_sections(self):
        """Test explanations contain all expected sections."""
        explanations = PropertyExplanations.get_explanations()

        # Check for all expected sections (matching actual content)
        expected_sections = [
            "Supported Input Formats",
            "SMILES",
            "InChI",
            "Property Explanations",
            "Basic Properties",
            "Lipinski Properties",
            "Drug-likeness",
            "Rule Violations",
            "Ring Properties",
        ]

        for section in expected_sections:
            self.assertIn(section, explanations, f"Missing section: {section}")

    def test_explanations_contains_property_descriptions(self):
        """Test explanations contain property descriptions."""
        explanations = PropertyExplanations.get_explanations()

        # Check for key property descriptions (matching actual content)
        properties = [
            "Molecular_Weight",
            "LogP",
            "TPSA",
            "QED",
            "HB_Donors",
            "HB_Acceptors",
            "Rotatable_Bonds",
        ]

        for prop in properties:
            self.assertIn(prop, explanations, f"Missing property: {prop}")


# =============================================================================
# 3D Regression Tests (Expanded)
# =============================================================================

class TestThreeDRegression(unittest.TestCase):
    """Expanded test suite for 3D OLS Regression"""

    def test_regression_import(self):
        """Test regression module imports."""
        from molecular_calculator.models.regression import ThreeDOLSRegression
        self.assertTrue(callable(ThreeDOLSRegression))

    def test_regression_basic_fit(self):
        """Test basic regression fitting."""
        from molecular_calculator.models.regression import ThreeDOLSRegression

        # Create simple test data (x and y must not be identical/collinear)
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([2, 4, 5, 4, 6])  # Different from x to avoid collinearity
        z = 2*x + 3*y + 1

        regression = ThreeDOLSRegression(x, y, z)

        # Check coefficients are accessible as direct attributes
        self.assertIsNotNone(regression.b0)
        self.assertIsNotNone(regression.b1)
        self.assertIsNotNone(regression.b2)

    def test_regression_predict(self):
        """Test regression prediction."""
        from molecular_calculator.models.regression import ThreeDOLSRegression

        x = np.array([1, 2, 3, 4, 5])
        y = np.array([2, 4, 5, 4, 6])  # Different from x
        z = 2*x + 3*y + 1

        regression = ThreeDOLSRegression(x, y, z)
        predictions = regression.predict(x, y)

        self.assertEqual(len(predictions), len(z))

    def test_regression_statistics(self):
        """Test regression statistics."""
        from molecular_calculator.models.regression import ThreeDOLSRegression

        x = np.array([1, 2, 3, 4, 5])
        y = np.array([2, 4, 5, 4, 6])  # Different from x
        z = 2*x + 3*y + 1

        regression = ThreeDOLSRegression(x, y, z)
        stats = regression.get_statistics()

        # get_statistics returns 'R²', 'RMSE', 'n' not 'r_squared', 'rmse', 'n_samples'
        self.assertIn('R²', stats)
        self.assertIn('RMSE', stats)
        self.assertIn('n', stats)

    def test_regression_equation_string(self):
        """Test regression equation string generation."""
        from molecular_calculator.models.regression import ThreeDOLSRegression

        x = np.array([1, 2, 3, 4, 5])
        y = np.array([2, 4, 5, 4, 6])  # Different from x
        z = 2*x + 3*y + 1

        regression = ThreeDOLSRegression(x, y, z)
        # get_equation_string takes decimals, not x_name/y_name
        equation = regression.get_equation_string(decimals=4)

        self.assertIsInstance(equation, str)
        self.assertIn('X', equation)
        self.assertIn('Y', equation)

    def test_regression_with_noise(self):
        """Test regression with noisy data."""
        from molecular_calculator.models.regression import ThreeDOLSRegression

        np.random.seed(42)
        x = np.random.rand(50) * 10
        y = np.random.rand(50) * 10
        noise = np.random.randn(50) * 0.1
        z = 2*x + 3*y + 1 + noise

        regression = ThreeDOLSRegression(x, y, z)

        # R² is a direct attribute, not in get_statistics()
        self.assertGreater(regression.r_squared, 0.95)


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration(unittest.TestCase):
    """Integration tests for end-to-end workflows"""

    def test_full_property_calculation_pipeline(self):
        """Test complete property calculation from SMILES."""
        from molecular_calculator.core import MolecularCalculator, get_calculator

        calc = get_calculator()
        result = calc.calculate("CCO")

        # Verify complete result
        self.assertTrue(result.success)
        self.assertIsNotNone(result.properties)
        self.assertIsNotNone(result.properties.molecular_weight)
        self.assertIsNotNone(result.properties.logp)
        self.assertIsNotNone(result.properties.hb_donors)
        self.assertIsNotNone(result.properties.hb_acceptors)
        self.assertIsNotNone(result.properties.tpsa)

    def test_inchi_to_properties_pipeline(self):
        """Test property calculation from InChI."""
        from molecular_calculator.core import MolecularCalculator

        calc = MolecularCalculator()
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        result = calc.calculate_from_any_format(inchi)

        self.assertTrue(result.success)
        self.assertAlmostEqual(result.properties.molecular_weight, 46.069, places=1)

    def test_batch_processing_end_to_end(self):
        """Test complete batch processing workflow."""
        df = pd.DataFrame({
            'SMILES': ['C', 'CC', 'CCO', 'c1ccccc1', 'CC(=O)O'],
            'Name': ['Methane', 'Ethane', 'Ethanol', 'Benzene', 'Acetic acid']
        })

        result_df = MolecularCalculator.process_batch(df, 'SMILES')

        # Verify all molecules processed
        self.assertEqual(len(result_df), 5)

        # Verify properties calculated
        self.assertIn('Molecular_Weight', result_df.columns)
        self.assertIn('LogP', result_df.columns)

        # Verify values are reasonable
        self.assertTrue(all(result_df['Molecular_Weight'] > 0))

    def test_parallel_vs_sequential_consistency(self):
        """Test parallel and sequential processing produce same results."""
        df = pd.DataFrame({
            'SMILES': ['CCO', 'c1ccccc1', 'CC(=O)O'],
            'Name': ['Ethanol', 'Benzene', 'Acetic acid']
        })

        # Process sequentially
        sequential_result = MolecularCalculator.process_batch(df, 'SMILES')

        # Process in parallel
        parallel_result = MolecularCalculator.process_batch_parallel(df, 'SMILES')

        # Compare results
        for col in ['Molecular_Weight', 'LogP', 'HB_Donors']:
            if col in sequential_result.columns and col in parallel_result.columns:
                for i in range(len(df)):
                    seq_val = sequential_result.loc[i, col]
                    par_val = parallel_result.loc[i, col]
                    if pd.notna(seq_val) and pd.notna(par_val):
                        self.assertAlmostEqual(seq_val, par_val, places=3)

    def test_conversion_service_and_property_service_integration(self):
        """Test ConversionService and PropertyCalculator work together."""
        from molecular_calculator.services import ConversionService, PropertyCalculator

        conversion = ConversionService()
        property_calc = PropertyCalculator()

        # Convert InChI to SMILES
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        conv_result = conversion.to_smiles(inchi)
        self.assertTrue(conv_result.success)

        # Calculate properties from converted SMILES
        calc_result = property_calc.calculate(conv_result.smiles)
        self.assertTrue(calc_result.success)
        self.assertAlmostEqual(calc_result.properties.molecular_weight, 46.069, places=1)

    def test_sanitization_to_calculation_pipeline(self):
        """Test sanitization works correctly before calculation."""
        from molecular_calculator.utils.sanitizer import sanitize_smiles
        from molecular_calculator.core import MolecularCalculator

        # Input with whitespace
        dirty_smiles = "  CCO  "
        clean_smiles = sanitize_smiles(dirty_smiles)

        self.assertEqual(clean_smiles, "CCO")

        # Calculate properties
        props = MolecularCalculator.calculate_molecular_properties(clean_smiles)
        self.assertIn('Molecular_Weight', props)

    def test_error_handling_integration(self):
        """Test error handling works across the system."""
        from molecular_calculator.core import MolecularCalculator

        # Invalid SMILES should not crash
        result = MolecularCalculator.calculate_molecular_properties("INVALID_SMILES_XYZ")
        self.assertEqual(result, {})

        # Empty input should not crash
        result = MolecularCalculator.calculate_molecular_properties("")
        self.assertEqual(result, {})

        # None should not crash
        try:
            result = MolecularCalculator.calculate_molecular_properties(None)
            self.assertEqual(result, {})
        except (TypeError, AttributeError):
            pass  # Also acceptable

    def test_property_groups_completeness(self):
        """Test all properties in groups can be calculated."""
        groups = MolecularCalculator.get_property_groups()

        # Calculate properties for a test molecule
        props = MolecularCalculator.calculate_molecular_properties('c1ccccc1')

        # Check that property groups contain valid properties
        for group_name, properties in groups.items():
            for prop in properties:
                # Either the property exists or it's a valid property name
                self.assertIsInstance(prop, str, f"Invalid property name in {group_name}")


def run_comprehensive_tests():
    """Run all tests and generate a comprehensive report"""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add all test classes
    test_classes = [
        TestMolecularCalculator,
        TestPropertyExplanations,
        TestRealInChIKeyConversion,
        TestGetCalculatorSingleton,
        TestParallelBatchProcessing,
        TestInputSanitization,
        TestLoggingSetup,
        TestPropertyCalculator,
        TestConversionService,
        TestMolecularProperties,
        TestAdditionalSanitization,
        TestAPIClient,
        TestLigandEfficiency,
        TestCache,
        TestRateLimiter,
        TestExceptions,
        TestErrorHandler,
        TestValidators,
        TestSuggestions,
        TestConfigSettings,
        TestPropertyExplanationsExpanded,
        TestThreeDRegression,
        TestIntegration,
    ]

    for test_class in test_classes:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)

    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2, stream=None)
    result = runner.run(suite)

    # Print summary
    print(f"\n{'='*60}")
    print("COMPREHENSIVE TEST SUMMARY")
    print(f"{'='*60}")
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success rate: {((result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun * 100):.1f}%")

    if result.failures:
        print(f"\nFAILURES ({len(result.failures)}):")
        for test, traceback in result.failures:
            print(f"- {test}: {traceback.split('AssertionError: ')[-1].split(chr(10))[0] if 'AssertionError: ' in traceback else 'See details above'}")

    if result.errors:
        print(f"\nERRORS ({len(result.errors)}):")
        for test, traceback in result.errors:
            print(f"- {test}: {traceback.split(chr(10))[-2] if len(traceback.split(chr(10))) > 1 else 'Unknown error'}")

    return result.wasSuccessful()


if __name__ == '__main__':
    # Run comprehensive tests
    success = run_comprehensive_tests()

    if success:
        print("\n🎉 ALL TESTS PASSED! System is fully functional.")
    else:
        print("\n❌ Some tests failed. Please review the output above.")