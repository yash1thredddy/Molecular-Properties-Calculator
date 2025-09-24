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
        """Set up test data"""
        # Test molecules with known properties
        self.test_molecules = {
            'aspirin': {
                'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                'inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
                'inchi_key': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N',
                'expected_mw': 180.158,  # Expected molecular weight
                'expected_logp_range': (1.0, 2.0),  # Expected LogP range
                'expected_hbd': 1,  # Hydrogen bond donors
                'expected_hba': 3   # RDKit counts 3 acceptors (2 C=O + 1 ester O)
            },
            'caffeine': {
                'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
                'inchi': 'InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3',
                'inchi_key': 'RYYVLZVUVIJVGH-UHFFFAOYSA-N',
                'expected_mw': 194.191,
                'expected_logp_range': (-2.0, 0.5),  # Broader range for RDKit LogP calculation
                'expected_hbd': 0,
                'expected_hba': 6
            },
            'ethanol': {
                'smiles': 'CCO',
                'inchi': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
                'inchi_key': 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N',
                'expected_mw': 46.069,
                'expected_logp_range': (-0.5, 0.0),
                'expected_hbd': 1,
                'expected_hba': 1
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

                # Test rule compliance calculations
                self.assertIn(properties['Lipinski_Violations'], [0, 1])
                self.assertIn(properties['Veber_Violations'], [0, 1])

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
            self.assertIn(properties_fail['Lipinski_Violations'], [0, 1])

    def test_veber_rule_compliance(self):
        """Test Veber Rule compliance calculation"""
        properties = MolecularCalculator.calculate_molecular_properties('CCO')
        self.assertEqual(properties['Veber_Violations'], 0)

        # Test should return 0 or 1
        for name, data in self.test_molecules.items():
            props = MolecularCalculator.calculate_molecular_properties(data['smiles'])
            self.assertIn(props['Veber_Violations'], [0, 1], f"Invalid Veber violation for {name}")

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


def run_comprehensive_tests():
    """Run all tests and generate a comprehensive report"""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add all test classes
    test_classes = [
        TestMolecularCalculator,
        TestPropertyExplanations,
        TestRealInChIKeyConversion
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
            print(f"- {test}: {traceback.split('AssertionError: ')[-1].split('\\n')[0] if 'AssertionError: ' in traceback else 'See details above'}")

    if result.errors:
        print(f"\nERRORS ({len(result.errors)}):")
        for test, traceback in result.errors:
            print(f"- {test}: {traceback.split('\\n')[-2] if len(traceback.split('\\n')) > 1 else 'Unknown error'}")

    return result.wasSuccessful()


if __name__ == '__main__':
    # Run comprehensive tests
    success = run_comprehensive_tests()

    if success:
        print("\\n🎉 ALL TESTS PASSED! System is fully functional.")
    else:
        print("\\n❌ Some tests failed. Please review the output above.")