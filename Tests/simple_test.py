#!/usr/bin/env python3
"""
Simple test to verify core functionality without heavy dependencies
"""

import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_basic_functionality():
    """Test basic functionality without dependencies"""
    print("Testing basic Python functionality...")

    # Test imports
    print("Testing imports...")
    import re
    import requests
    print("Basic imports successful")

    # Test input format detection logic
    print("Testing input format detection...")

    def detect_input_format(input_text):
        """Simplified version of format detection"""
        if input_text.startswith('InChI='):
            return 'inchi'
        elif re.match(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$', input_text):
            return 'inchi_key'
        else:
            return 'smiles'

    test_cases = [
        ('CCO', 'smiles'),
        ('InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', 'inchi'),
        ('LFQSCWFLJHTTHZ-UHFFFAOYSA-N', 'inchi_key'),
    ]

    for input_text, expected in test_cases:
        result = detect_input_format(input_text)
        assert result == expected, f"{input_text[:20]}... -> {result} (expected {expected})"
        print(f"{input_text[:20]}... -> {result}")

    # Test online API availability
    print("Testing online API availability...")
    try:
        response = requests.get('https://httpbin.org/status/200', timeout=5)
        if response.status_code == 200:
            print("Internet connection available")
        else:
            print("Internet connection issues")
    except Exception as e:
        print(f"Cannot test online APIs: {e}")

    print("Basic functionality tests passed!")

def test_rdkit_availability():
    """Test if RDKit can be imported and used"""
    print("\nTesting RDKit availability...")

    from rdkit import Chem
    from rdkit.Chem import QED, Descriptors
    print("RDKit imports successful")

    # Test basic RDKit functionality
    smiles = 'CCO'
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, "RDKit molecule creation failed"

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    qed_score = QED.qed(mol)

    assert mw > 0, "Molecular weight should be positive"
    assert isinstance(logp, float), "LogP should be a float"
    assert 0 <= qed_score <= 1, "QED should be between 0 and 1"

    print(f"RDKit calculations work: MW={mw:.1f}, LogP={logp:.2f}, QED={qed_score:.3f}")

def test_molecular_calculator_import():
    """Test if our molecular calculator can be imported"""
    print("\nTesting MolecularCalculator import...")

    from molecular_calculator import MolecularCalculator
    print("MolecularCalculator import successful")

    # Test basic methods exist
    methods_to_check = [
        'suppress_rdkit_warnings',
        'detect_input_format',
        'convert_to_smiles',
        'calculate_molecular_properties',
        'get_property_groups'
    ]

    for method in methods_to_check:
        assert hasattr(MolecularCalculator, method), f"Method {method} missing"
        print(f"Method {method} exists")

def test_property_calculations():
    """Test actual property calculations if RDKit is available"""
    print("\nTesting property calculations...")

    from molecular_calculator import MolecularCalculator

    # Test molecules
    test_molecules = {
        'ethanol': 'CCO',
        'aspirin': 'CC(=O)OC1=CC=CC=C1C(=O)O',
        'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
    }

    for name, smiles in test_molecules.items():
        properties = MolecularCalculator.calculate_molecular_properties(smiles)
        assert properties is not None, f"{name}: No properties calculated"
        assert 'Molecular_Weight' in properties, f"{name}: Missing Molecular_Weight"
        assert 'LogP' in properties, f"{name}: Missing LogP"

        mw = properties.get('Molecular_Weight', 'N/A')
        logp = properties.get('LogP', 'N/A')
        hbd = properties.get('HB_Donors', 'N/A')
        print(f"{name}: MW={mw}, LogP={logp}, HBD={hbd}")

def test_online_inchi_key_conversion():
    """Test InChI Key conversion if network is available"""
    import pytest
    print("\nTesting InChI Key conversion...")

    from molecular_calculator import MolecularCalculator

    # Test with ethanol InChI Key
    inchi_key = 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'

    try:
        result = MolecularCalculator.convert_inchi_key_to_smiles(inchi_key, timeout=10)
    except Exception as e:
        pytest.skip(f"Network unavailable for InChI Key conversion: {e}")

    if result is None:
        pytest.skip("InChI Key conversion failed (network issue)")

    print(f"InChI Key conversion successful: {inchi_key} -> {result}")

    # Verify the result by calculating properties
    properties = MolecularCalculator.calculate_molecular_properties(result)
    assert properties is not None, "Converted SMILES is invalid"
    assert 'Molecular_Weight' in properties, "Missing Molecular_Weight from converted SMILES"

    print(f"Converted SMILES is valid: MW={properties.get('Molecular_Weight', 'N/A')}")

def run_all_tests():
    """Run all available tests (for standalone script execution)"""
    print("="*60)
    print("SIMPLE FUNCTIONALITY TEST SUITE")
    print("="*60)

    tests = [
        ("Basic Functionality", test_basic_functionality),
        ("RDKit Availability", test_rdkit_availability),
        ("MolecularCalculator Import", test_molecular_calculator_import),
        ("Property Calculations", test_property_calculations),
        ("InChI Key Conversion", test_online_inchi_key_conversion),
    ]

    passed = 0
    failed = 0
    total = len(tests)

    for test_name, test_func in tests:
        print(f"\n{'-'*40}")
        print(f"Running: {test_name}")
        print(f"{'-'*40}")

        try:
            test_func()
            passed += 1
            print(f"PASSED: {test_name}")
        except AssertionError as e:
            failed += 1
            print(f"FAILED: {test_name} - {e}")
        except Exception as e:
            failed += 1
            print(f"ERROR: {test_name} - {e}")

    print(f"\n{'='*60}")
    print(f"TEST RESULTS: {passed}/{total} tests passed ({passed/total*100:.1f}%)")
    print(f"{'='*60}")

    if failed == 0:
        print("ALL TESTS PASSED! System is fully functional.")
        return True
    else:
        print(f"{failed} test(s) failed. See details above.")
        return False

if __name__ == '__main__':
    success = run_all_tests()
    sys.exit(0 if success else 1)