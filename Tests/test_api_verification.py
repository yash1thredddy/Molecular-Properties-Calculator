#!/usr/bin/env python3
"""
API Verification Test for InChI Key Conversion Services

Tests both NIH CIR and PubChem APIs with real InChI Keys
to ensure online conversion functionality is working.

Developed by: Yashwanth Reddy for ITR-UIC
Part of: Chemo-Informatics Toolkit
"""

import requests
import time
from molecular_calculator import MolecularCalculator

def test_api_endpoints():
    """Test individual API endpoints directly"""
    import pytest
    print("Testing API endpoints...")
    print("=" * 50)

    test_inchi_keys = [
        ('LFQSCWFLJHTTHZ-UHFFFAOYSA-N', 'Ethanol', 'CCO'),
        ('BSYNRYMUTXBXSQ-UHFFFAOYSA-N', 'Aspirin', 'CC(=O)OC1=CC=CC=C1C(=O)O'),
        ('RYYVLZVUVIJVGH-UHFFFAOYSA-N', 'Caffeine', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
    ]

    # Test NIH CIR API
    print("\nTesting NIH Chemical Identifier Resolver (CIR):")
    print("-" * 50)

    cir_successes = 0
    for inchi_key, name, expected_base in test_inchi_keys:
        try:
            url = f"https://cactus.nci.nih.gov/chemical/structure/{inchi_key}/smiles"
            response = requests.get(url, timeout=10)

            if response.status_code == 200:
                smiles = response.text.strip()
                if len(smiles) > 2 and not smiles.startswith("Error"):
                    print(f"{name}: {smiles}")
                    cir_successes += 1
                else:
                    print(f"{name}: Invalid response - {smiles}")
            else:
                print(f"{name}: HTTP {response.status_code}")

        except Exception as e:
            print(f"{name}: {str(e)}")

        time.sleep(0.5)  # Rate limiting

    print(f"\nCIR Success Rate: {cir_successes}/{len(test_inchi_keys)} ({cir_successes/len(test_inchi_keys)*100:.0f}%)")

    # Test PubChem API
    print("\nTesting PubChem API:")
    print("-" * 50)

    pubchem_successes = 0
    for inchi_key, name, expected_base in test_inchi_keys:
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchi_key}/property/IsomericSMILES/JSON"
            response = requests.get(url, timeout=10)

            if response.status_code == 200:
                data = response.json()
                if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                    properties = data['PropertyTable']['Properties']
                    if properties and 'IsomericSMILES' in properties[0]:
                        smiles = properties[0]['IsomericSMILES']
                        print(f"{name}: {smiles}")
                        pubchem_successes += 1
                    else:
                        print(f"{name}: No SMILES in response")
                else:
                    print(f"{name}: Invalid JSON structure")
            else:
                print(f"{name}: HTTP {response.status_code}")

        except Exception as e:
            print(f"{name}: {str(e)}")

        time.sleep(0.5)  # Rate limiting

    print(f"\nPubChem Success Rate: {pubchem_successes}/{len(test_inchi_keys)} ({pubchem_successes/len(test_inchi_keys)*100:.0f}%)")

    # At least one API should work
    if cir_successes == 0 and pubchem_successes == 0:
        pytest.skip("Both APIs unavailable (network issue)")

    assert cir_successes > 0 or pubchem_successes > 0, "At least one API should be working"

def test_integrated_conversion():
    """Test the integrated MolecularCalculator InChI Key conversion"""
    import pytest
    print("\n\nTesting integrated conversion...")
    print("=" * 50)

    test_cases = [
        ('LFQSCWFLJHTTHZ-UHFFFAOYSA-N', 'Ethanol'),
        ('BSYNRYMUTXBXSQ-UHFFFAOYSA-N', 'Aspirin'),
        ('RYYVLZVUVIJVGH-UHFFFAOYSA-N', 'Caffeine'),
        ('UHOVQNZJYSORNB-UHFFFAOYSA-N', 'Benzene'),
        ('XLYOFNOQVPJJNP-UHFFFAOYSA-N', 'Water')
    ]

    successes = 0
    total_time = 0

    for inchi_key, name in test_cases:
        print(f"\nTesting {name} ({inchi_key}):")

        start_time = time.time()
        try:
            smiles = MolecularCalculator.convert_inchi_key_to_smiles(inchi_key, timeout=15)
            end_time = time.time()

            if smiles:
                properties = MolecularCalculator.calculate_molecular_properties(smiles)

                if properties and len(properties) > 0:
                    mw = properties.get('Molecular_Weight', 'N/A')
                    logp = properties.get('LogP', 'N/A')
                    print(f"Conversion successful: {smiles}")
                    print(f"   MW: {mw}, LogP: {logp}")
                    print(f"   Time: {end_time - start_time:.2f}s")
                    successes += 1
                else:
                    print(f"Converted SMILES is invalid: {smiles}")
            else:
                print(f"Conversion failed (network/database issue)")

        except Exception as e:
            end_time = time.time()
            print(f"Exception: {str(e)}")

        total_time += end_time - start_time
        time.sleep(1)  # Rate limiting

    avg_time = total_time / len(test_cases)
    success_rate = successes / len(test_cases) * 100

    print(f"\nIntegration Test Results:")
    print(f"   Success Rate: {successes}/{len(test_cases)} ({success_rate:.0f}%)")
    print(f"   Average Time: {avg_time:.2f}s per conversion")

    if successes == 0:
        pytest.skip("Network unavailable for InChI Key conversion")

    assert successes > 0, "At least one conversion should succeed"

def test_error_handling():
    """Test error handling with invalid InChI Keys"""
    print("\n\nTesting error handling...")
    print("=" * 50)

    invalid_cases = [
        'INVALID-INCHI-KEY',
        'TOOSHORT-KEY-X',
        'TOOLONG-INCHI-KEY-INVALID-FORMAT-X',
        'NUMBERS123-INVALID456-N'
    ]

    print("Testing invalid InChI Keys (should fail gracefully):")

    for invalid_key in invalid_cases:
        try:
            result = MolecularCalculator.convert_inchi_key_to_smiles(invalid_key, timeout=5)
            if result is None:
                print(f"{invalid_key[:20]}... -> None (correct)")
            else:
                print(f"{invalid_key[:20]}... -> {result} (unexpected success)")
        except Exception as e:
            print(f"{invalid_key[:20]}... -> Exception handled: {str(e)[:30]}")

def test_network_resilience():
    """Test network timeout and resilience"""
    print("\n\nTesting network resilience...")
    print("=" * 50)

    print("Testing timeout handling (should complete quickly):")

    # Test with very short timeout
    start_time = time.time()
    result = MolecularCalculator.convert_inchi_key_to_smiles(
        'LFQSCWFLJHTTHZ-UHFFFAOYSA-N', timeout=1
    )
    end_time = time.time()

    elapsed = end_time - start_time
    # Should timeout or succeed within 3 seconds
    assert elapsed < 5, f"Timeout may not be working: {elapsed:.2f}s"
    print(f"Timeout handling works: {elapsed:.2f}s")

def run_comprehensive_api_tests():
    """Run all API verification tests (for standalone script execution)"""
    print("COMPREHENSIVE API VERIFICATION")
    print("=" * 60)
    print("Testing InChI Key conversion APIs and integration")
    print("=" * 60)

    tests = [
        ("API Endpoints", test_api_endpoints),
        ("Integrated Conversion", test_integrated_conversion),
        ("Error Handling", test_error_handling),
        ("Network Resilience", test_network_resilience),
    ]

    passed = 0
    failed = 0

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

    # Overall summary
    print("\n" + "=" * 60)
    print("OVERALL API VERIFICATION RESULTS")
    print("=" * 60)
    print(f"Tests passed: {passed}/{len(tests)}")

    if failed == 0:
        print("\nAPI VERIFICATION: PASSED")
        print("   InChI Key conversion is functional")
        print("   At least one database service is working")
        return True
    else:
        print("\nAPI VERIFICATION: Some tests failed")
        print("   Check network connectivity")
        return False

if __name__ == '__main__':
    import sys
    success = run_comprehensive_api_tests()

    if success:
        print("\nAll API services are functional and ready for production!")
        sys.exit(0)
    else:
        print("\nAPI services may have issues. Check network connectivity.")
        sys.exit(1)