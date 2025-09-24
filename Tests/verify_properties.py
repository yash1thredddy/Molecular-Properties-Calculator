#!/usr/bin/env python3
"""
Quick verification script to check actual RDKit values for test molecules
"""

from molecular_calculator import MolecularCalculator

def verify_test_molecules():
    """Check actual property values for test molecules"""

    test_molecules = {
        'ethanol': 'CCO',
        'aspirin': 'CC(=O)OC1=CC=CC=C1C(=O)O',
        'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
    }

    print("Actual RDKit Values for Test Molecules:")
    print("=" * 60)

    for name, smiles in test_molecules.items():
        print(f"\n{name.upper()} ({smiles})")
        print("-" * 40)

        properties = MolecularCalculator.calculate_molecular_properties(smiles)

        key_properties = [
            'Molecular_Weight', 'LogP', 'HB_Donors', 'HB_Acceptors',
            'TPSA', 'Rotatable_Bonds', 'QED'
        ]

        for prop in key_properties:
            if prop in properties:
                print(f"{prop}: {properties[prop]}")

        print(f"Total properties calculated: {len(properties)}")

def test_invalid_smiles():
    """Test invalid SMILES handling"""
    print("\n\nInvalid SMILES Testing:")
    print("=" * 60)

    invalid_cases = [
        'INVALID_SMILES_XYZ',
        '',
        'XYZ123',
        '[]',
        None
    ]

    for case in invalid_cases:
        try:
            result = MolecularCalculator.calculate_molecular_properties(case)
            print(f"'{case}' -> {len(result)} properties, empty dict: {result == {}}")
        except Exception as e:
            print(f"'{case}' -> Exception: {e}")

if __name__ == '__main__':
    verify_test_molecules()
    test_invalid_smiles()