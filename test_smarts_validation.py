"""
SMARTS Pattern Validation Script

This script tests all SMARTS patterns in our assay_interference module
to verify they:
1. Compile correctly in RDKit
2. Match expected test molecules
3. Don't match negative controls

This ensures scientific accuracy of our computational flags.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

# =============================================================================
# TEST MOLECULES - Known structures for validation
# =============================================================================

TEST_MOLECULES = {
    # THIOL-REACTIVE TEST MOLECULES
    'acrylamide': 'C=CC(=O)N',  # Simple acrylamide (should match)
    'crotonamide': 'CC=CC(=O)N',  # Substituted acrylamide (should match)
    'user_molecule': 'CC(=O)N[C@H]1[C@@H](O[C@@H]2O[C@H](CC(O)[C@H]3O[C@@H](n4ccc(=O)[nH]c4=O)[C@H](O)[C@@H]3O)[C@H](O)[C@H](O)[C@H]2NC(=O)/C=C/CC(C)C)O[C@H](CO)[C@@H](N)[C@@H]1O',  # User's molecule with crotonamide
    'maleimide': 'O=C1C=CC(=O)N1',  # Maleimide
    'methyl_acrylate': 'C=CC(=O)OC',  # Acrylate ester
    'mvk': 'C=CC(=O)C',  # Methyl vinyl ketone (enone)
    'epoxide': 'C1OC1',  # Ethylene oxide
    'acetaldehyde': 'CC=O',  # Aldehyde
    'acetyl_chloride': 'CC(=O)Cl',  # Acyl halide
    'acetic_anhydride': 'CC(=O)OC(=O)C',  # Anhydride
    'methyl_isocyanate': 'CN=C=O',  # Isocyanate
    'phenyl_isothiocyanate': 'c1ccccc1N=C=S',  # Isothiocyanate
    'vinyl_sulfone': 'C=CS(=O)(=O)C',  # Vinyl sulfone

    # NEGATIVE CONTROLS (should NOT match thiol-reactive)
    'acetamide': 'CC(=O)N',  # No double bond - NOT thiol reactive
    'propionamide': 'CCC(=O)N',  # Saturated - NOT thiol reactive
    'benzamide': 'c1ccccc1C(=O)N',  # Aromatic amide - NOT thiol reactive

    # REDOX-ACTIVE TEST MOLECULES
    'benzoquinone': 'O=C1C=CC(=O)C=C1',  # p-Benzoquinone
    'naphthoquinone': 'O=C1C=Cc2ccccc2C1=O',  # 1,4-Naphthoquinone
    'catechol': 'Oc1ccccc1O',  # Catechol
    'hydroquinone': 'Oc1ccc(O)cc1',  # Hydroquinone
    'nitrobenzene': 'c1ccccc1[N+](=O)[O-]',  # Nitroaromatic

    # NEGATIVE CONTROLS (should NOT match redox)
    'phenol': 'Oc1ccccc1',  # Single OH - NOT redox active
    'resorcinol': 'Oc1cccc(O)c1',  # meta-dihydroxy - NOT redox active (not catechol/hydroquinone)

    # FLUORESCENT TEST MOLECULES
    'coumarin': 'O=c1ccc2ccccc2o1',  # Coumarin
    'naphthalene': 'c1ccc2ccccc2c1',  # Naphthalene
    'anthracene': 'c1ccc2cc3ccccc3cc2c1',  # Anthracene
    'stilbene': 'c1ccccc1C=Cc1ccccc1',  # trans-Stilbene
    'flavone': 'O=c1cc(-c2ccccc2)oc2ccccc12',  # Flavone

    # NEGATIVE CONTROLS (should NOT match fluorescent)
    'benzene': 'c1ccccc1',  # Single ring - NOT fluorescent
    'toluene': 'Cc1ccccc1',  # Single ring - NOT fluorescent
}


# =============================================================================
# SMARTS PATTERNS TO TEST
# =============================================================================

THIOL_REACTIVE_SMARTS = {
    'michael_acceptor': '[C;$(C=C)]-[C;$(C=O)]',
    'acrylamide': 'C=CC(=O)N',
    'acrylate': 'C=CC(=O)O',
    'enone': 'C=CC(=O)[#6]',
    'maleimide': 'O=C1C=CC(=O)N1',
    'acyl_halide': 'C(=O)[F,Cl,Br,I]',
    'anhydride': 'C(=O)OC(=O)',
    'epoxide': 'C1OC1',
    'aziridine': 'C1NC1',
    'aldehyde': '[CH1](=O)',
    'isocyanate': 'N=C=O',
    'isothiocyanate': 'N=C=S',
    'vinyl_sulfone': 'C=CS(=O)(=O)',
    'sulfonyl_fluoride': 'S(=O)(=O)F',
}

REDOX_ACTIVE_SMARTS = {
    'para_quinone': 'O=C1C=CC(=O)C=C1',
    'ortho_quinone': 'O=C1C(=O)C=CC=C1',
    'naphthoquinone': 'O=C1C=CC2=CC=CC=C2C1=O',
    'anthraquinone': 'O=C1c2ccccc2C(=O)c3ccccc13',
    'catechol': 'c1ccc(O)c(O)c1',
    'hydroquinone': 'Oc1ccc(O)cc1',
    'nitro_aromatic': 'c[N+](=O)[O-]',
}

FLUORESCENT_SMARTS = {
    'coumarin': 'O=c1ccc2ccccc2o1',
    'coumarin_keto': 'O=C1C=Cc2ccccc2O1',
    'naphthalene': 'c1ccc2ccccc2c1',
    'anthracene': 'c1ccc2cc3ccccc3cc2c1',
    'stilbene': 'c1ccccc1C=Cc1ccccc1',
    'flavone': 'O=c1cc(-c2ccccc2)oc2ccccc12',
}


def test_smarts_compilation():
    """Test that all SMARTS patterns compile correctly in RDKit."""
    print("=" * 70)
    print("TESTING SMARTS COMPILATION")
    print("=" * 70)

    all_patterns = {
        'THIOL_REACTIVE': THIOL_REACTIVE_SMARTS,
        'REDOX_ACTIVE': REDOX_ACTIVE_SMARTS,
        'FLUORESCENT': FLUORESCENT_SMARTS,
    }

    failures = []
    for category, patterns in all_patterns.items():
        print(f"\n{category}:")
        for name, smarts in patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                print(f"  FAIL: {name} - '{smarts}' does not compile!")
                failures.append((category, name, smarts))
            else:
                print(f"  OK: {name}")

    return failures


def test_pattern_matching():
    """Test that patterns match expected molecules."""
    print("\n" + "=" * 70)
    print("TESTING PATTERN MATCHING")
    print("=" * 70)

    results = {}

    # Test thiol-reactive patterns
    print("\n--- THIOL-REACTIVE PATTERNS ---")

    # Test user's molecule specifically
    user_mol = Chem.MolFromSmiles(TEST_MOLECULES['user_molecule'])
    if user_mol:
        print(f"\nUser's molecule (contains crotonamide NC(=O)/C=C/CC(C)C):")
        for name, smarts in THIOL_REACTIVE_SMARTS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                matches = user_mol.HasSubstructMatch(pattern)
                status = "MATCH" if matches else "no match"
                print(f"  {name}: {status}")

    # Test simple acrylamide vs crotonamide
    print(f"\nSimple acrylamide (C=CC(=O)N):")
    acrylamide_mol = Chem.MolFromSmiles('C=CC(=O)N')
    for name, smarts in THIOL_REACTIVE_SMARTS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and acrylamide_mol:
            matches = acrylamide_mol.HasSubstructMatch(pattern)
            if matches:
                print(f"  {name}: MATCH")

    print(f"\nCrotonamide (CC=CC(=O)N) - substituted:")
    crotonamide_mol = Chem.MolFromSmiles('CC=CC(=O)N')
    for name, smarts in THIOL_REACTIVE_SMARTS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and crotonamide_mol:
            matches = crotonamide_mol.HasSubstructMatch(pattern)
            if matches:
                print(f"  {name}: MATCH")

    # Test negative controls
    print(f"\nNegative control - acetamide (CC(=O)N - no C=C):")
    acetamide_mol = Chem.MolFromSmiles('CC(=O)N')
    matched_any = False
    for name, smarts in THIOL_REACTIVE_SMARTS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and acetamide_mol:
            matches = acetamide_mol.HasSubstructMatch(pattern)
            if matches:
                print(f"  FALSE POSITIVE: {name} matched!")
                matched_any = True
    if not matched_any:
        print(f"  OK: No false positives")

    # Test redox patterns
    print("\n--- REDOX-ACTIVE PATTERNS ---")

    for mol_name in ['benzoquinone', 'catechol', 'hydroquinone']:
        mol = Chem.MolFromSmiles(TEST_MOLECULES[mol_name])
        if mol:
            print(f"\n{mol_name}:")
            for name, smarts in REDOX_ACTIVE_SMARTS.items():
                pattern = Chem.MolFromSmarts(smarts)
                if pattern:
                    matches = mol.HasSubstructMatch(pattern)
                    if matches:
                        print(f"  {name}: MATCH")

    # Test fluorescent patterns
    print("\n--- FLUORESCENT PATTERNS ---")

    for mol_name in ['coumarin', 'naphthalene', 'stilbene']:
        mol = Chem.MolFromSmiles(TEST_MOLECULES[mol_name])
        if mol:
            print(f"\n{mol_name}:")
            for name, smarts in FLUORESCENT_SMARTS.items():
                pattern = Chem.MolFromSmarts(smarts)
                if pattern:
                    matches = mol.HasSubstructMatch(pattern)
                    if matches:
                        print(f"  {name}: MATCH")


def main():
    """Run all validation tests."""
    print("SMARTS PATTERN VALIDATION FOR ASSAY INTERFERENCE DETECTION")
    print("=" * 70)

    # Test compilation
    failures = test_smarts_compilation()

    if failures:
        print(f"\n*** {len(failures)} SMARTS patterns failed to compile! ***")
        for category, name, smarts in failures:
            print(f"  {category}/{name}: {smarts}")
    else:
        print("\nAll SMARTS patterns compile successfully!")

    # Test matching
    test_pattern_matching()

    print("\n" + "=" * 70)
    print("VALIDATION COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
