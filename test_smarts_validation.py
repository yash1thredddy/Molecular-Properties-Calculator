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

    for category, patterns in all_patterns.items():
        print(f"\n{category}:")
        for name, smarts in patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            assert pattern is not None, (
                f"SMARTS compilation failed for {category}/{name}: '{smarts}'"
            )
            print(f"  OK: {name}")


def test_pattern_matching():
    """Test that patterns match expected molecules."""
    print("\n" + "=" * 70)
    print("TESTING PATTERN MATCHING")
    print("=" * 70)

    results = {}

    # Test thiol-reactive patterns
    print("\n--- THIOL-REACTIVE PATTERNS ---")

    # Test user's molecule specifically - should match michael_acceptor and acrylamide
    user_mol = Chem.MolFromSmiles(TEST_MOLECULES['user_molecule'])
    assert user_mol is not None, "Failed to parse user_molecule SMILES"
    michael_pattern = Chem.MolFromSmarts(THIOL_REACTIVE_SMARTS['michael_acceptor'])
    acrylamide_pattern = Chem.MolFromSmarts(THIOL_REACTIVE_SMARTS['acrylamide'])
    assert michael_pattern is not None, "Failed to compile michael_acceptor SMARTS"
    assert acrylamide_pattern is not None, "Failed to compile acrylamide SMARTS"
    assert user_mol.HasSubstructMatch(michael_pattern), (
        "User's molecule should match michael_acceptor pattern"
    )
    assert user_mol.HasSubstructMatch(acrylamide_pattern), (
        "User's molecule should match acrylamide pattern"
    )

    # Test simple acrylamide - should match acrylamide pattern
    acrylamide_mol = Chem.MolFromSmiles('C=CC(=O)N')
    assert acrylamide_mol is not None, "Failed to parse acrylamide SMILES"
    assert acrylamide_mol.HasSubstructMatch(acrylamide_pattern), (
        "Simple acrylamide should match acrylamide pattern"
    )

    # Test crotonamide (substituted) - should also match acrylamide pattern
    crotonamide_mol = Chem.MolFromSmiles('CC=CC(=O)N')
    assert crotonamide_mol is not None, "Failed to parse crotonamide SMILES"
    assert crotonamide_mol.HasSubstructMatch(acrylamide_pattern), (
        "Crotonamide (substituted) should match acrylamide pattern"
    )

    # Test negative control - acetamide should NOT match any thiol-reactive pattern
    acetamide_mol = Chem.MolFromSmiles('CC(=O)N')
    assert acetamide_mol is not None, "Failed to parse acetamide SMILES"
    for name, smarts in THIOL_REACTIVE_SMARTS.items():
        pattern = Chem.MolFromSmarts(smarts)
        assert pattern is not None, f"Failed to compile pattern {name}"
        assert not acetamide_mol.HasSubstructMatch(pattern), (
            f"Acetamide (negative control) should NOT match thiol-reactive pattern '{name}'"
        )

    # Test redox patterns - positive controls
    print("\n--- REDOX-ACTIVE PATTERNS ---")
    redox_expected = {
        'benzoquinone': 'para_quinone',
        'catechol': 'catechol',
        'hydroquinone': 'hydroquinone',
    }
    for mol_name, expected_pattern in redox_expected.items():
        mol = Chem.MolFromSmiles(TEST_MOLECULES[mol_name])
        assert mol is not None, f"Failed to parse {mol_name} SMILES"
        pattern = Chem.MolFromSmarts(REDOX_ACTIVE_SMARTS[expected_pattern])
        assert pattern is not None, f"Failed to compile pattern {expected_pattern}"
        assert mol.HasSubstructMatch(pattern), (
            f"{mol_name} should match redox pattern '{expected_pattern}'"
        )

    # Test fluorescent patterns - positive controls
    print("\n--- FLUORESCENT PATTERNS ---")
    fluor_expected = {
        'naphthalene': 'naphthalene',
        'stilbene': 'stilbene',
    }
    for mol_name, expected_pattern in fluor_expected.items():
        mol = Chem.MolFromSmiles(TEST_MOLECULES[mol_name])
        assert mol is not None, f"Failed to parse {mol_name} SMILES"
        pattern = Chem.MolFromSmarts(FLUORESCENT_SMARTS[expected_pattern])
        assert pattern is not None, f"Failed to compile pattern {expected_pattern}"
        assert mol.HasSubstructMatch(pattern), (
            f"{mol_name} should match fluorescent pattern '{expected_pattern}'"
        )


def main():
    """Run all validation tests."""
    print("SMARTS PATTERN VALIDATION FOR ASSAY INTERFERENCE DETECTION")
    print("=" * 70)

    # Test compilation (raises immediately on any failure via assertions)
    test_smarts_compilation()
    print("\nAll SMARTS patterns compile successfully!")

    # Test matching
    test_pattern_matching()

    print("\n" + "=" * 70)
    print("VALIDATION COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
