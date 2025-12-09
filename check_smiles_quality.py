"""
SMILES Data Quality Checker

This utility helps identify invalid or problematic SMILES strings in your dataset.
Run this on your data files before processing to catch issues early.

Usage:
    python check_smiles_quality.py your_data.csv
"""

import pandas as pd
import sys
from rdkit import Chem
from rdkit import RDLogger

# Suppress RDKit warnings for cleaner output
RDLogger.DisableLog('rdApp.*')


def validate_smiles(smiles_string):
    """
    Validate a SMILES string and return status information.
    
    Returns:
        dict: Status information including validity, warnings, and suggestions
    """
    result = {
        'valid': False,
        'warnings': [],
        'errors': [],
        'suggestions': []
    }
    
    # Check for empty/null
    if not smiles_string or (isinstance(smiles_string, str) and smiles_string.strip() == ''):
        result['errors'].append('Empty or null SMILES')
        result['suggestions'].append('Remove row or provide valid SMILES')
        return result
    
    # Convert to string if needed
    smiles_str = str(smiles_string).strip()
    
    # Check length
    if len(smiles_str) > 1000:
        result['warnings'].append(f'Unusually long SMILES ({len(smiles_str)} characters)')
    
    # Check for suspicious characters
    suspicious_chars = [char for char in smiles_str if not char.isalnum() and char not in '()[]=#-+@/\\*.%']
    if suspicious_chars:
        result['warnings'].append(f'Contains unusual characters: {set(suspicious_chars)}')
    
    # Check for incomplete brackets
    if smiles_str.count('(') != smiles_str.count(')'):
        result['errors'].append('Unmatched parentheses')
        result['suggestions'].append('Check if SMILES is truncated or corrupted')
    
    if smiles_str.count('[') != smiles_str.count(']'):
        result['errors'].append('Unmatched square brackets')
        result['suggestions'].append('Check if SMILES is truncated or corrupted')
    
    # Try to parse with RDKit
    try:
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            result['errors'].append('RDKit cannot parse this SMILES')
            result['suggestions'].append('Verify SMILES syntax is correct')
        else:
            result['valid'] = True
            # Additional checks on valid molecule
            num_atoms = mol.GetNumAtoms()
            if num_atoms == 0:
                result['warnings'].append('Molecule has zero atoms')
            elif num_atoms > 500:
                result['warnings'].append(f'Very large molecule ({num_atoms} atoms)')
    except Exception as e:
        result['errors'].append(f'Parse error: {str(e)}')
        result['suggestions'].append('Check SMILES syntax')
    
    return result


def check_dataset(filepath, smiles_column=None):
    """
    Check all SMILES in a dataset file.
    
    Args:
        filepath: Path to CSV file
        smiles_column: Name of SMILES column (auto-detected if None)
    """
    print(f"\n{'='*70}")
    print(f"🔍 SMILES Data Quality Check")
    print(f"{'='*70}")
    print(f"File: {filepath}\n")
    
    # Read data
    try:
        df = pd.read_csv(filepath)
        print(f"✅ Loaded {len(df)} rows, {len(df.columns)} columns\n")
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        return
    
    # Detect SMILES column
    if smiles_column is None:
        smiles_candidates = ['SMILES', 'smiles', 'Smiles', 'SMILES_String', 'smiles_string']
        for col in smiles_candidates:
            if col in df.columns:
                smiles_column = col
                break
        
        if smiles_column is None:
            # Look for column with SMILES-like data
            for col in df.columns:
                if df[col].dtype == 'object':
                    sample = str(df[col].iloc[0]) if len(df) > 0 else ''
                    if any(char in sample for char in ['C', 'c', 'N', 'O', '(', ')']):
                        smiles_column = col
                        break
    
    if smiles_column is None:
        print("❌ Could not detect SMILES column. Please specify with --column argument")
        print(f"Available columns: {', '.join(df.columns)}")
        return
    
    print(f"📊 Analyzing SMILES column: '{smiles_column}'")
    print(f"{'='*70}\n")
    
    # Validate each SMILES
    issues = []
    valid_count = 0
    warning_count = 0
    error_count = 0
    
    for idx, row in df.iterrows():
        smiles = row[smiles_column]
        result = validate_smiles(smiles)
        
        if result['valid']:
            valid_count += 1
            if result['warnings']:
                warning_count += 1
        else:
            error_count += 1
        
        # Collect rows with issues
        if result['errors'] or result['warnings']:
            issues.append({
                'row': idx,
                'smiles': str(smiles)[:80] + ('...' if len(str(smiles)) > 80 else ''),
                'result': result
            })
    
    # Summary
    print(f"📈 Summary:")
    print(f"  ✅ Valid SMILES: {valid_count} ({valid_count/len(df)*100:.1f}%)")
    print(f"  ⚠️  With warnings: {warning_count} ({warning_count/len(df)*100:.1f}%)")
    print(f"  ❌ Invalid/Errors: {error_count} ({error_count/len(df)*100:.1f}%)")
    print(f"\n{'='*70}\n")
    
    # Show issues
    if issues:
        print(f"🔍 Detailed Issues (showing first 10):\n")
        for i, issue in enumerate(issues[:10]):
            print(f"Row {issue['row']}:")
            print(f"  SMILES: {issue['smiles']}")
            
            if issue['result']['errors']:
                for error in issue['result']['errors']:
                    print(f"  ❌ {error}")
            
            if issue['result']['warnings']:
                for warning in issue['result']['warnings']:
                    print(f"  ⚠️  {warning}")
            
            if issue['result']['suggestions']:
                print(f"  💡 Suggestions:")
                for suggestion in issue['result']['suggestions']:
                    print(f"     - {suggestion}")
            print()
        
        if len(issues) > 10:
            print(f"... and {len(issues) - 10} more issues\n")
    else:
        print("✅ No issues found! All SMILES are valid.\n")
    
    print(f"{'='*70}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python check_smiles_quality.py <csv_file> [smiles_column]")
        print("\nExample:")
        print("  python check_smiles_quality.py data.csv")
        print("  python check_smiles_quality.py data.csv SMILES")
        sys.exit(1)
    
    filepath = sys.argv[1]
    smiles_col = sys.argv[2] if len(sys.argv) > 2 else None
    
    check_dataset(filepath, smiles_col)

