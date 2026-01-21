"""Property explanations and documentation.

Provides detailed explanations of all molecular properties calculated by the application.
"""


class PropertyExplanations:
    """Class containing property explanations and documentation."""

    @staticmethod
    def get_explanations() -> str:
        """
        Get detailed explanations of all molecular properties.

        Returns:
            Markdown formatted string with property explanations
        """
        return """
        ### Supported Input Formats
        - **SMILES**: Simplified Molecular Input Line Entry System
        - **InChI**: International Chemical Identifier
        - **InChI Key**: Hashed version of InChI (database lookup required)

        ### File Formats
        - **CSV**: Comma-separated values
        - **XLSX**: Excel spreadsheet format

        ### Property Explanations

        #### Basic Properties
        - **Molecular_Weight**: Molecular weight in Daltons (Da)
        - **Heavy_Atom_Count**: Number of non-hydrogen atoms
        - **Atom_Count**: Total number of atoms (including hydrogens)
        - **Bond_Count**: Total number of bonds
        - **Formal_Charge**: Net formal charge of the molecule

        #### Lipinski Properties (Rule of Five)
        - **LogP**: Partition coefficient (lipophilicity, -2 to 6 typical range)
        - **HB_Donors**: Hydrogen bond donors (≤5 for drug-likeness)
        - **HB_Acceptors**: Hydrogen bond acceptors (≤10 for drug-likeness)
        - **TPSA**: Topological polar surface area in Ų (≤140 for oral bioavailability)
        - **10xPSA_MW**: PSA/MW ratio scaled by 10 (lower values indicate better membrane permeability)
        - **NPOLoNHA**: Polar atoms / Heavy atoms ratio (NPOL/NHA, indicates molecular polarity fraction)
        - **Rotatable_Bonds**: Number of rotatable bonds (≤10 for oral bioavailability)

        #### Drug-likeness
        - **QED**: Quantitative Estimate of Drug-likeness (0-1 scale, higher is better)

        #### Rule Violations (Binary: 0=Compliant, 1=Violates)
        - **Lipinski_Violations**:
          - 0 = Passes Lipinski Rule of Five (MW≤500, LogP≤5, HBD≤5, HBA≤10)
          - 1 = Violates at least one Lipinski rule
        - **Veber_Violations**:
          - 0 = Passes Veber Rule (TPSA≤140, RotBonds≤10)
          - 1 = Violates at least one Veber rule

        #### Ring Properties
        - **Aromatic_Rings**: Number of aromatic rings
        - **Aliphatic_Rings**: Number of non-aromatic rings
        - **Saturated_Rings**: Number of saturated rings
        - **Ring_Count**: Total number of rings
        - **Heteroatoms**: Number of non-carbon, non-hydrogen atoms

        #### Complexity
        - **BertzCT**: Bertz complexity index (higher = more complex)
        - **Chi0**: Chi connectivity index 0 (molecular connectivity)
        - **Chi1**: Chi connectivity index 1 (molecular connectivity)

        #### Additional Descriptors
        - **CrippenLogP**: Crippen's LogP calculation method
        - **CrippenMR**: Crippen's molar refractivity
        - **LabuteASA**: Labute's approximate surface area

        ### Usage Tips
        - No properties are selected by default - choose what you need
        - Perfect for files with existing calculations - add only missing properties
        - Common SMILES column names are automatically detected
        - Invalid molecules will result in empty property values
        - All calculations use RDKit library
        - Results are compatible with StarDrop, Pipeline Pilot, and other software

        ### InChI Key Conversion
        - **Enabled by default** for convenience
        - Uses **NIH Chemical Identifier Resolver** (primary)
        - **PubChem API** fallback for redundancy
        - Can be disabled in Settings for privacy/offline use
        - Requires internet connection for database lookup

        ### Ligand Efficiency Indices (LEI) - Requires pKi Values
        Based on AtlasCBS methodology (Cele Abad-Zapatero, A. Cortes-Cabrera 2013)

        #### Available LEIs:
        - **NSEI**: Normalized Surface Efficiency Index (pKi / Polar Atoms)
        - **NBEI**: Normalized Binding Efficiency Index (pKi / Heavy Atoms)
        - **BEI**: Binding Efficiency Index (pKi / (MW/1000))
        - **SEI**: Surface Efficiency Index (pKi / (TPSA/100))
        - **nBEI**: Alternative Binding Efficiency Index (-log10(Ki / Heavy Atoms))
        - **mBEI**: Molecular Binding Efficiency Index (-log10(Ki / MW))
        - **LEH**: Ligand Efficiency Hopkins (-ΔG / Heavy Atoms) [Hopkins et al. DDT, 2004]
        - **LEP**: Ligand Efficiency Polar (-ΔG / Polar Atoms)

        #### LEI Requirements:
        - **pKi column required** in uploaded file
        - Optional columns (will be calculated from SMILES if missing):
          - MW (Molecular Weight)
          - TPSA (Topological Polar Surface Area)
          - Heavy_Atom_Count
        - System automatically detects available columns
        - Only calculates missing values from SMILES when needed

        ### Assay Interference Flags
        Detection of compounds with known assay interference mechanisms.
        Based on Bisson et al. (2016) and Baell & Holloway (2010).

        #### Interference Mechanisms Detected:
        - **PAINS**: Pan-Assay Interference Substructures - compounds that show promiscuous activity
        - **Aggregator**: Colloidal aggregation risk - compounds that form aggregates causing false positives
        - **Redox**: Redox-active groups (catechols, quinones, thiols) that interfere with redox-based assays
        - **Fluorescence**: Autofluorescent scaffolds that interfere with fluorescence-based assays
        - **Thiol**: Thiol-reactive electrophiles (Michael acceptors, maleimides) that covalently modify proteins

        #### Interpretation:
        - Flags indicate potential assay interference, not necessarily invalid activity
        - Compounds with structural evidence (e.g., PDB structures) may have genuine activity despite flags
        - Use as orthogonal information alongside other drug-likeness metrics

        ### Development Information
        Developed by: **Yashwanth Reddy** for **ITR-UIC**
        Part of: **Chemo-Informatics Toolkit**
        """

    @staticmethod
    def get_property_groups() -> dict:
        """
        Get property groups for UI display.

        Returns:
            Dictionary mapping group names to property lists
        """
        return {
            "Basic Properties": [
                "Molecular_Weight",
                "Heavy_Atom_Count",
                "Atom_Count",
                "Bond_Count",
                "Formal_Charge"
            ],
            "Lipinski Properties": [
                "LogP",
                "HB_Donors",
                "HB_Acceptors",
                "TPSA",
                "NPOLoNHA",
                "Rotatable_Bonds"
            ],
            "Drug-likeness": [
                "QED",
                "10xPSA_MW"
            ],
            "Rule Violations": [
                "Lipinski_Violations",
                "Veber_Violations"
            ],
            "Ring Properties": [
                "Aromatic_Rings",
                "Aliphatic_Rings",
                "Saturated_Rings",
                "Ring_Count",
                "Heteroatoms"
            ],
            "Complexity": [
                "BertzCT",
                "Chi0",
                "Chi1"
            ],
            "Additional Descriptors": [
                "CrippenLogP",
                "CrippenMR",
                "LabuteASA"
            ]
        }
