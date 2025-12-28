# Molecular Properties Calculator

A web-based tool for calculating chemical and molecular properties from SMILES, InChI, and InChI Key structures. Built for researchers, medicinal chemists, and anyone working with molecular data.

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://molecular-properties-calculator-2.streamlit.app/)

## What It Does

Drop in a molecule (or thousands of them), pick the properties you need, and get results you can export. Works with single molecules for quick lookups or batch files when you're processing compound libraries.

**Supported input formats:**
- SMILES notation
- InChI strings
- InChI Keys (converted automatically via NIH/PubChem databases)

**File uploads:**
- CSV and Excel files
- Automatic column detection for common naming patterns

## Properties You Can Calculate

**Physical Properties**
- Molecular Weight, Atom Counts, Bond Counts, Formal Charge

**Lipinski Descriptors**
- LogP, Hydrogen Bond Donors/Acceptors, TPSA, Rotatable Bonds

**Drug-likeness**
- QED Score (0-1 scale measuring how drug-like a molecule is)
- Lipinski Rule of Five violations
- Veber oral bioavailability violations

**Ring Analysis**
- Aromatic, Aliphatic, and Saturated ring counts
- Heteroatom count

**Complexity Indices**
- Bertz Complexity, Chi Connectivity Indices
- Crippen descriptors, Labute ASA

**Ligand Efficiency**
- LE, LLE, LELP, SILE, and other efficiency metrics
- Requires activity data (pIC50, pKi, etc.)

## Getting Started

### Run Locally

```bash
git clone https://github.com/yash1threddy/molecular-properties-calculator.git
cd molecular-properties-calculator
pip install -r requirements.txt
streamlit run app.py
```

Open `http://localhost:8501` in your browser.

### Deploy to Streamlit Cloud

Fork this repo, connect it to Streamlit Cloud, and deploy. Dependencies install automatically.

## How to Use

### Single Molecule

1. Enter your SMILES, InChI, or InChI Key
2. Select which properties to calculate
3. View results and download as CSV

### Batch Processing

1. Upload a CSV or Excel file containing molecular structures
2. The app detects your SMILES column automatically
3. Select properties and hit process
4. Download results with all calculated properties appended

Your input file just needs a column with molecular structures. Common column names like `SMILES`, `smiles`, `SMI`, or `CANONICAL_SMILES` are detected automatically.

**Example input:**
```csv
Name,SMILES
Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O
Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
```

## 3D Regression Analysis

The app includes a 3D OLS regression tool for exploring relationships between molecular properties. Pick any three numeric properties as X, Y, and Z axes, and fit a regression plane through your data. Useful for SAR analysis and property optimization studies.

## Visualization

Interactive charts built with Plotly:
- Scatter plots, histograms, box plots, violin plots
- Correlation heatmaps and pair plots
- 3D scatter with regression surfaces
- Click any data point to view its 2D structure

## Rule Compliance

**Lipinski Rule of Five** — A molecule passes (violation = 0) if:
- Molecular Weight ≤ 500
- LogP ≤ 5
- H-bond Donors ≤ 5
- H-bond Acceptors ≤ 10

**Veber Rules** — For oral bioavailability (violation = 0) if:
- TPSA ≤ 140 Å²
- Rotatable Bonds ≤ 10

## Technical Notes

**InChI Key Conversion**
- Primary: NIH Chemical Identifier Resolver
- Fallback: PubChem REST API
- Requires internet connection (can be disabled in settings)

**Performance**
- Handles files with thousands of molecules
- Progress tracking for batch jobs
- Rate limiting on external API calls

**Compatibility**
- Output works with StarDrop, Pipeline Pilot, KNIME, and other cheminformatics tools
- Clean numeric values, no problematic text columns

## Requirements

- Python 3.8+
- RDKit for cheminformatics calculations
- Streamlit for the web interface
- See `requirements.txt` for full list

## Contributing

Found a bug or have a feature request? Open an issue on GitHub. Pull requests welcome.

## License

CC BY-NC-ND 4.0 (Creative Commons Attribution-NonCommercial-NoDerivatives 4.0)

You may share this software with attribution, but commercial use and modifications are not permitted without permission. See [LICENSE](LICENSE) for details.

## Acknowledgments

Built with RDKit, Streamlit, Plotly, and Pandas.

---

Developed by Yashwanth Reddy | ITR-UIC | Dr. Guido Pauli's Research Group
