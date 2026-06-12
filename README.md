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

**Ligand Efficiency Indices (LEI)**
- 8 indices based on AtlasCBS methodology: NSEI, NBEI, BEI, SEI, nBEI, mBEI, LEH, LEP
- Requires activity data (pKi or pIC50)
- Auto-detects existing columns (MW, TPSA, heavy atoms) to avoid redundant calculations

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

## GMM Analysis

Group your molecules into natural clusters using a Gaussian Mixture Model. Upload a dataset, pick the numeric columns to analyze, and the app discovers how your compounds separate into groups -- useful for spotting chemical series, property clusters, or outliers across a library.

**Two modes:**
- **Single property** -- overlays the discovered groups on a density curve, showing how one property (e.g. molecular weight) splits into subpopulations.
- **Multiple properties** -- clusters molecules across several properties at once and assigns each one to a group.

**How it works for you:**
- **Automatic group count** -- the app picks the optimal number of groups using BIC, so there's nothing to tune. Want manual control? A slider lets you set it yourself, and a BIC/AIC plot shows how each choice scores.
- **Works on any numeric data** -- run it on columns you already have, or let the app calculate molecular properties (MW, LogP, TPSA, QED, and more) from a SMILES column first.
- **Preprocessing built in** -- optional standardization and log-transforms handle skewed or differently-scaled properties.
- **Reproducible** -- a random seed keeps results identical across runs; reset or reshuffle it whenever you want.

**What you get:**
- A per-group summary table in real units, so "Group 1" is something you can actually read
- Each molecule's assigned group plus a confidence score, with low-confidence members flagged as potential outliers
- A plain-language interpretation of what the groups represent
- The fully labeled dataset as a downloadable CSV

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

## Assay Interference Detection

Screens compounds for known assay interference mechanisms across five categories:

| Flag | What It Catches | How |
|------|----------------|-----|
| **PAINS** | Pan-assay interference substructures | RDKit FilterCatalog (480 published SMARTS patterns from Baell & Holloway 2010) |
| **Aggregator** | Colloidal aggregation risk | Shoichet Lab heuristics (cLogP, MW, aromatic ring thresholds) |
| **Thiol-reactive** | Electrophilic compounds that react with cysteine residues | 14 SMARTS patterns for Michael acceptors, acylating agents, SN2 electrophiles |
| **Redox-active** | Redox cycling compounds (quinones, catechols, etc.) | 10 SMARTS patterns for redox-active functional groups |
| **Autofluorescent** | Fluorophore scaffolds that interfere with fluorescence assays | 13 SMARTS patterns for coumarins, xanthenes, anthracenes, stilbenes, etc. |

These are **flags for human review**, not automatic exclusions. A PAINS hit doesn't mean the compound is bad -- it means you should look more carefully at your assay data. Each flag includes the matched pattern name and a risk level so you can prioritize your review.

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

## Architecture

The codebase uses a modular package structure under `molecular_calculator/`:

- **`core/`** -- Core calculation engine (property computation, SMILES handling)
- **`services/`** -- Business logic (assay interference detection, ligand efficiency, API clients)
- **`models/`** -- Data models, regression, and property explanations
- **`ui/`** -- Streamlit pages and reusable UI components (charts, viewers, file uploaders)
- **`utils/`** -- Shared utilities (validation, caching, error handling, rate limiting)
- **`config/`** -- App settings and logging configuration

The entry point is `app.py`. If you're contributing or extending the calculator, this is the layout you'll be working with.

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
