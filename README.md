# 🧪 Molecular Properties Calculator

A comprehensive Streamlit web application for calculating chemical and molecular properties from molecular structures (SMILES, InChI, InChI Key).

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](YOUR_STREAMLIT_CLOUD_URL_HERE)

## 🚀 Features

### 📝 Input Formats
- **SMILES** (Simplified Molecular Input Line Entry System)
- **InChI** (International Chemical Identifier)
- **InChI Key** (Hashed InChI - limited conversion support)

### 📊 File Support
- **CSV** files
- **XLSX** (Excel) files
- Automatic format detection
- Flexible column naming (SMILES, Smiles, smiles, etc.)

### 🔬 Calculated Properties

#### Basic Properties
- Molecular Weight (Da)
- Heavy Atom Count
- Total Atom Count
- Bond Count
- Formal Charge

#### Lipinski Rule of Five
- LogP (Partition Coefficient)
- Hydrogen Bond Donors
- Hydrogen Bond Acceptors
- Topological Polar Surface Area (TPSA)
- Rotatable Bonds

#### Drug-likeness & Compliance
- **QED Score** (Quantitative Estimate of Drug-likeness, 0-1 scale)
- **Lipinski Violations** (0=compliant, 1=violates Rule of Five)
- **Veber Violations** (0=compliant, 1=violates oral bioavailability rules)

#### Ring Properties
- Aromatic Rings
- Aliphatic Rings
- Saturated Rings
- Total Ring Count
- Heteroatoms

#### Molecular Complexity
- Bertz Complexity Index
- Chi Connectivity Indices (Chi0, Chi1)

#### Additional Descriptors
- Crippen LogP & Molar Refractivity
- Labute Approximate Surface Area

### 🎯 Usage Modes

#### 1. Single Molecule Analysis
- Input individual molecules
- Interactive property selection
- Real-time calculation
- Visualization with charts
- Export results as CSV

#### 2. Batch Processing
- Upload CSV/XLSX files
- Process hundreds/thousands of molecules
- Progress tracking
- Summary statistics
- Export enhanced datasets

## 🛠️ Installation & Setup

### Prerequisites
- Python 3.8+
- pip package manager

### Local Installation

1. **Clone the repository**
```bash
git clone https://github.com/YOUR_USERNAME/molecular-properties-calculator.git
cd molecular-properties-calculator
```

2. **Install dependencies**
```bash
pip install -r requirements.txt
```

3. **Run the application**
```bash
streamlit run molecular_properties_app.py
```

4. **Open in browser**
Navigate to `http://localhost:8501`

### 🌐 Streamlit Cloud Deployment

This app is optimized for Streamlit Cloud deployment:

1. Fork this repository
2. Connect to Streamlit Cloud
3. Deploy with one click
4. All dependencies are automatically handled

## 📋 Usage Guide

### Single Molecule Analysis

1. **Select "Single Molecule" mode**
2. **Enter molecular structure** (SMILES, InChI, or InChI Key)
3. **Choose properties** using checkboxes:
   - Use group checkboxes to select all properties in a category
   - Use individual checkboxes to select specific properties
   - Or check "Calculate All Properties"
4. **Click "Calculate Properties"**
5. **View results** in table and chart format
6. **Download results** as CSV

### Batch Processing

1. **Select "Batch Processing" mode**
2. **Upload your file** (CSV or XLSX)
3. **Verify SMILES column** is detected correctly
4. **Select properties** to calculate
5. **Click "Process Batch"**
6. **Monitor progress** with progress bar
7. **Review results** and summary statistics
8. **Download enhanced dataset**

### 📁 Input File Format

Your CSV/XLSX should contain molecular structures in one column. Common column names that are auto-detected:

- `SMILES`, `Smiles`, `smiles`
- `SMI`, `smi`
- `CANONICAL_SMILES`

**Example CSV:**
```csv
Name,SMILES,Activity
Aspirin,CC(=O)OC1=CC=CC=C1C(=O)O,Active
Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C,Active
```

## 🎨 Key Benefits

### ✅ User-Friendly
- **No default selections** - choose only what you need
- **Flexible property selection** - individual or group selection
- **Clear explanations** - detailed property descriptions in-app
- **Progress tracking** - real-time batch processing updates

### ✅ Professional
- **Software compatibility** - works with StarDrop, Pipeline Pilot, KNIME
- **Clean data output** - numeric values, no problematic text/boolean columns
- **Error handling** - graceful handling of invalid molecules
- **Performance optimized** - efficient batch processing

### ✅ Research-Ready
- **Comprehensive properties** - covers most common molecular descriptors
- **Rule compliance** - Lipinski and Veber rules built-in
- **Export capabilities** - CSV download for further analysis
- **Visualization** - interactive charts for data exploration

## 🔬 Property Explanations

### Rule Violations (Binary Indicators)

#### Lipinski Violations
- **0** = Passes all Lipinski Rule of Five criteria:
  - Molecular Weight ≤ 500 Da
  - LogP ≤ 5
  - H-bond Donors ≤ 5
  - H-bond Acceptors ≤ 10
- **1** = Violates at least one criterion

#### Veber Violations
- **0** = Passes oral bioavailability rules:
  - TPSA ≤ 140 Ų
  - Rotatable Bonds ≤ 10
- **1** = Violates at least one criterion

### QED Score
- Range: 0.0 to 1.0
- Higher values = more drug-like
- Based on molecular properties and their distributions in approved drugs
- Considers: MW, LogP, HBD, HBA, PSA, rotatable bonds, aromatic rings, alerts

## 📊 Example Use Cases

### Drug Discovery
- **Lead optimization** - calculate properties for compound series
- **Library filtering** - identify drug-like compounds
- **ADMET prediction** - assess absorption and permeability

### Chemical Analysis
- **Descriptor calculation** - generate features for QSAR models
- **Diversity analysis** - characterize chemical libraries
- **Property profiling** - understand structure-property relationships

### Data Enhancement
- **Missing properties** - add calculated descriptors to existing datasets
- **Format standardization** - ensure consistent property calculations
- **Quality control** - validate molecular structures and properties

## 🐛 Troubleshooting

### Common Issues

**File Upload Problems:**
- Ensure file is CSV or XLSX format
- Check that SMILES column exists
- Verify no special characters in column names

**Calculation Errors:**
- Invalid SMILES will result in empty property values
- Some complex descriptors may not calculate for all molecules
- Check molecule structure validity

**Performance:**
- Large files (>10K molecules) may take several minutes
- Consider splitting very large datasets
- Close other browser tabs for better performance

## 🤝 Contributing

We welcome contributions! Please feel free to:

1. Report bugs or request features via GitHub Issues
2. Submit pull requests for improvements
3. Share feedback and suggestions
4. Contribute additional molecular descriptors

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **RDKit** - Open-source cheminformatics toolkit
- **Streamlit** - Web app framework
- **Plotly** - Interactive visualizations
- **Pandas** - Data manipulation and analysis

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/YOUR_USERNAME/molecular-properties-calculator/issues)
- **Documentation**: See the in-app "Information & Property Explanations" section
- **Examples**: Check the `examples/` folder for sample files

---

**Made with ❤️ for the cheminformatics community**

---

## 🚀 Quick Start

```bash
# Install and run locally
pip install -r requirements.txt
streamlit run molecular_properties_app.py
```

Or try the [live demo](YOUR_STREAMLIT_CLOUD_URL_HERE) on Streamlit Cloud!