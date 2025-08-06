# HeartMAP: Heart Multi-chamber Analysis Platform

[![PyPI version](https://badge.fury.io/py/heartmap.svg)](https://badge.fury.io/py/heartmap)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![CI Status](https://github.com/Tumo505/HeartMap/workflows/CI/badge.svg)](https://github.com/Tumo505/HeartMap/actions)

> **A comprehensive Python package for analyzing cell-cell communication across all four chambers of the human heart using single-cell RNA sequencing data.**

## 🫀 What is HeartMAP?

HeartMAP (Heart Multi-chamber Analysis Platform) is a specialized bioinformatics package designed to decode the complex cellular interactions within the human heart. Unlike general single-cell analysis tools, HeartMAP is purpose-built for cardiac biology, offering chamber-specific insights that are crucial for understanding heart function, disease, and therapeutic opportunities.

### 🎯 Why HeartMAP?

- **🏥 Clinical Relevance**: Designed for translational cardiac research with direct therapeutic implications
- **🔬 Scientific Innovation**: First comprehensive framework for multi-chamber cardiac communication analysis
- **🚀 Production Ready**: Fully tested, documented, and deployed package with multiple interfaces
- **📊 Accessible**: Works on standard hardware (8GB+ RAM) with optimized computational efficiency
- **🧪 Validated**: Tested on real human heart datasets with reproducible results

## 📦 Installation

### Quick Install

```bash
# Install from PyPI
pip install heartmap

# Install with all optional features
pip install heartmap[all]

# Install specific feature sets
pip install heartmap[communication]  # Cell communication analysis
pip install heartmap[api]            # REST API and web interface
pip install heartmap[dev]            # Development tools
```

### Verify Installation

```bash
# Test the installation
python -c "import heartmap; print('✅ HeartMAP installed successfully!')"

# Check available pipelines
python -c "from heartmap import BasicPipeline, AdvancedCommunicationPipeline, MultiChamberPipeline; print('✅ All pipelines available!')"
```

## 🚀 Quick Start

### 1. Command Line Interface (Simplest)

```bash
# Analyze your heart data with one command
heartmap your_data.h5ad

# Comprehensive analysis with custom output
heartmap your_data.h5ad --analysis-type comprehensive --output-dir results/
```

### 2. Python API (Most Flexible)

```python
from heartmap import Config
from heartmap.pipelines import ComprehensivePipeline

# Load configuration (or use defaults)
config = Config.default()

# Create and run analysis pipeline
pipeline = ComprehensivePipeline(config)
results = pipeline.run('your_data.h5ad', 'results/')

print("Analysis complete! Check the 'results/' directory.")
```

### 3. Web Interface (Most User-Friendly)

```bash
# Start the web interface
python -c "from heartmap.api import run_web_interface; run_web_interface()"
# Open http://localhost:7860 in your browser
```

## 🔬 What Can HeartMAP Do?

### Core Analysis Pipelines

#### 1. **Basic Pipeline** - Foundation Analysis
```python
from heartmap.pipelines import BasicPipeline

pipeline = BasicPipeline(config)
results = pipeline.run('data.h5ad')
```

**What it does:**
- Quality control and data preprocessing
- Cell type identification and clustering
- Basic visualization (UMAP, t-SNE)
- Fundamental cardiac cell annotation

**Use when:** You need reliable cell type identification and quality metrics

#### 2. **Advanced Communication Pipeline** - Cellular Interactions
```python
from heartmap.pipelines import AdvancedCommunicationPipeline

pipeline = AdvancedCommunicationPipeline(config)
results = pipeline.run('data.h5ad')
```

**What it does:**
- Ligand-receptor interaction analysis
- Communication hub identification
- Pathway enrichment analysis
- Temporal communication dynamics

**Use when:** You want to understand how heart cells communicate

#### 3. **Multi-Chamber Pipeline** - Spatial Analysis
```python
from heartmap.pipelines import MultiChamberPipeline

pipeline = MultiChamberPipeline(config)
results = pipeline.run('data.h5ad')
```

**What it does:**
- Chamber-specific marker identification
- Cross-chamber correlation analysis
- Chamber-specific pathway analysis
- Comparative chamber biology

**Use when:** You need chamber-specific insights (RA, RV, LA, LV)

#### 4. **Comprehensive Pipeline** - Complete Analysis
```python
from heartmap.pipelines import ComprehensivePipeline

pipeline = ComprehensivePipeline(config)
results = pipeline.run('data.h5ad')
```

**What it does:**
- All of the above in one integrated workflow
- Comprehensive visualizations and reports
- Complete multi-chamber communication atlas

**Use when:** You want the full HeartMAP analysis experience

## 📊 Real-World Applications

### 🏥 Clinical Research

```python
# Identify chamber-specific therapeutic targets
from heartmap import Config
from heartmap.pipelines import MultiChamberPipeline

config = Config.default()
config.analysis.focus_chambers = ['LV', 'RV']  # Focus on ventricles

pipeline = MultiChamberPipeline(config)
results = pipeline.run('patient_data.h5ad')

# Extract chamber-specific markers
lv_markers = results['chamber_markers']['LV']
rv_markers = results['chamber_markers']['RV']
```

### 🔬 Drug Discovery

```python
# Analyze communication pathways for drug targets
from heartmap.pipelines import AdvancedCommunicationPipeline

pipeline = AdvancedCommunicationPipeline(config)
results = pipeline.run('heart_disease_data.h5ad')

# Identify communication hubs for targeting
hubs = results['communication_hubs']
pathways = results['pathway_enrichment']
```

### 📚 Educational Research

```python
# Compare healthy vs diseased heart patterns
healthy_results = pipeline.run('healthy_heart.h5ad')
diseased_results = pipeline.run('diseased_heart.h5ad')

# Export for publications
from heartmap.utils import ResultsExporter
exporter = ResultsExporter(config)
exporter.export_publication_figures(healthy_results, 'figures/')
```

## ⚙️ Configuration

HeartMAP is highly configurable to work with different hardware and research needs:

### Memory Optimization

```python
from heartmap import Config

# For 8GB RAM systems
config = Config.default()
config.data.max_cells_subset = 30000
config.data.max_genes_subset = 3000

# For 16GB+ RAM systems
config.data.max_cells_subset = 50000
config.data.max_genes_subset = 5000

# For testing/development
config.data.test_mode = True
config.data.max_cells_subset = 5000
```

### Analysis Parameters

```python
# Customize analysis parameters
config.analysis.resolution = 0.8  # Clustering resolution
config.analysis.n_neighbors = 15  # UMAP neighbors
config.analysis.min_genes = 200   # Quality control

# Chamber-specific analysis
config.analysis.focus_chambers = ['LV', 'LA']  # Left heart only
config.analysis.chamber_comparison = True      # Enable comparisons
```

### Output Configuration

```python
# Customize outputs
config.output.save_figures = True
config.output.figure_format = 'png'
config.output.save_intermediate = False
config.output.generate_report = True
```

## 🎯 Use Cases by Research Domain

### 💊 Pharmaceutical Research
- **Drug Target Discovery**: Identify chamber-specific therapeutic targets
- **Safety Assessment**: Understand chamber-specific drug effects
- **Biomarker Development**: Discover chamber-specific disease markers

### 🏥 Clinical Cardiology
- **Precision Medicine**: Chamber-specific treatment strategies
- **Disease Mechanism**: Understand chamber-specific pathology
- **Patient Stratification**: Classify patients by communication patterns

### 🔬 Basic Research
- **Cardiac Development**: Understand chamber formation and maturation
- **Evolutionary Biology**: Compare cardiac communication across species
- **Systems Biology**: Model cardiac cellular networks

### 📊 Computational Biology
- **Method Development**: Benchmark new single-cell methods
- **Data Integration**: Combine multiple cardiac datasets
- **Network Analysis**: Study cellular communication networks

## 📈 Performance & Scalability

### Tested Configurations

| Hardware | Dataset Size | Memory Usage | Runtime | Status |
|----------|-------------|--------------|---------|---------|
| 8GB RAM | 30K cells | ~6GB | 15 min | ✅ Recommended |
| 16GB RAM | 50K cells | ~12GB | 25 min | ✅ Optimal |
| 32GB RAM | 100K cells | ~24GB | 45 min | ✅ High-throughput |

### Optimization Tips

```python
# For large datasets
config.data.use_highly_variable_genes = True
config.data.highly_variable_top_n = 3000

# For faster analysis
config.analysis.skip_advanced_qc = True
config.visualization.create_plots = False

# For memory efficiency
config.data.process_in_chunks = True
config.data.chunk_size = 10000
```

## 🔌 API Interfaces

### REST API

```bash
# Start the API server
python -m heartmap.api.rest

# Analyze data via HTTP
curl -X POST "http://localhost:8000/analyze" \
     -F "file=@data.h5ad" \
     -F "analysis_type=comprehensive"
```

### Web Interface

```python
from heartmap.api import run_web_interface

# Start Gradio web interface
run_web_interface(port=7860, share=True)
```

### CLI Tools

```bash
# Check available analyses
heartmap --list-analyses

# Get help for specific analysis
heartmap --help comprehensive

# Validate your data format
heartmap --validate data.h5ad
```

## 📊 Output Structure

HeartMAP generates comprehensive results in organized directories:

```
results/
├── figures/                    # All visualizations
│   ├── umap_clusters.png      # Basic clustering
│   ├── communication_hubs.png # Communication analysis
│   └── chamber_heatmap.png    # Chamber-specific patterns
├── tables/                     # All data tables
│   ├── chamber_markers.csv    # Marker genes per chamber
│   ├── communication_pairs.csv # Cell-cell interactions
│   └── pathway_enrichment.csv # Enriched pathways
├── processed_data/             # Processed datasets
│   └── heartmap_annotated.h5ad # Fully annotated data
└── reports/                    # Analysis reports
    ├── summary_report.html     # Interactive HTML report
    └── methods_report.md       # Methods and parameters
```

## 🧪 Validation & Testing

### Test Your Installation

```python
# Quick functionality test
from heartmap.tests import run_package_tests
run_package_tests()

# Test with demo data
from heartmap.demo import run_demo_analysis
run_demo_analysis()
```

### Validate Your Data

```python
from heartmap.data import DataValidator

validator = DataValidator()
is_valid, issues = validator.validate_h5ad('your_data.h5ad')

if not is_valid:
    print("Data issues found:", issues)
else:
    print("✅ Data format is compatible with HeartMAP")
```

## 🤝 Community & Support

### Getting Help

- **📖 Documentation**: [HeartMAP Wiki](https://github.com/Tumo505/HeartMap/wiki)
- **💬 Discussions**: [GitHub Discussions](https://github.com/Tumo505/HeartMap/discussions)
- **🐛 Issues**: [GitHub Issues](https://github.com/Tumo505/HeartMap/issues)
- **📧 Email**: 28346416@mylife.unisa.ac.za

### Contributing

```bash
# Development setup
git clone https://github.com/Tumo505/HeartMap.git
cd HeartMap
pip install -e .[dev]

# Run tests before contributing
python -m pytest tests/
python -m flake8 src/heartmap/
python -m mypy src/heartmap/
```

## 📚 Citation

If you use HeartMAP in your research, please cite:

```bibtex
@software{heartmap2025,
  title={HeartMAP: Heart Multi-chamber Analysis Platform},
  author={Kgabeng, Tumo and Wang, Lulu and Ngwangwa, Harry and Pandelani, Thanyani},
  year={2025},
  url={https://github.com/Tumo505/HeartMap},
  version={1.0.0}
}
```

## 📄 License

HeartMAP is released under the Apache 2.0 License. See [LICENSE](LICENSE) for details.

---

## 🎉 Get Started Today!

```bash
# Install HeartMAP
pip install heartmap

# Analyze your heart data
heartmap your_data.h5ad

# Explore results
ls results/
```

**HeartMAP: Unlocking the secrets of cardiac cellular communication, one chamber at a time.** 🫀✨
