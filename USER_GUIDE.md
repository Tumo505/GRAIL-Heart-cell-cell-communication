# HeartMAP User Guide

> **Complete guide to using the HeartMAP package for cardiac single-cell analysis**

## 📚 Table of Contents

1. [Getting Started](#getting-started)
2. [Understanding Your Data](#understanding-your-data)
3. [Choosing the Right Analysis](#choosing-the-right-analysis)
4. [Step-by-Step Tutorials](#step-by-step-tutorials)
5. [Interpreting Results](#interpreting-results)
6. [Advanced Usage](#advanced-usage)
7. [Troubleshooting](#troubleshooting)

---

## 🚀 Getting Started

### Installation

```bash
# Install HeartMAP
pip install heartmap

# Verify installation
python -c "import heartmap; print('✅ HeartMAP ready!')"
```

### Your First Analysis

```bash
# Download example data (replace with your data)
# heartmap-data --download example

# Run basic analysis
heartmap your_heart_data.h5ad

# Check results
ls results/
```

**That's it!** HeartMAP will automatically:
- Quality control your data
- Identify cell types
- Analyze chamber patterns
- Generate visualizations
- Create a summary report

---

## 🔬 Understanding Your Data

### Data Requirements

HeartMAP works with single-cell RNA-seq data in AnnData format (`.h5ad` files).

**Required:**
- Gene expression matrix (cells × genes)
- Cell and gene annotations

**Optional but Helpful:**
- Chamber information (`obs['chamber']`)
- Cell type annotations (`obs['cell_type']`)
- Quality metrics (`obs['n_genes']`, `obs['n_counts']`)

### Data Format Check

```python
import scanpy as sc

# Load your data
adata = sc.read_h5ad('your_data.h5ad')

# Check data structure
print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")
print(f"Available metadata: {list(adata.obs.columns)}")

# Check for chamber information
if 'chamber' in adata.obs.columns:
    print(f"Chambers: {adata.obs['chamber'].unique()}")
else:
    print("No chamber information found")
```

### Preparing Your Data

If your data needs preparation:

```python
import pandas as pd
import scanpy as sc

# Load data
adata = sc.read_h5ad('raw_data.h5ad')

# Add chamber information if missing
# (Replace with your chamber assignment logic)
chamber_mapping = {
    'Atrial_cells': 'RA',  # Right Atrium
    'Ventricular_cells': 'LV',  # Left Ventricle
    # ... add your mappings
}

if 'chamber' not in adata.obs.columns:
    # Example: assign based on cell type
    adata.obs['chamber'] = adata.obs['cell_type'].map(chamber_mapping)
    
# Save prepared data
adata.write('prepared_data.h5ad')
```

---

## 🎯 Choosing the Right Analysis

### Analysis Types Overview

| Analysis Type | When to Use | What You Get | Time |
|---------------|-------------|--------------|------|
| **Basic** | First look at data, quality check | Cell types, basic QC, clustering | 5-10 min |
| **Communication** | Understand cell interactions | Cell-cell communication networks | 10-15 min |
| **Multi-Chamber** | Chamber-specific insights | Chamber markers, comparisons | 15-20 min |
| **Comprehensive** | Complete HeartMAP analysis | Everything above + integrated report | 20-30 min |

### Decision Tree

```
Do you need chamber-specific insights?
├─ YES: Multi-Chamber or Comprehensive
└─ NO: Basic or Communication

Is this your first analysis of this dataset?
├─ YES: Start with Basic
└─ NO: Choose based on research question

Do you have limited time/resources?
├─ YES: Basic
└─ NO: Comprehensive (recommended)

Are you interested in cell-cell interactions?
├─ YES: Communication or Comprehensive
└─ NO: Basic or Multi-Chamber
```

---

## 📖 Step-by-Step Tutorials

### Tutorial 1: First-Time User (Basic Analysis)

**Goal:** Understand your heart data basics

**Step 1: Run Analysis**
```bash
heartmap your_data.h5ad --analysis-type basic --output-dir results/basic/
```

**Step 2: Check Results**
```bash
ls results/basic/
# figures/umap_clusters.png    - Cell clustering
# figures/quality_metrics.png  - Data quality
# tables/cell_annotations.csv  - Cell type assignments
# report.html                  - Summary report
```

**Step 3: Interpret Results**
```python
import pandas as pd
import scanpy as sc

# Load annotated data
adata = sc.read_h5ad('results/basic/annotated_data.h5ad')

# Check cell type distribution
cell_types = adata.obs['cell_type'].value_counts()
print("Cell types found:")
print(cell_types)

# Visualize
sc.pl.umap(adata, color='cell_type', save='_cell_types.pdf')
```

### Tutorial 2: Understanding Communication (Communication Analysis)

**Goal:** Discover how heart cells communicate

**Step 1: Run Communication Analysis**
```python
from heartmap import Config
from heartmap.pipelines import AdvancedCommunicationPipeline

config = Config.default()
pipeline = AdvancedCommunicationPipeline(config)
results = pipeline.run('your_data.h5ad', 'results/communication/')
```

**Step 2: Explore Communication Results**
```python
# Load communication results
interactions = pd.read_csv('results/communication/tables/interactions.csv')
hubs = pd.read_csv('results/communication/tables/communication_hubs.csv')

# Top interacting cell type pairs
top_interactions = interactions.nlargest(10, 'interaction_strength')
print("Top cell-cell interactions:")
print(top_interactions[['source_type', 'target_type', 'ligand', 'receptor']])

# Communication hub cells
print("\nTop communication hubs:")
print(hubs.nlargest(10, 'hub_score')[['cell_type', 'hub_score']])
```

**Step 3: Visualize Networks**
```python
from heartmap.utils import Visualizer

viz = Visualizer(config)

# Plot communication network
viz.plot_interaction_network(
    interactions,
    save_path='results/communication/figures/network.png'
)
```

### Tutorial 3: Chamber-Specific Analysis

**Goal:** Compare different heart chambers

**Step 1: Run Chamber Analysis**
```bash
heartmap your_data.h5ad --analysis-type multi-chamber --output-dir results/chambers/
```

**Step 2: Analyze Chamber Differences**
```python
# Load chamber-specific results
chamber_markers = {}
chambers = ['RA', 'RV', 'LA', 'LV']

for chamber in chambers:
    file_path = f'results/chambers/tables/markers_{chamber}.csv'
    if os.path.exists(file_path):
        chamber_markers[chamber] = pd.read_csv(file_path)

# Compare chamber marker genes
print("Top markers per chamber:")
for chamber, markers in chamber_markers.items():
    top_genes = markers.nlargest(5, 'score')['gene'].tolist()
    print(f"{chamber}: {', '.join(top_genes)}")
```

**Step 3: Cross-Chamber Analysis**
```python
# Load cross-chamber correlations
correlations = pd.read_csv('results/chambers/tables/cross_chamber_correlations.csv')

# Find most similar chambers
correlation_matrix = correlations.pivot(index='chamber1', columns='chamber2', values='correlation')
print("Chamber similarity matrix:")
print(correlation_matrix)
```

### Tutorial 4: Complete Analysis (Comprehensive)

**Goal:** Full HeartMAP analysis with all features

**Step 1: Configure Analysis**
```python
from heartmap import Config

config = Config.default()

# Optimize for your system
config.data.max_cells_subset = 50000  # Adjust based on RAM
config.data.max_genes_subset = 5000

# Customize analysis
config.analysis.resolution = 0.8
config.analysis.focus_chambers = ['LV', 'RV']  # Focus on ventricles

# Save configuration
config.to_yaml('my_config.yaml')
```

**Step 2: Run Comprehensive Analysis**
```bash
heartmap your_data.h5ad \
    --analysis-type comprehensive \
    --config my_config.yaml \
    --output-dir results/comprehensive/ \
    --verbose
```

**Step 3: Explore All Results**
```python
# Load comprehensive results
results_dir = 'results/comprehensive/'

# Basic results
basic_data = sc.read_h5ad(f'{results_dir}/basic/annotated_data.h5ad')

# Communication results
interactions = pd.read_csv(f'{results_dir}/communication/tables/interactions.csv')
hubs = pd.read_csv(f'{results_dir}/communication/tables/hubs.csv')

# Chamber results
chamber_stats = pd.read_csv(f'{results_dir}/chambers/tables/chamber_statistics.csv')

# View comprehensive report
import webbrowser
webbrowser.open(f'{results_dir}/comprehensive_report.html')
```

---

## 📊 Interpreting Results

### Understanding Output Files

```
results/
├── figures/                     # All visualizations
│   ├── umap_clusters.png       # Cell type clustering
│   ├── chamber_composition.png # Chamber cell distribution
│   ├── communication_hubs.png  # Communication hub cells
│   └── pathway_enrichment.png  # Enriched biological pathways
├── tables/                      # Analysis results as CSV files
│   ├── cell_annotations.csv    # Cell type assignments
│   ├── chamber_markers.csv     # Chamber-specific marker genes
│   ├── interactions.csv        # Cell-cell communication pairs
│   └── quality_metrics.csv     # Data quality statistics
├── processed_data/              # Processed datasets
│   └── heartmap_annotated.h5ad # Fully annotated data
└── reports/                     # Analysis summaries
    ├── summary_report.html      # Interactive HTML report
    └── methods_report.md        # Analysis methods used
```

### Key Metrics to Check

#### Data Quality
```python
# Load quality metrics
qc = pd.read_csv('results/tables/quality_metrics.csv')

print("Data Quality Summary:")
print(f"Cells analyzed: {qc['n_cells_final'].iloc[0]}")
print(f"Genes analyzed: {qc['n_genes_final'].iloc[0]}")
print(f"Median genes per cell: {qc['median_genes_per_cell'].iloc[0]}")
print(f"Median counts per cell: {qc['median_counts_per_cell'].iloc[0]}")
```

#### Cell Type Distribution
```python
# Check cell type balance
cell_types = pd.read_csv('results/tables/cell_annotations.csv')
type_counts = cell_types['cell_type'].value_counts()

print("Cell Type Distribution:")
for cell_type, count in type_counts.items():
    percentage = (count / len(cell_types)) * 100
    print(f"{cell_type}: {count} cells ({percentage:.1f}%)")
```

#### Chamber Analysis
```python
# Chamber composition
if 'chamber' in cell_types.columns:
    chamber_composition = pd.crosstab(
        cell_types['chamber'], 
        cell_types['cell_type'], 
        normalize='index'
    ) * 100
    
    print("Chamber Composition (% of cells):")
    print(chamber_composition.round(1))
```

### Biological Interpretation

#### Communication Hubs
```python
# Identify important communication hubs
hubs = pd.read_csv('results/tables/communication_hubs.csv')
top_hubs = hubs.nlargest(5, 'hub_score')

print("Key Communication Hubs:")
for _, hub in top_hubs.iterrows():
    print(f"- {hub['cell_type']}: Score {hub['hub_score']:.3f}")
    print(f"  Role: {hub['biological_role']}")
```

#### Pathway Enrichment
```python
# Check enriched pathways
if os.path.exists('results/tables/pathway_enrichment.csv'):
    pathways = pd.read_csv('results/tables/pathway_enrichment.csv')
    significant = pathways[pathways['padj'] < 0.05]
    
    print("Significantly Enriched Pathways:")
    for _, pathway in significant.nlargest(5, 'enrichment_score').iterrows():
        print(f"- {pathway['pathway_name']}: Score {pathway['enrichment_score']:.2f}")
```

---

## 🔧 Advanced Usage

### Memory Optimization

For large datasets or limited RAM:

```python
from heartmap import Config

config = Config.default()

# For 8GB RAM systems
config.data.max_cells_subset = 20000
config.data.max_genes_subset = 3000
config.data.process_in_chunks = True

# For 16GB+ RAM systems
config.data.max_cells_subset = 50000
config.data.max_genes_subset = 5000

# Enable memory-efficient processing
config.analysis.use_memory_efficient_methods = True
config.output.save_intermediate = False  # Don't save intermediate files
```

### Batch Processing

Analyze multiple datasets:

```python
import os
from pathlib import Path
from heartmap.pipelines import ComprehensivePipeline

config = Config.default()
pipeline = ComprehensivePipeline(config)

# Process multiple files
data_dir = Path('data/')
results_dir = Path('batch_results/')

for data_file in data_dir.glob('*.h5ad'):
    print(f"Processing {data_file.name}...")
    
    output_dir = results_dir / data_file.stem
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        results = pipeline.run(str(data_file), str(output_dir))
        print(f"✅ Completed {data_file.name}")
    except Exception as e:
        print(f"❌ Failed {data_file.name}: {e}")
```

### Custom Analysis Workflows

Create custom analysis pipelines:

```python
from heartmap.pipelines import BasePipeline
from heartmap.data import DataProcessor
from heartmap.utils import Visualizer

class CustomHeartPipeline(BasePipeline):
    def __init__(self, config):
        super().__init__(config)
        self.processor = DataProcessor(config)
        self.visualizer = Visualizer(config)
    
    def run(self, data_path, output_dir=None):
        # Custom preprocessing
        adata = self.processor.process_from_raw(data_path)
        
        # Custom analysis steps
        # ... your analysis code here ...
        
        # Custom visualizations
        if output_dir:
            self.visualizer.create_custom_plots(adata, output_dir)
        
        return {'adata': adata, 'custom_results': results}

# Use custom pipeline
custom_pipeline = CustomHeartPipeline(config)
results = custom_pipeline.run('data.h5ad', 'custom_results/')
```

### Integration with Other Tools

```python
# Export for CellChat
from heartmap.utils import ResultsExporter

exporter = ResultsExporter(config)
exporter.export_for_cellchat(results, 'cellchat_input/')

# Export for Seurat
exporter.export_for_seurat(results, 'seurat_input.h5ad')

# Export for scanpy
exporter.export_for_scanpy(results, 'scanpy_input.h5ad')
```

---

## 🆘 Troubleshooting

### Common Issues

#### Memory Errors
```
MemoryError: Unable to allocate array
```

**Solution:**
```python
# Reduce dataset size
config.data.max_cells_subset = 20000
config.data.max_genes_subset = 2000

# Enable chunked processing
config.data.process_in_chunks = True
config.data.chunk_size = 5000
```

#### Import Errors
```
ImportError: No module named 'heartmap'
```

**Solution:**
```bash
# Reinstall HeartMAP
pip uninstall heartmap
pip install heartmap

# Check installation
python -c "import heartmap; print(heartmap.__version__)"
```

#### Data Format Issues
```
ValueError: Data format not recognized
```

**Solution:**
```python
# Validate data format
from heartmap.data import DataValidator

validator = DataValidator()
is_valid, issues = validator.validate_h5ad('your_data.h5ad')

if not is_valid:
    print("Issues found:")
    for issue in issues:
        print(f"- {issue}")
```

#### Analysis Failures
```
AnalysisError: Pipeline failed at step X
```

**Solution:**
```python
# Enable debug mode
config.debug.verbose = True
config.debug.save_intermediate = True

# Use test mode for troubleshooting
config.data.test_mode = True
config.data.max_cells_subset = 1000
```

### Getting Help

1. **Check logs**: Look in `results/logs/` for detailed error messages
2. **Use test mode**: Set `config.data.test_mode = True` for debugging
3. **Validate data**: Run `heartmap --validate your_data.h5ad`
4. **Check documentation**: Visit [GitHub Wiki](https://github.com/Tumo505/HeartMap/wiki)
5. **Ask for help**: Post issues on [GitHub](https://github.com/Tumo505/HeartMap/issues)

### Performance Tips

```python
# Speed up analysis
config.analysis.skip_advanced_qc = True
config.visualization.create_plots = False  # Skip plots during development

# Optimize for your system
import psutil
available_memory = psutil.virtual_memory().available // (1024**3)  # GB

if available_memory >= 16:
    config.data.max_cells_subset = 50000
elif available_memory >= 8:
    config.data.max_cells_subset = 30000
else:
    config.data.max_cells_subset = 20000
```

---

## 🎓 Learning Resources

### Example Notebooks

```bash
# Download example notebooks
heartmap-examples --download

# Available notebooks:
# - basic_analysis.ipynb
# - communication_analysis.ipynb
# - chamber_comparison.ipynb
# - comprehensive_workflow.ipynb
```

---

## 🎯 Next Steps

1. **Start Simple**: Begin with basic analysis to understand your data
2. **Explore Features**: Try different analysis types based on your research questions
3. **Customize**: Adapt configurations for your specific needs
4. **Integrate**: Combine HeartMAP with other tools in your workflow
5. **Share**: Contribute your discoveries back to the community

**Happy analyzing! 🫀✨**

---

*For more help, visit the [HeartMAP GitHub repository](https://github.com/Tumo505/HeartMap) or email 28346416@mylife.unisa.ac.za*
