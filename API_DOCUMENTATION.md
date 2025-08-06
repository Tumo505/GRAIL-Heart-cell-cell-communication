# HeartMAP API Documentation

## 📚 Table of Contents

- [Python API Reference](#python-api-reference)
- [Command Line Interface](#command-line-interface)
- [REST API](#rest-api)
- [Web Interface](#web-interface)
- [Configuration API](#configuration-api)
- [Data Processing API](#data-processing-api)
- [Visualization API](#visualization-api)

---

## 🐍 Python API Reference

### Core Classes

#### `heartmap.Config`

Configuration management for HeartMAP analysis.

```python
from heartmap import Config

# Create default configuration
config = Config.default()

# Load from file
config = Config.from_yaml('config.yaml')
config = Config.from_json('config.json')

# Customize configuration
config.data.max_cells_subset = 50000
config.analysis.resolution = 0.8
config.output.save_figures = True

# Save configuration
config.to_yaml('my_config.yaml')
config.to_json('my_config.json')
```

**Key Properties:**
- `config.data.*` - Data processing parameters
- `config.analysis.*` - Analysis parameters
- `config.output.*` - Output configuration
- `config.visualization.*` - Plotting parameters

#### `heartmap.pipelines.BasicPipeline`

Foundation single-cell heart analysis pipeline.

```python
from heartmap.pipelines import BasicPipeline

# Initialize pipeline
pipeline = BasicPipeline(config)

# Run analysis
results = pipeline.run(
    data_path='heart_data.h5ad',
    output_dir='results/basic/'
)

# Access results
adata = results['adata']  # Annotated data
clustering = results['clustering']  # Cell clusters
qc_metrics = results['qc_metrics']  # Quality control
```

**Methods:**
- `run(data_path, output_dir=None)` - Execute analysis pipeline
- `save_results(output_dir)` - Save analysis results
- `load_results(output_dir)` - Load previous results

#### `heartmap.pipelines.AdvancedCommunicationPipeline`

Advanced cell-cell communication analysis.

```python
from heartmap.pipelines import AdvancedCommunicationPipeline

pipeline = AdvancedCommunicationPipeline(config)
results = pipeline.run('heart_data.h5ad', 'results/communication/')

# Access communication results
interactions = results['interactions']  # Cell-cell interactions
hubs = results['communication_hubs']  # Communication hub cells
pathways = results['pathway_enrichment']  # Enriched pathways
networks = results['networks']  # Communication networks
```

#### `heartmap.pipelines.MultiChamberPipeline`

Chamber-specific analysis for cardiac data.

```python
from heartmap.pipelines import MultiChamberPipeline

pipeline = MultiChamberPipeline(config)
results = pipeline.run('heart_data.h5ad', 'results/chambers/')

# Access chamber-specific results
chamber_markers = results['chamber_markers']  # Markers per chamber
correlations = results['cross_chamber_correlations']  # Chamber relationships
chamber_stats = results['chamber_statistics']  # Chamber cell counts
```

#### `heartmap.pipelines.ComprehensivePipeline`

Complete HeartMAP analysis combining all pipelines.

```python
from heartmap.pipelines import ComprehensivePipeline

pipeline = ComprehensivePipeline(config)
results = pipeline.run('heart_data.h5ad', 'results/comprehensive/')

# Access all results
basic_results = results['basic']
communication_results = results['communication']
chamber_results = results['multi_chamber']
```

### Data Processing

#### `heartmap.data.DataProcessor`

Data preprocessing and quality control.

```python
from heartmap.data import DataProcessor

processor = DataProcessor(config)

# Process raw data
adata = processor.process_from_raw('raw_data.h5ad')

# Quality control
adata_filtered = processor.quality_control(adata)

# Normalization
adata_norm = processor.normalize(adata_filtered)

# Feature selection
adata_features = processor.select_features(adata_norm)
```

#### `heartmap.data.DataLoader`

Data loading utilities.

```python
from heartmap.data import DataLoader

loader = DataLoader(config)

# Load various formats
adata = loader.load_h5ad('data.h5ad')
adata = loader.load_csv('data.csv', obs_file='metadata.csv')
adata = loader.load_mtx('matrix.mtx', features='features.tsv', barcodes='barcodes.tsv')

# Validate data
is_valid, issues = loader.validate_format(adata)
```

### Visualization

#### `heartmap.utils.Visualizer`

Scientific visualization for cardiac analysis.

```python
from heartmap.utils import Visualizer

viz = Visualizer(config)

# Basic plots
viz.plot_umap(adata, color='cell_type', save_path='umap.png')
viz.plot_quality_metrics(adata, save_path='qc.png')

# Communication plots
viz.plot_communication_hubs(adata, hub_scores, save_path='hubs.png')
viz.plot_interaction_network(interactions, save_path='network.png')

# Chamber-specific plots
viz.plot_chamber_heatmap(chamber_markers, save_path='chambers.png')
viz.plot_cross_chamber_correlations(correlations, save_path='correlations.png')

# Comprehensive dashboard
viz.create_comprehensive_dashboard(adata, results, save_dir='figures/')
```

### Utilities

#### `heartmap.utils.ResultsExporter`

Export analysis results in various formats.

```python
from heartmap.utils import ResultsExporter

exporter = ResultsExporter(config)

# Export tables
exporter.export_csv(results['markers'], 'markers.csv')
exporter.export_excel(results, 'analysis_results.xlsx')

# Export figures
exporter.export_publication_figures(results, 'publication_figures/')

# Generate reports
exporter.generate_html_report(results, 'report.html')
exporter.generate_pdf_report(results, 'report.pdf')
```

---

## 💻 Command Line Interface

### Basic Usage

```bash
# Analyze data with defaults
heartmap data.h5ad

# Specify analysis type
heartmap data.h5ad --analysis-type basic
heartmap data.h5ad --analysis-type communication
heartmap data.h5ad --analysis-type multi-chamber
heartmap data.h5ad --analysis-type comprehensive

# Custom output directory
heartmap data.h5ad --output-dir results/my_analysis/

# Use custom configuration
heartmap data.h5ad --config my_config.yaml

# Memory optimization
heartmap data.h5ad --max-cells 30000 --max-genes 3000
```

### Advanced Options

```bash
# List available analyses
heartmap --list-analyses

# Validate data format
heartmap --validate data.h5ad

# Check configuration
heartmap --check-config config.yaml

# Version information
heartmap --version

# Verbose output
heartmap data.h5ad --verbose

# Dry run (check parameters without running)
heartmap data.h5ad --dry-run
```

### Examples

```bash
# Research pipeline
heartmap patient_data.h5ad \
    --analysis-type comprehensive \
    --output-dir results/patient_study/ \
    --config research_config.yaml \
    --verbose

# High-throughput analysis
heartmap large_dataset.h5ad \
    --analysis-type basic \
    --max-cells 50000 \
    --output-dir results/htp_analysis/ \
    --threads 8

# Development/testing
heartmap test_data.h5ad \
    --analysis-type basic \
    --max-cells 5000 \
    --test-mode \
    --output-dir test_results/
```

---

## 🌐 REST API

### Starting the API Server

```bash
# Start with defaults
python -m heartmap.api.rest

# Custom port and host
python -m heartmap.api.rest --host 0.0.0.0 --port 8080
```

### Endpoints

#### `POST /analyze`

Analyze single-cell data.

**Request:**
```bash
curl -X POST "http://localhost:8000/analyze" \
     -F "file=@data.h5ad" \
     -F "analysis_type=comprehensive" \
     -F "output_format=json"
```

**Response:**
```json
{
    "status": "success",
    "analysis_id": "abc123",
    "results": {
        "summary": {
            "n_cells": 45000,
            "n_genes": 3500,
            "analysis_completed": true
        },
        "annotation_summary": {...},
        "communication_summary": {...},
        "chamber_summary": {...}
    },
    "download_urls": {
        "figures": "/download/abc123/figures.zip",
        "tables": "/download/abc123/tables.zip",
        "report": "/download/abc123/report.html"
    }
}
```

#### `GET /models`

List available analysis pipelines.

**Response:**
```json
{
    "available_models": [
        "basic",
        "advanced_communication", 
        "multi_chamber",
        "comprehensive"
    ],
    "descriptions": {
        "basic": "Foundation single-cell analysis",
        "advanced_communication": "Cell-cell communication analysis",
        "multi_chamber": "Chamber-specific analysis",
        "comprehensive": "Complete HeartMAP analysis"
    }
}
```

#### `POST /config`

Update analysis configuration.

**Request:**
```json
{
    "data": {
        "max_cells_subset": 40000,
        "max_genes_subset": 4000
    },
    "analysis": {
        "resolution": 0.8
    }
}
```

#### `GET /status/{analysis_id}`

Check analysis status.

**Response:**
```json
{
    "analysis_id": "abc123",
    "status": "running",
    "progress": 65,
    "current_step": "communication_analysis",
    "estimated_time_remaining": "5 minutes"
}
```

---

## 🖥️ Web Interface

### Starting the Web Interface

```python
from heartmap.api import run_web_interface

# Start with defaults
run_web_interface()

# Custom configuration
run_web_interface(
    port=7860,
    share=True,  # Create public link
    debug=False
)
```

### Features

- **File Upload**: Drag-and-drop .h5ad files
- **Analysis Selection**: Choose from available pipelines
- **Parameter Configuration**: Adjust analysis parameters
- **Real-time Progress**: Monitor analysis progress
- **Results Download**: Download figures, tables, and reports
- **Interactive Plots**: Explore results interactively

---

## ⚙️ Configuration API

### Configuration Structure

```python
from heartmap import Config

config = Config.default()

# Data configuration
config.data.max_cells_subset = 50000      # Maximum cells to analyze
config.data.max_genes_subset = 5000       # Maximum genes to analyze
config.data.test_mode = False             # Test mode flag
config.data.random_seed = 42              # Reproducibility seed

# Analysis configuration
config.analysis.resolution = 0.8          # Clustering resolution
config.analysis.n_neighbors = 15          # UMAP neighbors
config.analysis.min_genes = 200           # Minimum genes per cell
config.analysis.focus_chambers = None     # Specific chambers to analyze

# Output configuration
config.output.save_figures = True         # Save visualization figures
config.output.figure_format = 'png'       # Figure format
config.output.save_intermediate = True    # Save intermediate results
config.output.generate_report = True      # Generate analysis report

# Visualization configuration
config.visualization.figure_size = (10, 8)  # Default figure size
config.visualization.dpi = 300               # Figure resolution
config.visualization.color_palette = 'tab10' # Color palette
```

### Dynamic Configuration

```python
# Update configuration at runtime
config.update({
    'data': {'max_cells_subset': 30000},
    'analysis': {'resolution': 1.0},
    'output': {'save_figures': False}
})

# Conditional configuration
if memory_limited:
    config.data.max_cells_subset = 20000
    config.data.process_in_chunks = True

# Environment-specific configuration
if config.environment == 'production':
    config.output.save_intermediate = False
    config.visualization.create_plots = False
```

---

## 🔧 Data Processing API

### Quality Control

```python
from heartmap.data import DataProcessor

processor = DataProcessor(config)

# Cell filtering
adata_cells = processor.filter_cells(
    adata,
    min_genes=200,      # Minimum genes per cell
    max_genes=5000,     # Maximum genes per cell
    min_counts=1000,    # Minimum counts per cell
    max_mito_percent=20 # Maximum mitochondrial percentage
)

# Gene filtering
adata_genes = processor.filter_genes(
    adata_cells,
    min_cells=10,       # Minimum cells expressing gene
    min_counts=50       # Minimum total counts
)

# Quality metrics
qc_metrics = processor.calculate_qc_metrics(adata)
```

### Normalization and Scaling

```python
# Normalization
adata_norm = processor.normalize(
    adata,
    target_sum=1e4,           # Target counts per cell
    log_transform=True,       # Log transformation
    highly_variable_genes=True # Select highly variable genes
)

# Scaling
adata_scaled = processor.scale(
    adata_norm,
    max_value=10,     # Maximum scaled value
    zero_center=True  # Center on zero
)

# Principal component analysis
adata_pca = processor.run_pca(
    adata_scaled,
    n_comps=50,       # Number of components
    svd_solver='arpack' # Solver method
)
```

### Feature Selection

```python
# Highly variable genes
processor.find_highly_variable_genes(
    adata,
    n_top_genes=3000,     # Number of genes to select
    flavor='seurat_v3',   # Selection method
    subset=True           # Keep only selected genes
)

# Chamber-specific features
chamber_features = processor.find_chamber_specific_genes(
    adata,
    chamber_key='chamber',
    min_fold_change=2.0,
    max_pvalue=0.05
)
```

---

## 📊 Visualization API

### Basic Plots

```python
from heartmap.utils import Visualizer

viz = Visualizer(config)

# UMAP plot
viz.plot_umap(
    adata,
    color='cell_type',        # Coloring variable
    palette='tab10',          # Color palette
    size=50,                  # Point size
    alpha=0.7,               # Transparency
    save_path='umap.png'     # Save location
)

# Quality control plots
viz.plot_qc_metrics(
    adata,
    metrics=['n_genes', 'n_counts', 'mito_percent'],
    save_path='qc_metrics.png'
)

# Gene expression plot
viz.plot_gene_expression(
    adata,
    genes=['NPPA', 'MYH7', 'CD36'],
    save_path='gene_expression.png'
)
```

### Communication Plots

```python
# Communication hubs
viz.plot_communication_hubs(
    adata,
    hub_scores,
    threshold=0.7,
    save_path='communication_hubs.png'
)

# Interaction network
viz.plot_interaction_network(
    interactions,
    layout='spring',          # Network layout
    node_size_col='hub_score', # Node size variable
    edge_width_col='strength', # Edge width variable
    save_path='network.png'
)

# Pathway enrichment
viz.plot_pathway_enrichment(
    enrichment_results,
    top_n=20,
    save_path='pathways.png'
)
```

### Chamber-Specific Plots

```python
# Chamber composition
viz.plot_chamber_composition(
    adata,
    chamber_key='chamber',
    cell_type_key='cell_type',
    save_path='chamber_composition.png'
)

# Chamber marker heatmap
viz.plot_chamber_heatmap(
    chamber_markers,
    top_n_genes=10,
    save_path='chamber_markers.png'
)

# Cross-chamber correlations
viz.plot_cross_chamber_correlations(
    correlations,
    method='spearman',
    save_path='correlations.png'
)
```

### Advanced Visualizations

```python
# Multi-panel figure
viz.create_multi_panel_figure(
    adata,
    panels=[
        {'type': 'umap', 'color': 'cell_type'},
        {'type': 'umap', 'color': 'chamber'},
        {'type': 'heatmap', 'genes': marker_genes},
        {'type': 'violin', 'genes': ['NPPA', 'MYH7']}
    ],
    save_path='multi_panel.png'
)

# Interactive plots
viz.create_interactive_plot(
    adata,
    plot_type='umap',
    color='cell_type',
    save_path='interactive_umap.html'
)

# Publication-ready figures
viz.create_publication_figure(
    adata,
    results,
    figure_type='main',
    save_path='publication_main.png'
)
```

---

## 🔍 Example Workflows

### Complete Analysis Workflow

```python
from heartmap import Config
from heartmap.pipelines import ComprehensivePipeline
from heartmap.utils import Visualizer, ResultsExporter

# 1. Setup
config = Config.default()
config.data.max_cells_subset = 50000

# 2. Analysis
pipeline = ComprehensivePipeline(config)
results = pipeline.run('heart_data.h5ad', 'results/')

# 3. Visualization
viz = Visualizer(config)
viz.create_comprehensive_dashboard(
    results['adata'], 
    results, 
    'results/figures/'
)

# 4. Export
exporter = ResultsExporter(config)
exporter.generate_html_report(results, 'results/report.html')
exporter.export_publication_figures(results, 'results/publication/')
```

### Chamber-Specific Analysis

```python
from heartmap.pipelines import MultiChamberPipeline

# Focus on specific chambers
config.analysis.focus_chambers = ['LV', 'RV']  # Ventricles only

pipeline = MultiChamberPipeline(config)
results = pipeline.run('heart_data.h5ad', 'results/ventricles/')

# Access ventricle-specific results
lv_markers = results['chamber_markers']['LV']
rv_markers = results['chamber_markers']['RV']
lv_rv_correlation = results['cross_chamber_correlations']['LV_vs_RV']
```

### Communication Analysis

```python
from heartmap.pipelines import AdvancedCommunicationPipeline

pipeline = AdvancedCommunicationPipeline(config)
results = pipeline.run('heart_data.h5ad', 'results/communication/')

# Identify key communication hubs
hub_cells = results['communication_hubs']
top_hubs = hub_cells.nlargest(10, 'hub_score')

# Analyze specific pathways
cardiac_pathways = results['pathway_enrichment'][
    results['pathway_enrichment']['pathway'].str.contains('cardiac')
]
```

---

## 🐛 Error Handling

### Common Exceptions

```python
from heartmap.exceptions import (
    HeartMAPError,
    DataFormatError,
    ConfigurationError,
    AnalysisError
)

try:
    pipeline = ComprehensivePipeline(config)
    results = pipeline.run('data.h5ad')
except DataFormatError as e:
    print(f"Data format issue: {e}")
except ConfigurationError as e:
    print(f"Configuration problem: {e}")
except AnalysisError as e:
    print(f"Analysis failed: {e}")
except HeartMAPError as e:
    print(f"HeartMAP error: {e}")
```

### Validation

```python
from heartmap.data import DataValidator

# Validate input data
validator = DataValidator()
is_valid, issues = validator.validate_h5ad('data.h5ad')

if not is_valid:
    print("Data validation issues:")
    for issue in issues:
        print(f"- {issue}")
else:
    print("Data is valid for HeartMAP analysis")
```

---

This API documentation provides comprehensive coverage of all HeartMAP functionality. For more examples and tutorials, visit the [HeartMAP GitHub repository](https://github.com/Tumo505/HeartMap).
