# HeartMAP: Heart Multi-chamber Analysis Platform

## 🫀 Project Overview

HeartMAP is a comprehensive analysis platform for mapping cell-cell communication across all four chambers of the human heart using single-cell RNA sequencing data. The project integrates multiple models and analysis pipelines to provide insights into chamber-specific biology, cross-chamber signaling, and therapeutic targets.

## 📁 Project Structure

```
HeartMAP/
├── manuscript/                    # Manuscript and documentation
│   └── HeartMAP_manuscript.md     # Main manuscript
├── data/                          # Data files
│   ├── raw/                       # Raw data files
│   └── processed/                 # Processed data files
├── analysis/                      # Analysis scripts
│   ├── basic_pipeline/            # Basic single-cell analysis
│   ├── advanced_communication/    # Advanced communication analysis
│   └── multi_chamber_atlas/       # Multi-chamber atlas analysis
├── results/                       # Analysis results
└── figures/                       # Publication-ready figures
```

## 🔬 Analysis Components

### 1. Basic Pipeline Analysis
- Foundation single-cell analysis (preprocessing, QC, clustering, annotation, basic communication)

### 2. Advanced Communication Analysis
- Temporal, pathway, and hub analysis of cell-cell communication

### 3. Multi-Chamber Atlas (HeartMAP Core)
- Chamber-specific and cross-chamber communication analysis

## 📊 Key Findings

### Chamber Distribution
- **RA (Right Atrium):** 28.4% of cells
- **LV (Left Ventricle):** 27.0% of cells
- **LA (Left Atrium):** 26.4% of cells
- **RV (Right Ventricle):** 18.2% of cells

### Chamber-Specific Markers
- **RA:** NPPA, MIR100HG, MYL7, MYL4, PDE4D
- **RV:** NEAT1, MYH7, FHL2, C15orf41, PCDH7
- **LA:** NPPA, ELN, MYL7, EBF2, RORA
- **LV:** CD36, LINC00486, FHL2, RP11-532N4.2, MYH7

### Cross-Chamber Correlations
- **RV vs LV:** r = 0.985 (highest correlation)
- **RA vs LA:** r = 0.960
- **LA vs LV:** r = 0.870 (lowest correlation)

## 🎯 Clinical Implications

- Chamber-specific treatments may improve outcomes
- Chamber-specific drug targets may be more effective
- Chamber-specific analysis may reveal disease mechanisms

## 🚀 Next Steps

1. Validate markers and communication patterns with literature and experiments
2. Integrate spatial transcriptomics and disease data
3. Develop chamber-specific clinical applications

