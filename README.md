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

## 🔒 Data Integrity

### SHA-256 Checksums
To ensure the integrity of raw data files, the HeartMAP project uses **SHA-256 checksums**. This guarantees that the raw data files have not been modified or corrupted during storage or transfer.

### Where SHA-256 is Used:
1. **Raw Data Verification**:
   - Before preprocessing, the pipeline verifies the integrity of raw data files in the `data/raw/` directory using a `checksums.txt` file.
   - Example:
     ```bash
     python utils/sha256_checksum.py verify data/raw data/raw/checksums.txt
     ```

2. **Checksum Generation**:
   - SHA-256 checksums are generated for all raw data files and stored in `checksums.txt`.
   - Example:
     ```bash
     python utils/sha256_checksum.py generate data/raw data/raw/checksums.txt
     ```

### Why SHA-256?
Using SHA-256 ensures:
- Data integrity during storage and transfer.
- Detection of accidental or malicious modifications to raw data files.
- Reproducibility by verifying that the same raw data is used across different runs.

## 🔄 Reproducibility

The HeartMAP project ensures reproducibility in all stochastic processes by using fixed random seeds. This guarantees consistent results across different runs of the pipeline. Below are the key areas where random seeds are applied:

1. **Random Sampling**:
   - Fixed seed (`seed = 42`) is used for randomly sampling cells during data scaling and preprocessing.
   - Example:
     ```python
     np.random.seed(42)
     cell_indices = np.random.choice(adata.n_obs, size=50000, replace=False)
     ```

2. **Mock Data Generation**:
   - Fixed seed is used to generate mock communication interactions for testing and visualization.
   - Example:
     ```python
     np.random.seed(42)
     mock_interactions = pd.DataFrame({
         'source': np.random.choice(cell_types.index, n_interactions),
         'target': np.random.choice(cell_types.index, n_interactions),
         'score': np.random.uniform(0, 1, n_interactions)
     })
     ```

3. **Clustering**:
   - Fixed seed (`random_state=42`) is used for K-means clustering as a fallback when Leiden clustering is unavailable.
   - Example:
     ```python
     kmeans = KMeans(n_clusters=n_clusters, random_state=42)
     ```

4. **LIANA Analysis**:
   - Fixed seed is passed to the LIANA analysis function for reproducibility in ligand-receptor interaction analysis.
   - Example:
     ```python
     li.mt.rank_aggregate.by_sample(
         adata,
         groupby=cell_type_col,
         resource_name='consensus',
         n_perms=100,
         seed=42,
         verbose=True
     )
     ```

### Why Fixed Seeds?
Using fixed random seeds ensures:
- Consistent results across different runs of the pipeline.
- Easier debugging and validation of results.
- Reproducibility for scientific publications and collaborative workflows.

## 🚀 Next Steps

1. Validate markers and communication patterns with literature and experiments
2. Integrate spatial transcriptomics and disease data
3. Develop chamber-specific clinical applications

