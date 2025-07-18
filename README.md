[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)

# GRAIL-Heart: Graph-based Reconstruction of Artificial Intercellular Links

This project, "**cell_comm**", is a component of the broader `GRAIL-Heart (Graph-based Reconstruction of Artificial Intercellular Links)` project, developed under the [UNISA Biomedical Engineering Research Group](https://www.unisa.ac.za/sites/corporate/default/Colleges/Science,-Engineering-&-Technology/Schools,-departments-&-institutes/School-of-Engineering-and-Built-Environment/Department-of-Mechanical-Bioresources-and-Biomedical-Engineering), under the `Department of Mechanical, Bioresources and Biomedical Engineering` at the [University of South Africa (UNISA)](https://www.unisa.ac.za).

GRAIL-Heart is a novel computational framework that leverages deep learning—specifically, a hybrid of Graph Neural Networks (GNNs) and Recurrent Neural Networks (RNNs)—to infer functional signalling networks and predict cardiomyocyte differentiation efficiency from spatial multi-omics data. The project integrates spatial transcriptomics, proteomics, and epigenomics data from developing iPSC-derived cardiomyocytes (iPSC-CMs) to identify critical molecular markers and spatial correlations that influence differentiation outcomes. By combining spatial (GNN) and temporal (RNN) modeling, GRAIL-Heart provides a comprehensive approach to understanding and optimizing cardiac tissue regeneration.

## About This Project (`cell_comm`)

This particular project, **cell_comm**, focuses on analyzing cell-cell communication patterns in human heart tissue using single-cell RNA sequencing data. It serves as a foundational step for the GRAIL-Heart framework by mapping intercellular signalling networks and tissue architecture at single-cell resolution. The insights gained here inform the design and training of GNN models for simulating intercellular interactions and RNN models for capturing temporal gene expression dynamics.

### Key Features

- **Cell-Cell Communication Analysis:**  
  This project identifies and characterizes intercellular signalling networks in cardiac tissue using advanced computational tools and single-cell transcriptomics data.

- **Spatial Multi-Omics Integration:**  
  The analysis is designed to be extensible to spatial transcriptomics, proteomics, and epigenomics, supporting the broader GRAIL-Heart goal of multi-modal data integration.

- **Deep Learning Foundations:**  
  Results from this project directly inform the construction of GNNs for spatial modeling and RNNs for temporal modeling, ultimately supporting a hybrid GNN-RNN framework for predicting iPSC-CM differentiation efficiency.

- **Regenerative Medicine Impact:**  
  By uncovering the molecular and spatial determinants of cardiomyocyte differentiation, this work advances the precision of stem cell-based cardiac regeneration therapies.

## Data Source

For this project, we use data from the following study:

> **He, S., Wang, L.H., Liu, Y. et al. Single-cell transcriptome profiling of an adult human cell atlas of 15 major organs. Genome Biol 21, 294 (2020). [https://doi.org/10.1186/s13059-020-02210-0](https://doi.org/10.1186/s13059-020-02210-0)**  
> [Read the article PDF](https://rdcu.be/ewSLg)  
> [Springer Article Link](https://link.springer.com/article/10.1186/s13059-020-02210-0#availability-of-data-and-materials)

- **Dataset:** `healthy_human_4chamber_map_unnormalized_V3.h5ad`
- **Source:** [Single Cell Portal (SCP498)](https://singlecell.broadinstitute.org/single_cell/study/SCP498)
- **Description:** This dataset provides comprehensive single-cell RNA-seq profiles capturing the transcriptional and cellular diversity of the human heart, as described in the above publication.

## Project Structure

- `data/` – Raw and processed data
- `results/` – Output figures and tables
- `scripts/` – Analysis scripts (quality control, annotation, communication analysis)
- `README.md` – Project overview and instructions
- `requirements.txt` – Python dependencies

## Installation and Reproduction Steps

To reproduce the analysis in this project:

1. **Clone or Download the Repository**
   ```bash
   git clone https://github.com/Tumo505/GRAIL-Heart-cell-cell-communication-.git
   cd grail-heart
   ```

2. **Create and Activate a Conda Environment**
   ```bash
   conda create -n cell_comm python=3.9
   conda activate cell_comm
   ```

3. **Install Required Dependencies**
   ```bash
   pip install -r requirements.txt
   ```

4. **Download the Dataset**
   - Register and download `healthy_human_4chamber_map_unnormalized_V3.h5ad` from the [Single Cell Portal (SCP498)](https://singlecell.broadinstitute.org/single_cell/study/SCP498).
   - Place the file in the `data/raw/` directory.

5. **Run the Analysis Pipeline**
   - **Quality Control:**  
     ```bash
     python scripts/02_quality_control.py
     ```
   - **Cell Annotation:**  
     ```bash
     python scripts/03_cell_annotation.py
     ```
   - **Cell-Cell Communication Analysis:**  
     ```bash
     python scripts/04_communication_analysis.py
     ```

6. **Results**
   - Processed data, figures, and tables are saved in the `results/` and `data/processed/` directories.

## Acknowledgements

This project is developed as part of the GRAIL-Heart project at the `UNISA Biomedical Engineering Research Group` under the `Department of Mechanical, Bioresources and Biomedical Engineering`. It leverages open-source tools and publicly available datasets to advance the field of cardiac regenerative medicine.

For questions or collaboration, please contact Tumo Kgabeng: [28346416@mylife.unisa.ac.za](mailto:28346416@mylife.unisa.ac.za).

## License

This project is licensed under the Apache License 2.0. See the [LICENSE](LICENSE) file for details.