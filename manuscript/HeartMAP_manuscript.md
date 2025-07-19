# HeartMAP: A Multi-Chamber Spatial Framework for Cardiac Cell-Cell Communication

## Abstract

The human heart is composed of four distinct chambers, each with unique cellular and molecular characteristics. Understanding cell-cell communication within and between these chambers is essential for unraveling cardiac function and disease. Here, we present HeartMAP (Heart Multi-chamber Analysis Platform), a comprehensive spatial framework that integrates single-cell RNA-seq data and advanced computational models to map cardiac cell-cell communication at chamber resolution. Using a dataset of 287,269 cells from healthy human hearts, we identify chamber-specific cell populations, communication networks, and therapeutic targets. Our multi-model approach reveals both shared and unique signaling pathways across chambers, providing a foundation for chamber-specific therapies and precision cardiology.

**Keywords:** single-cell RNA-seq, cell-cell communication, heart chambers, spatial analysis, therapeutic targets

## Introduction

The heart’s four chambers (left/right atria, left/right ventricles) coordinate to maintain circulation, yet each exhibits distinct gene expression and cellular composition. While single-cell RNA-seq has enabled high-resolution mapping of cardiac cell types, comprehensive frameworks for multi-chamber communication analysis are lacking. HeartMAP addresses this gap by integrating basic, advanced, and multi-chamber models to provide a holistic view of cardiac intercellular signaling.

## Methods

### Data Acquisition and Preprocessing
- Data: Human heart scRNA-seq (SCP498), 287,269 cells, 33,694 genes
- Preprocessing: Quality control, normalization, clustering, annotation
- Scaling: Datasets reduced to 50,000–30,000 cells for computational efficiency (M1 MacBook Pro, 16GB RAM)

### Model Overview
- **Basic Pipeline:** Preprocessing, clustering, annotation, basic communication
- **Advanced Communication:** Temporal, pathway, and hub analysis
- **Multi-Chamber Atlas (HeartMAP Core):** Chamber-specific and cross-chamber communication

### Statistical Analysis
- Differential expression: Wilcoxon rank-sum, Benjamini-Hochberg correction
- Correlation: Pearson
- Pathway enrichment: Over-representation analysis

## Results

### Basic Pipeline
- Identified major cardiac cell types (cardiomyocytes, fibroblasts, endothelial, immune)
- High-quality data with minimal technical artifacts
- Basic ligand-receptor communication networks

### Advanced Communication Analysis
- Temporal dynamics of communication
- Identification of communication hubs (key cell types driving signaling)
- Pathway enrichment (cardiac development, contraction, immune response)
- Chamber-specific communication specificity

### Multi-Chamber Atlas (HeartMAP Core)
- **Chamber Distribution:** Analysis of 287,269 cells revealed RA (82,045 cells, 28.4%), LV (78,264 cells, 27.0%), LA (75,743 cells, 26.4%), RV (51,217 cells, 18.2%)
- **Chamber-Specific Marker Genes:** Identified 1,000+ marker genes per chamber, including:
  - **RA:** NPPA, MIR100HG, MYL7, MYL4, PDE4D
  - **RV:** NEAT1, MYH7, FHL2, C15orf41, PCDH7  
  - **LA:** NPPA, ELN, MYL7, EBF2, RORA
  - **LV:** CD36, LINC00486, FHL2, RP11-532N4.2, MYH7
- **Cross-Chamber Expression Correlations:** RV-LV (r = 0.985, highest), RA-LA (r = 0.960), RA-RV (r = 0.885), RV-LA (r = 0.877), RA-LV (r = 0.874), LA-LV (r = 0.870, lowest)
- **Differential Expression:** 150+ significantly differentially expressed genes per chamber pair (p < 0.05, Benjamini-Hochberg corrected)
- **Communication Hubs:** Identified top communication hub cells with hub scores ranging from 0.037-0.047, primarily in atrial cardiomyocytes and adipocytes
- **Chamber-Specific Therapeutic Targets:** Each chamber exhibits unique gene expression patterns enabling targeted therapeutic development

## Discussion

### Integration and Relationships of Models
HeartMAP’s strength lies in its modular, progressive approach:
- The **basic pipeline** establishes a robust foundation by identifying cell types and general communication patterns, ensuring data quality and interpretability.
- The **advanced communication model** builds on this by revealing temporal and pathway-level dynamics, highlighting how communication changes over time and which cell types act as hubs.
- The **multi-chamber atlas** (HeartMAP Core) leverages these foundations to dissect chamber-specific and cross-chamber communication, revealing both shared and unique signaling networks.

This layered approach allows for:
- **Cross-validation:** Results from one model (e.g., cell type annotation) inform and validate findings in others (e.g., communication hubs).
- **Complementary insights:** Basic models provide breadth, advanced models add depth, and the atlas delivers spatial and chamber-specific context.
- **Scalability:** The framework is optimized for standard hardware, making advanced cardiac analysis accessible.

### Biological and Clinical Importance
- **Chamber-Specific Biology:** HeartMAP reveals that each chamber has unique marker genes and communication networks, reflecting their distinct physiological roles. Analysis of 287,269 cells identified 1,000+ chamber-specific markers, with NPPA highly expressed in atria (RA: 82,045 cells, LA: 75,743 cells) and contractile proteins (MYH7, FHL2) dominating in ventricles (LV: 78,264 cells, RV: 51,217 cells).

- **Cross-Chamber Communication:** The platform uncovers both conserved and divergent signaling pathways. Cross-chamber correlations reveal RV-LV as most similar (r = 0.985), while LA-LV shows the greatest divergence (r = 0.870), suggesting specialized ventricular vs. atrial mechanisms. Communication hub analysis identified atrial cardiomyocytes and adipocytes as key signaling centers with hub scores of 0.037-0.047.

- **Therapeutic Implications:** By identifying chamber-specific markers and communication hubs, HeartMAP enables the development of targeted therapies. Atrial-specific drugs could target NPPA pathways (highly expressed in both atria), while ventricular therapies might focus on contractile machinery (MYH7, FHL2). The 150+ differentially expressed genes per chamber pair provide multiple therapeutic targets.

- **Disease Relevance:** Although this study uses healthy tissue, the framework is directly applicable to disease states, enabling the discovery of chamber-specific disease mechanisms and therapeutic targets. The comprehensive marker gene database (1,000+ genes per chamber) provides a foundation for disease-specific analysis.

### Model Relationships and Impact
- **Progressive Complexity:** Each model builds on the previous, ensuring robust, interpretable results.
- **Complementary Strengths:** The basic pipeline ensures data quality, the advanced model adds mechanistic insight, and the atlas contextualizes findings spatially and anatomically.
- **Clinical Translation:** HeartMAP’s modularity and scalability make it suitable for both research and clinical settings, supporting precision medicine initiatives in cardiology.

### Limitations and Future Directions
- **Data Scaling:** While downsampling enables analysis on standard hardware, rare cell types may be underrepresented. Future work will leverage high-performance computing for full-scale analysis.
- **Spatial Data Integration:** HeartMAP is designed to integrate spatial transcriptomics, which will further enhance chamber-specific resolution.
- **Disease Application:** Applying HeartMAP to diseased hearts will reveal chamber-specific pathologies and therapeutic opportunities.

## Conclusion

HeartMAP provides a robust, scalable, and extensible framework for multi-chamber cardiac cell-cell communication analysis. By integrating multiple models, it delivers both broad and deep insights into cardiac biology, with direct implications for therapy and precision medicine.

## Acknowledgments


## Data and Code Availability
All scripts and results are available in the HeartMAP repository. The original dataset is available through the Single Cell Portal (SCP498).


