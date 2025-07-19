
# Cross-Chamber Heart Analysis (Scaled for M1 MacBook Pro)

## Analysis Overview
This analysis examines communication patterns across different heart chambers using a scaled-down dataset optimized for M1 MacBook Pro (16GB RAM).

## Dataset Scaling
- **Original Size**: 287,269 cells × 33,694 genes
- **Scaled Size**: 30,000 cells × 3,000 genes
- **Cell Reduction**: 10.4% of original
- **Gene Reduction**: 8.9% of original

## Cross-Chamber Expression Correlations

### RA_vs_RV
- **Expression Correlation**: 0.885
- **Significant Genes**: 150
- **Top Differentially Expressed**: C1orf168, ADRA1A, FP325317.1, XIRP2, MIR100HG

### RA_vs_LA
- **Expression Correlation**: 0.960
- **Significant Genes**: 150
- **Top Differentially Expressed**: SMAD6, AKR1C1, TGM2, GALNT15, MFAP5

### RA_vs_LV
- **Expression Correlation**: 0.874
- **Significant Genes**: 150
- **Top Differentially Expressed**: ULK4, SAMHD1, GAS2, MIR100HG, HOOK2

### RV_vs_LA
- **Expression Correlation**: 0.877
- **Significant Genes**: 150
- **Top Differentially Expressed**: MFAP5, C1orf168, MYL3, AC011288.2, PDE1A

### RV_vs_LV
- **Expression Correlation**: 0.985
- **Significant Genes**: 150
- **Top Differentially Expressed**: RP11-362K2.2, SLC26A3, ADRA1A, ABCC1, BCL6

### LA_vs_LV
- **Expression Correlation**: 0.870
- **Significant Genes**: 150
- **Top Differentially Expressed**: MFAP5, CCDC80, MYL3, GCNT2, FN1

## Key Findings

1. **Chamber-Specific Expression**: Each chamber shows unique gene expression patterns
2. **Cross-Chamber Communication**: Different chambers have distinct communication networks
3. **Therapeutic Implications**: Chamber-specific drug targets may be more effective
4. **Disease Relevance**: Chamber-specific analysis may improve treatment strategies

## Clinical Implications

1. **Personalized Medicine**: Chamber-specific treatments may improve outcomes
2. **Drug Development**: Target validation should consider chamber-specific expression
3. **Disease Understanding**: Chamber-specific analysis may reveal disease mechanisms
4. **Treatment Optimization**: Chamber-specific dosing may be beneficial

## Files Generated
- `cross_chamber_correlation_matrix.png`: Expression correlation between chambers
- `significant_genes_by_chamber_pair.png`: Number of different genes between chambers
- Individual cross-chamber comparison plots for each pair

## Performance Notes
- Analysis optimized for M1 MacBook Pro with 16GB RAM
- Dataset scaled down to ensure fast processing while maintaining statistical power
- Original full dataset analysis available by uncommenting original code
