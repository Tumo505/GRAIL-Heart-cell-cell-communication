import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def find_variable_genes(adata, n_top_genes=2000):
    """Find highly variable genes"""
    print("Finding highly variable genes...")
    
    # This computes statistics (mean, variance, dispersion) for each gene and identifies the set of highly variable genes, which are most informative for downstream analysis.
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=n_top_genes)
    
    # Plot highly variable genes
    sc.pl.highly_variable_genes(adata)
    plt.savefig(Path("results/figures/highly_variable_genes.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    return adata

def perform_pca(adata, n_comps=50):
    """Perform Principal Component Analysis(PCA)"""
    print("Performing PCA...")
    
    # Keep only highly variable genes for PCA
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    
    # Scale data
    sc.pp.scale(adata, max_value=10)
    
    # PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)
    
    # Plot PCA
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
    plt.savefig(Path("results/figures/pca_variance.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    return adata

def compute_neighborhood_graph(adata, n_neighbors=10, n_pcs=40):
    """Compute neighborhood graph"""
    print("Computing neighborhood graph...")
    
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    
    return adata

def perform_clustering(adata, resolution=0.5):
    """Perform Leiden clustering"""
    print("Performing clustering...")
    
    sc.tl.leiden(adata, resolution=resolution)
    
    return adata

def compute_umap(adata):
    """Compute UMAP embedding"""
    print("Computing UMAP...")
    
    sc.tl.umap(adata)
    
    return adata

def plot_clustering_results(adata, save_path):
    """Plot clustering results"""
    print("Plotting clustering results...")
    
    # Plot UMAP with clusters
    sc.pl.umap(adata, color=['leiden'], legend_loc='on data',
               title='Leiden Clustering', frameon=False, save='_leiden_clustering.png')
    
    # If cell type annotations exist, plot them too
    if 'cell_type' in adata.obs.columns:
        sc.pl.umap(adata, color=['cell_type'], legend_loc='on data',
                   title='Cell Types', frameon=False, save='_cell_types.png')

def identify_marker_genes(adata, groupby='leiden', n_genes=25):
    """Identify marker genes for each cluster"""
    print("Identifying marker genes...")
    
    sc.tl.rank_genes_groups(adata, groupby, method='wilcoxon')
    
    # Plot marker genes
    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False)
    plt.savefig(Path("results/figures/marker_genes.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save marker genes to file
    marker_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    marker_genes.to_csv(Path("results/tables/marker_genes.csv"))
    
    return adata

def main():
    # Load filtered and normalized data
    data_path = Path("data/processed/heart_data_filtered_normalized.h5ad")
    adata = sc.read_h5ad(data_path)
    
    # Create results directories
    Path("results/figures").mkdir(parents=True, exist_ok=True)
    Path("results/tables").mkdir(parents=True, exist_ok=True)
    
    # Find highly variable genes using a subset for speed (temporary)
    adata_subset = adata[:10000, :1000]  # Subset to 10,000 cells, 1,000 genes (the full dataset is ~260,000 cells and ~28,000 genes)
    adata_subset = find_variable_genes(adata_subset)

    # Use the subset for downstream steps (temporary)
    adata = adata_subset
    
    # Perform PCA
    adata = perform_pca(adata)
    
    # Compute neighborhood graph
    adata = compute_neighborhood_graph(adata)
    
    # Perform clustering
    adata = perform_clustering(adata)
    
    # Compute UMAP
    adata = compute_umap(adata)
    
    # Plot results
    plot_clustering_results(adata, Path("results/figures"))
    
    # Identify marker genes
    adata = identify_marker_genes(adata)
    
    # Save annotated data
    adata.write(Path("data/processed/heart_data_annotated.h5ad"))
    
    print("Cell annotation complete!")

if __name__ == "__main__":
    main()