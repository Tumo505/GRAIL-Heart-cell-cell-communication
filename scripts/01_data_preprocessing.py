import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set up scanpy
sc.settings.verbosity = 3
sc.set_figure_params(dpi=80, facecolor='white')

def load_and_preprocess_data(data_path):
    """Load and preprocess the heart dataset"""
    print("Loading data...")
    adata = sc.read_h5ad(data_path)
    
    # Basic info
    print(f"Dataset shape: {adata.shape}")
    print(f"Number of cells: {adata.n_obs}")
    print(f"Number of genes: {adata.n_vars}")
    
    # Make gene names unique
    adata.var_names_make_unique()
    
    # Store raw data
    adata.raw = adata
    
    return adata

def basic_filtering(adata, min_genes=200, min_cells=3):
    """Basic filtering of cells and genes"""
    print("Performing basic filtering...")
    
    # Filter cells with too few genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    # Filter genes present in too few cells
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    print(f"After filtering: {adata.shape}")
    return adata

def calculate_qc_metrics(adata):
    """Calculate quality control metrics"""
    print("Calculating QC metrics...")
    
    # Mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    
    # Ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    
    # Hemoglobin genes
    adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]')
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], 
                              percent_top=None, log1p=False, inplace=True)
    
    return adata

def main():
    # Paths
    data_path = Path("data/raw/healthy_human_4chamber_map_unnormalized_V3.h5ad")
    output_path = Path("data/processed/")
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Load data
    adata = load_and_preprocess_data(data_path)
    
    # Basic filtering
    adata = basic_filtering(adata)
    
    # Calculate QC metrics
    adata = calculate_qc_metrics(adata)
    
    # Save preprocessed data
    adata.write(output_path / "heart_data_preprocessed.h5ad")
    
    print("Preprocessing complete!")
    print(f"Final dataset shape: {adata.shape}")

if __name__ == "__main__":
    main()