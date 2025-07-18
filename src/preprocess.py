import scanpy as sc  
import numpy as np  
  
def preprocess_sc_data(adata_path, n_top_genes=1000):  
    adata = sc.read_h5ad(adata_path)  
    sc.pp.filter_cells(adata, min_genes=200)  
    sc.pp.filter_genes(adata, min_cells=3)  
    sc.pp.normalize_total(adata, target_sum=1e4)  
    sc.pp.log1p(adata)  
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True)  
    adata = adata[:, adata.var.highly_variable]  
    cell_features = adata.X.copy()  # as numpy array  
    try:  
        cell_labels = np.array(adata.obs['cell_type'])  # Replace with your cell type label column  
    except KeyError:  
        cell_labels = np.zeros(cell_features.shape[0])  
    return cell_features, cell_labels  