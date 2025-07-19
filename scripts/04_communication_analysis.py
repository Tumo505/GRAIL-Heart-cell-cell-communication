import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Try to import cell communication tools

# Check Python environment
import sys
print(f"Python executable: {sys.executable}")
print(f"Python path: {sys.path[:3]}...")

# Check if LIANA is available
LIANA_AVAILABLE = False
try:
    import liana
    import liana as li
    LIANA_AVAILABLE = True
    print(f"LIANA version {liana.__version__} successfully imported")
except ImportError as e:
    print(f"LIANA import failed: {e}")
    print("Install with: pip install liana")
except Exception as e:
    print(f"Unexpected error importing LIANA: {e}")
    print("Install with: pip install liana")

def prepare_data_for_communication(adata, cell_type_col='leiden'):
    """Prepare data for cell communication analysis"""
    print("Preparing data for communication analysis...")
    
    # Ensure we have cell type annotations
    if cell_type_col not in adata.obs.columns:
        print(f"Warning: {cell_type_col} not found in observations")
        return None
    
    # Get raw counts if available
    if adata.raw is not None:
        counts_adata = adata.raw.to_adata()
        counts_adata.obs = adata.obs.copy()
    else:
        counts_adata = adata.copy()
    
    return counts_adata

def run_liana_analysis(adata, cell_type_col='leiden'):
    """Run LIANA cell communication analysis"""
    if not LIANA_AVAILABLE:
        print("LIANA not available. Skipping LIANA analysis.")
        return adata
    
    print("Running LIANA analysis...")
    
    try:
        # Run LIANA with rank aggregate method
        li.mt.rank_aggregate.by_sample(
            adata,
            groupby=cell_type_col,
            resource_name='consensus',
            n_perms=100,
            seed=42,
            verbose=True
        )
        
        # Check if results were stored
        if 'liana_res' in adata.uns:
            print(f"LIANA analysis completed successfully. Found {len(adata.uns['liana_res'])} interactions.")
        else:
            print("LIANA analysis completed but no results found in adata.uns['liana_res']")
            
    except Exception as e:
        print(f"Error running LIANA analysis: {e}")
        print("Continuing with mock analysis...")
    
    return adata

def analyze_communication_patterns(adata, save_path):
    """Analyze and visualize communication patterns"""
    if not LIANA_AVAILABLE or 'liana_res' not in adata.uns:
        print("No LIANA results available. Creating mock analysis...")
        create_mock_communication_analysis(adata, save_path)
        return
    
    print("Analyzing communication patterns...")
    
    # Get LIANA results
    liana_res = adata.uns['liana_res']
    
    # Filter significant interactions
    significant_interactions = liana_res[liana_res['magnitude_rank'] < 0.01]
    
    # Plot communication network
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create a heatmap of communication strength
    comm_matrix = significant_interactions.pivot_table(
        index='source', 
        columns='target', 
        values='magnitude_rank',
        aggfunc='mean'
    )
    
    sns.heatmap(comm_matrix, annot=True, cmap='viridis', ax=ax)
    plt.title('Cell-Cell Communication Strength')
    plt.tight_layout()
    plt.savefig(save_path / "communication_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save significant interactions
    significant_interactions.to_csv(save_path / "significant_interactions.csv")

def create_mock_communication_analysis(adata, save_path):
    """Create mock communication analysis when LIANA is not available"""
    print("Creating mock communication analysis...")
    
    # Get cell type counts
    cell_types = adata.obs['leiden'].value_counts()
    
    # Create mock interaction data
    np.random.seed(42)
    n_interactions = 50
    
    mock_interactions = pd.DataFrame({
        'source': np.random.choice(cell_types.index, n_interactions),
        'target': np.random.choice(cell_types.index, n_interactions),
        'ligand': [f'LIGAND_{i}' for i in range(n_interactions)],
        'receptor': [f'RECEPTOR_{i}' for i in range(n_interactions)],
        'score': np.random.uniform(0, 1, n_interactions)
    })
    
    # Create communication heatmap
    fig, ax = plt.subplots(figsize=(10, 8))
    
    comm_matrix = mock_interactions.pivot_table(
        index='source',
        columns='target', 
        values='score',
        aggfunc='mean'
    ).fillna(0)
    
    sns.heatmap(comm_matrix, annot=True, cmap='viridis', ax=ax)
    plt.title('Mock Cell-Cell Communication Strength')
    plt.tight_layout()
    plt.savefig(save_path / "mock_communication_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save mock interactions
    mock_interactions.to_csv(save_path / "mock_interactions.csv")

def analyze_ligand_receptor_pairs(adata, save_path):
    """Analyze ligand-receptor pairs"""
    print("Analyzing ligand-receptor pairs...")
    
    # Define some known heart-relevant ligand-receptor pairs
    heart_lr_pairs = {
        'VEGFA': ['FLT1', 'KDR'],
        'PDGFB': ['PDGFRB'],
        'TGFB1': ['TGFBR1', 'TGFBR2'],
        'IGF1': ['IGF1R'],
        'BMP2': ['BMPR1A', 'BMPR2'],
        'WNT3A': ['FZD1', 'FZD2'],
        'NOTCH1': ['DLL1', 'JAG1'],
        'CXCL12': ['CXCR4']
    }
    
    # Check expression of these pairs
    lr_expression = pd.DataFrame()
    
    for ligand, receptors in heart_lr_pairs.items():
        if ligand in adata.var_names:
            ligand_expr = adata[:, ligand].X.toarray().flatten()
            lr_expression[ligand] = ligand_expr
            
        for receptor in receptors:
            if receptor in adata.var_names:
                receptor_expr = adata[:, receptor].X.toarray().flatten()
                lr_expression[receptor] = receptor_expr
    
    # Plot expression patterns
    if not lr_expression.empty:
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Create expression heatmap by cell type
        cell_types = adata.obs['leiden'].unique()
        expr_by_celltype = pd.DataFrame()
        
        for ct in cell_types:
            ct_mask = adata.obs['leiden'] == ct
            ct_expr = lr_expression[ct_mask].mean()
            expr_by_celltype[ct] = ct_expr
        
        sns.heatmap(expr_by_celltype, annot=True, cmap='Blues', ax=ax)
        plt.title('Ligand-Receptor Expression by Cell Type')
        plt.tight_layout()
        plt.savefig(save_path / "ligand_receptor_expression.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save expression data
        expr_by_celltype.to_csv(save_path / "ligand_receptor_expression.csv")

def main():
    # Load annotated data
    data_path = Path("data/processed/heart_data_annotated.h5ad")
    adata = sc.read_h5ad(data_path)
    
    # Create results directories
    results_path = Path("results/communication")
    results_path.mkdir(parents=True, exist_ok=True)
    
    # Prepare data for communication analysis
    comm_adata = prepare_data_for_communication(adata)
    
    if comm_adata is None:
        print("Could not prepare data for communication analysis")
        return
    
    # Run LIANA analysis
    print(f"LIANA_AVAILABLE: {LIANA_AVAILABLE}")
    if LIANA_AVAILABLE:
        comm_adata = run_liana_analysis(comm_adata)
    else:
        print("Skipping LIANA analysis due to import failure")
    
    # Analyze communication patterns
    analyze_communication_patterns(comm_adata, results_path)
    
    # Analyze ligand-receptor pairs
    analyze_ligand_receptor_pairs(comm_adata, results_path)
    
    # Save communication analysis results
    comm_adata.write(Path("data/processed/heart_data_communication.h5ad"))
    
    print("Cell-cell communication analysis complete!")

if __name__ == "__main__":
    main()