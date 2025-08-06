"""
Analysis pipelines for HeartMAP
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List, Tuple
import warnings
from pathlib import Path

try:
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import anndata as ad
    import matplotlib.pyplot as plt
    import seaborn as sns
    DEPS_AVAILABLE = True
except ImportError:
    DEPS_AVAILABLE = False
    warnings.warn("Some dependencies not available. Install requirements for full functionality.")

from ..config import Config
from ..data import DataProcessor
from ..models import CellAnnotationModel, CommunicationModel, MultiChamberModel, HeartMapModel
from ..utils import Visualizer, ResultsExporter


class BasePipeline(ABC):
    """Base class for analysis pipelines"""
    
    def __init__(self, config: Config):
        self.config = config
        self.data_processor = DataProcessor(config)
        self.visualizer = Visualizer(config)
        self.exporter = ResultsExporter(config)
        self.results = {}
    
    @abstractmethod
    def run(self, data_path: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """Run the complete pipeline"""
        pass
    
    def save_results(self, output_dir: str) -> None:
        """Save pipeline results"""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        self.exporter.export_results(self.results, output_path)


class BasicPipeline(BasePipeline):
    """Basic single-cell analysis pipeline"""
    
    def __init__(self, config: Config):
        super().__init__(config)
        self.model = CellAnnotationModel(config)
    
    def run(self, data_path: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """Run basic analysis pipeline"""
        if not DEPS_AVAILABLE:
            raise ImportError("Required dependencies not available")
        
        print("=== Running Basic Pipeline ===")
        
        # Load and process data
        print("1. Loading and processing data...")
        adata = self.data_processor.process_from_raw(data_path)
        
        # Fit annotation model
        print("2. Performing cell annotation...")
        self.model.fit(adata)
        
        # Get results
        results = self.model.predict(adata)
        adata.obs['leiden'] = results['cluster_labels']
        
        # Generate visualizations
        print("3. Generating visualizations...")
        if output_dir:
            viz_dir = Path(output_dir) / "figures"
            viz_dir.mkdir(parents=True, exist_ok=True)
            
            # UMAP plot
            sc.pl.umap(adata, color=['leiden'], legend_loc='on data',
                      title='Cell Type Clusters', show=False)
            plt.savefig(viz_dir / "umap_clusters.png", dpi=300, bbox_inches='tight')
            plt.close()
            
            # QC metrics
            self.visualizer.plot_qc_metrics(adata, viz_dir)
        
        # Store results
        self.results = {
            'adata': adata,
            'model': self.model,
            'results': results
        }
        
        # Save results
        if output_dir:
            self.save_results(output_dir)
            # Save processed data
            adata.write(Path(output_dir) / "annotated_data.h5ad")
            # Save model
            self.model.save(Path(output_dir) / "annotation_model.pkl")
        
        print("Basic pipeline completed!")
        return self.results


class AdvancedCommunicationPipeline(BasePipeline):
    """Advanced cell-cell communication analysis pipeline"""
    
    def __init__(self, config: Config):
        super().__init__(config)
        self.model = CommunicationModel(config)
    
    def run(self, data_path: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """Run advanced communication analysis pipeline"""
        if not DEPS_AVAILABLE:
            raise ImportError("Required dependencies not available")
        
        print("=== Running Advanced Communication Pipeline ===")
        
        # Load processed data (should have cell annotations)
        print("1. Loading annotated data...")
        adata = sc.read_h5ad(data_path)
        
        if 'leiden' not in adata.obs.columns:
            raise ValueError("Input data must have cell type annotations. Run BasicPipeline first.")
        
        # Fit communication model
        print("2. Analyzing cell-cell communication...")
        self.model.fit(adata)
        
        # Get results
        results = self.model.predict(adata)
        
        # Generate visualizations
        print("3. Generating communication visualizations...")
        if output_dir:
            viz_dir = Path(output_dir) / "figures"
            viz_dir.mkdir(parents=True, exist_ok=True)
            
            self.visualizer.plot_communication_heatmap(
                results['communication_scores'], viz_dir
            )
            self.visualizer.plot_hub_scores(
                adata, results['hub_scores'], viz_dir
            )
            self.visualizer.plot_pathway_scores(
                results['pathway_scores'], viz_dir
            )
        
        # Store results
        self.results = {
            'adata': adata,
            'model': self.model,
            'results': results
        }
        
        # Save results
        if output_dir:
            self.save_results(output_dir)
            self.model.save(Path(output_dir) / "communication_model.pkl")
        
        print("Advanced communication pipeline completed!")
        return self.results


class MultiChamberPipeline(BasePipeline):
    """Multi-chamber heart analysis pipeline"""
    
    def __init__(self, config: Config):
        super().__init__(config)
        self.model = MultiChamberModel(config)
    
    def run(self, data_path: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """Run multi-chamber analysis pipeline"""
        if not DEPS_AVAILABLE:
            raise ImportError("Required dependencies not available")
        
        print("=== Running Multi-Chamber Pipeline ===")
        
        # Load data
        print("1. Loading data...")
        adata = sc.read_h5ad(data_path)
        
        # Fit multi-chamber model
        print("2. Analyzing multi-chamber patterns...")
        self.model.fit(adata)
        
        # Get results
        results = self.model.predict(adata)
        
        # Generate visualizations
        print("3. Generating multi-chamber visualizations...")
        if output_dir:
            viz_dir = Path(output_dir) / "figures"
            viz_dir.mkdir(parents=True, exist_ok=True)
            
            self.visualizer.plot_chamber_composition(
                adata, viz_dir
            )
            self.visualizer.plot_chamber_markers(
                results['chamber_markers'], viz_dir
            )
            self.visualizer.plot_cross_chamber_correlations(
                results['cross_chamber_correlations'], viz_dir
            )
        
        # Store results
        self.results = {
            'adata': adata,
            'model': self.model,
            'results': results
        }
        
        # Save results
        if output_dir:
            self.save_results(output_dir)
            self.model.save(Path(output_dir) / "multi_chamber_model.pkl")
        
        print("Multi-chamber pipeline completed!")
        return self.results


class ComprehensivePipeline(BasePipeline):
    """Comprehensive HeartMAP analysis pipeline"""
    
    def __init__(self, config: Config):
        super().__init__(config)
        self.model = HeartMapModel(config)
    
    def run(self, data_path: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """Run comprehensive HeartMAP analysis"""
        if not DEPS_AVAILABLE:
            raise ImportError("Required dependencies not available")
        
        print("=== Running Comprehensive HeartMAP Pipeline ===")
        
        # Load and process data
        print("1. Loading and processing data...")
        adata = self.data_processor.process_from_raw(data_path)
        
        # Fit complete model
        print("2. Fitting comprehensive HeartMAP model...")
        self.model.fit(adata)
        
        # Get results
        results = self.model.predict(adata)
        
        # Update adata with all results
        adata.obs['leiden'] = results['annotation']['cluster_labels']
        adata.obs['hub_score'] = results['communication']['hub_scores']
        
        # Generate comprehensive visualizations
        print("3. Generating comprehensive visualizations...")
        if output_dir:
            viz_dir = Path(output_dir) / "figures"
            viz_dir.mkdir(parents=True, exist_ok=True)
            
            # Create comprehensive dashboard
            self.visualizer.create_comprehensive_dashboard(adata, results, viz_dir)
        
        # Store results
        self.results = {
            'adata': adata,
            'model': self.model,
            'results': results
        }
        
        # Save results
        if output_dir:
            self.save_results(output_dir)
            adata.write(Path(output_dir) / "heartmap_complete.h5ad")
            self.model.save(Path(output_dir) / "heartmap_model.pkl")
            
            # Generate comprehensive report
            self.exporter.generate_comprehensive_report(self.results, output_dir)
        
        print("Comprehensive HeartMAP pipeline completed!")
        return self.results
