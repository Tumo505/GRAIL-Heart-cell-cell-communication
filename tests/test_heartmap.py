"""
Tests for HeartMAP functionality
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import numpy as np

try:
    import scanpy as sc
    import pandas as pd
    import anndata as ad
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False

import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from heartmap.config import Config
from heartmap.data import DataProcessor


class TestConfig(unittest.TestCase):
    """Test configuration functionality"""
    
    def test_default_config(self):
        """Test default configuration creation"""
        config = Config.default()
        self.assertIsInstance(config.data.min_genes, int)
        self.assertIsInstance(config.analysis.resolution, float)
    
    def test_config_dict_conversion(self):
        """Test config to/from dict conversion"""
        config = Config.default()
        config_dict = config.to_dict()
        
        self.assertIsInstance(config_dict, dict)
        self.assertIn('data', config_dict)
        self.assertIn('analysis', config_dict)
        
        # Test loading from dict
        new_config = Config.from_dict(config_dict)
        self.assertEqual(config.data.min_genes, new_config.data.min_genes)
    
    def test_config_file_operations(self):
        """Test config file save/load"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            config = Config.default()
            yaml_path = Path(tmp_dir) / "test.yaml"
            json_path = Path(tmp_dir) / "test.json"
            
            # Save
            config.save_yaml(str(yaml_path))
            config.save_json(str(json_path))
            
            # Load
            if yaml_path.exists():
                loaded_yaml = Config.from_yaml(str(yaml_path))
                self.assertEqual(config.data.min_genes, loaded_yaml.data.min_genes)
            
            if json_path.exists():
                loaded_json = Config.from_json(str(json_path))
                self.assertEqual(config.data.min_genes, loaded_json.data.min_genes)


class TestDataProcessor(unittest.TestCase):
    """Test data processing functionality"""
    
    @unittest.skipUnless(SCANPY_AVAILABLE, "Scanpy not available")
    def test_create_test_dataset(self):
        """Test test dataset creation"""
        # Create synthetic data
        n_cells, n_genes = 1000, 500
        X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes))
        adata = ad.AnnData(X)
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
        adata.var_names = [f"gene_{i}" for i in range(n_genes)]
        
        config = Config.default()
        processor = DataProcessor(config)
        
        # Create test subset
        test_adata = processor.create_test_dataset(adata, n_cells=100)
        
        self.assertEqual(test_adata.n_obs, 100)
        self.assertEqual(test_adata.n_vars, n_genes)
    
    @unittest.skipUnless(SCANPY_AVAILABLE, "Scanpy not available")
    def test_basic_preprocessing(self):
        """Test basic preprocessing"""
        # Create synthetic data
        n_cells, n_genes = 1000, 500
        X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes))
        adata = ad.AnnData(X)
        adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
        adata.var_names = [f"gene_{i}" for i in range(n_genes)]
        
        config = Config.default()
        processor = DataProcessor(config)
        
        # Basic preprocessing
        processed = processor.loader.preprocess_basic(adata)
        
        self.assertIsNotNone(processed.raw)
        self.assertTrue(processed.var_names.is_unique)


class TestPipelines(unittest.TestCase):
    """Test pipeline functionality"""
    
    def test_pipeline_creation(self):
        """Test pipeline creation"""
        from heartmap.pipelines import BasicPipeline
        
        config = Config.default()
        pipeline = BasicPipeline(config)
        
        self.assertIsNotNone(pipeline.config)
        self.assertIsNotNone(pipeline.data_processor)


class TestAPI(unittest.TestCase):
    """Test API functionality"""
    
    def test_cli_interface_creation(self):
        """Test CLI interface creation"""
        from heartmap.api import CLIInterface
        
        cli = CLIInterface()
        self.assertIsNotNone(cli)


def run_tests():
    """Run all tests"""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    test_classes = [
        TestConfig,
        TestDataProcessor,
        TestPipelines,
        TestAPI
    ]
    
    for test_class in test_classes:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    if not success:
        exit(1)
