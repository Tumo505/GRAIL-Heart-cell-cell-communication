"""
Configuration management for HeartMAP
"""

from pathlib import Path
from typing import Dict, Any, Optional, Union
import yaml
import json
from dataclasses import dataclass, asdict


@dataclass
class DataConfig:
    """Data processing configuration"""
    min_genes: int = 200
    min_cells: int = 3
    max_cells_subset: Optional[int] = None
    max_genes_subset: Optional[int] = None
    target_sum: float = 1e4
    n_top_genes: int = 2000
    random_seed: int = 42
    test_mode: bool = False


@dataclass
class AnalysisConfig:
    """Analysis configuration"""
    n_components_pca: int = 50
    n_neighbors: int = 10
    n_pcs: int = 40
    resolution: float = 0.5
    n_marker_genes: int = 25
    use_leiden: bool = True
    use_liana: bool = True


@dataclass
class ModelConfig:
    """Model configuration"""
    model_type: str = "comprehensive"
    save_intermediate: bool = True
    use_gpu: bool = False
    batch_size: Optional[int] = None
    max_memory_gb: Optional[float] = None


@dataclass
class PathConfig:
    """Path configuration"""
    data_dir: str = "data"
    raw_data_dir: str = "data/raw"
    processed_data_dir: str = "data/processed"
    results_dir: str = "results"
    figures_dir: str = "figures"
    # models_dir: str = "models"


@dataclass
class Config:
    """Main configuration class"""
    data: DataConfig
    analysis: AnalysisConfig
    model: ModelConfig
    paths: PathConfig

    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'Config':
        """Create config from dictionary"""
        return cls(
            data=DataConfig(**config_dict.get('data', {})),
            analysis=AnalysisConfig(**config_dict.get('analysis', {})),
            model=ModelConfig(**config_dict.get('model', {})),
            paths=PathConfig(**config_dict.get('paths', {}))
        )

    @classmethod
    def from_yaml(cls, yaml_path: str) -> 'Config':
        """Load config from YAML file"""
        with open(yaml_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        return cls.from_dict(config_dict)

    @classmethod
    def from_json(cls, json_path: str) -> 'Config':
        """Load config from JSON file"""
        with open(json_path, 'r') as f:
            config_dict = json.load(f)
        return cls.from_dict(config_dict)

    @classmethod
    def default(cls) -> 'Config':
        """Create default configuration"""
        return cls(
            data=DataConfig(),
            analysis=AnalysisConfig(),
            model=ModelConfig(),
            paths=PathConfig()
        )

    @classmethod
    def from_file(cls, file_path: str) -> 'Config':
        """Load config from file (YAML or JSON)"""
        return load_config(file_path)

    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary"""
        return asdict(self)

    def save_yaml(self, yaml_path: str) -> None:
        """Save config to YAML file"""
        with open(yaml_path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False)

    def save_json(self, json_path: str) -> None:
        """Save config to JSON file"""
        with open(json_path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    def update_paths(self, base_dir: str) -> None:
        """Update all paths relative to base directory"""
        base_path = Path(base_dir)
        self.paths.data_dir = str(base_path / "data")
        self.paths.raw_data_dir = str(base_path / "data" / "raw")
        self.paths.processed_data_dir = str(base_path / "data" / "processed")
        self.paths.results_dir = str(base_path / "results")
        self.paths.figures_dir = str(base_path / "figures")
        # self.paths.models_dir = str(base_path / "models")

    def create_directories(self) -> None:
        """Create all configured directories"""
        dirs_to_create = [
            self.paths.data_dir,
            self.paths.raw_data_dir,
            self.paths.processed_data_dir,
            self.paths.results_dir,
            self.paths.figures_dir,
            # self.paths.models_dir
        ]

        for dir_path in dirs_to_create:
            Path(dir_path).mkdir(parents=True, exist_ok=True)


def load_config(config_path: Optional[Union[str, Path]] = None) -> Config:
    """Load configuration from file or return default"""
    if config_path is None:
        return Config.default()

    config_path = Path(config_path)
    if config_path.suffix.lower() in ['.yml', '.yaml']:
        return Config.from_yaml(str(config_path))
    elif config_path.suffix.lower() == '.json':
        return Config.from_json(str(config_path))
    else:
        raise ValueError(f"Unsupported config file format: {config_path.suffix}")


# Export configuration classes
__all__ = [
    'DataConfig',
    'AnalysisConfig',
    'ModelConfig',
    'PathConfig',
    'Config',
    'load_config'
]
