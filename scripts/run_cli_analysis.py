#!/usr/bin/env python
"""
Command line analysis example
"""

import sys
from pathlib import Path

# Add src to path for development
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / 'src'))

from heartmap.api import CLIInterface

def main():
    cli = CLIInterface()
    
    # Example: replace with your data path
    data_path = "data/raw/healthy_human_4chamber_map_unnormalized_V3.h5ad"
    
    if Path(data_path).exists():
        cli.run_analysis(
            data_path=data_path,
            analysis_type="comprehensive",
            output_dir="results/comprehensive",
            config_path="config.yaml"
        )
    else:
        print(f"Data file not found: {data_path}")
        print("Please update the data_path variable with your actual data file.")

if __name__ == "__main__":
    main()
