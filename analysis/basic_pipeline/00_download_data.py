import os
import requests
from pathlib import Path

def download_data():
    """Download the heart dataset"""
    # Use project root as reference
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "data/raw"
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Note: You'll need to download this manually from the website
    # as it requires registration
    url = "https://singlecell.broadinstitute.org/single_cell/study/SCP498"
    
    print(f"Please download the dataset manually from: {url}")
    print("Look for: healthy_human_4chamber_map_unnormalized_V3.h5ad")
    print(f"Save it to: {data_dir}/healthy_human_4chamber_map_unnormalized_V3.h5ad")

if __name__ == "__main__":
    download_data()