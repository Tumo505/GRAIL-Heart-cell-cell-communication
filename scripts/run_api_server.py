#!/usr/bin/env python
"""
API server example
"""

import sys
from pathlib import Path

# Add src to path for development
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / 'src'))

from heartmap.api import create_api

def main():
    try:
        # Create API instance
        api = create_api(config_path="config.yaml")
        
        # Run server
        print("Starting HeartMAP API server...")
        print("Access at: http://localhost:8000")
        print("API docs at: http://localhost:8000/docs")
        
        api.run(host="0.0.0.0", port=8000, debug=True)
        
    except ImportError as e:
        print("FastAPI dependencies not available.")
        print("Install with: pip install -e .[api]")
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
