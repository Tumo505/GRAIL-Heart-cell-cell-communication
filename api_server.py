#!/usr/bin/env python3
"""
FastAPI server for HeartMAP
Run with: uvicorn api_server:app --host 0.0.0.0 --port 8000
"""

import os
import sys
from pathlib import Path

# Add src to Python path for development
sys.path.insert(0, str(Path(__file__).parent / "src"))

from heartmap.api import create_api

# Create the FastAPI app
app_instance = create_api()
app = app_instance.app

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000, log_level="info")
