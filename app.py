#!/usr/bin/env python
"""
Gradio app for Hugging Face Spaces deployment
"""

import gradio as gr
import tempfile
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, 'src')

try:
    from heartmap import Config
    from heartmap.pipelines import ComprehensivePipeline
    HEARTMAP_AVAILABLE = True
except ImportError:
    HEARTMAP_AVAILABLE = False

def analyze_heart_data(uploaded_file, analysis_type, max_cells, max_genes):
    """Analyze uploaded heart data"""
    
    if not HEARTMAP_AVAILABLE:
        return "HeartMAP not available. Please install dependencies."
    
    if uploaded_file is None:
        return "Please upload a file."
    
    try:
        # Create config
        config = Config.default()
        config.data.max_cells_subset = int(max_cells) if max_cells else None
        config.data.max_genes_subset = int(max_genes) if max_genes else None
        
        # Create pipeline
        if analysis_type == "comprehensive":
            pipeline = ComprehensivePipeline(config)
        else:
            return f"Analysis type {analysis_type} not implemented in demo"
        
        # For demo, just return configuration info
        return f"""Analysis configured successfully!

Configuration:
- Analysis type: {analysis_type}
- Max cells: {max_cells}
- Max genes: {max_genes}
- File: {uploaded_file.name}

Note: This is a demo interface. 
In the full version, your data would be processed here.
"""
        
    except Exception as e:
        return f"Error: {str(e)}"

# Create Gradio interface
demo = gr.Interface(
    fn=analyze_heart_data,
    inputs=[
        gr.File(label="Upload single-cell data (.h5ad)", file_types=[".h5ad"]),
        gr.Dropdown(
            choices=["comprehensive", "basic", "advanced", "multi_chamber"],
            value="comprehensive",
            label="Analysis Type"
        ),
        gr.Number(label="Max Cells (for memory optimization)", value=50000, precision=0),
        gr.Number(label="Max Genes (for memory optimization)", value=5000, precision=0),
    ],
    outputs=gr.Textbox(label="Analysis Results"),
    title="HeartMAP: Heart Multi-chamber Analysis Platform",
    description="""
    Upload single-cell RNA-seq data from heart tissue for comprehensive analysis.
    
    **Features:**
    - Cell type annotation
    - Cell-cell communication analysis  
    - Multi-chamber analysis
    - Chamber-specific marker identification
    
    **Input format:** AnnData (.h5ad) files
    
    **Demo Note:** This demo shows the interface. Full analysis requires local installation.
    """,
    examples=[
        [None, "comprehensive", 10000, 2000],
        [None, "basic", 5000, 1000],
    ]
)

if __name__ == "__main__":
    demo.launch()
