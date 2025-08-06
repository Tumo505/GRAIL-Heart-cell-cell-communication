"""
API interface for HeartMAP
"""

from typing import Dict, Any, Optional, Union, List
import warnings
from pathlib import Path
import tempfile

from ..config import Config, load_config

try:
    from fastapi import FastAPI, UploadFile, File, HTTPException
    from pydantic import BaseModel as PydanticBaseModel
    import uvicorn
    FASTAPI_AVAILABLE = True
except ImportError:
    FASTAPI_AVAILABLE = False
    warnings.warn("FastAPI not available. Install with: pip install fastapi uvicorn")

    # Create a dummy BaseModel for when FastAPI/Pydantic is not available
    class PydanticBaseModel:  # type: ignore
        def __init__(self, **kwargs):
            for key, value in kwargs.items():
                setattr(self, key, value)

from ..pipelines import (
    BasicPipeline,
    AdvancedCommunicationPipeline,
    MultiChamberPipeline,
    ComprehensivePipeline
)


class AnalysisRequest(PydanticBaseModel):
    """Request model for analysis"""
    analysis_type: str = "comprehensive"  # basic, advanced, multi_chamber, comprehensive
    config_overrides: Optional[Dict[str, Any]] = None
    output_format: str = "json"  # json, csv, h5ad


class AnalysisResponse(PydanticBaseModel):
    """Response model for analysis"""
    status: str
    message: str
    results: Optional[Dict[str, Any]] = None
    output_files: Optional[List[str]] = None


class HeartMapAPI:
    """REST API for HeartMAP analysis"""

    def __init__(self, config: Union[str, Config, None] = None):
        if isinstance(config, Config):
            self.config = config
        else:
            self.config = load_config(config)
        self.app = FastAPI(
            title="HeartMAP API",
            description="Heart Multi-chamber Analysis Platform API",
            version="1.0.0"
        ) if FASTAPI_AVAILABLE else None

        if FASTAPI_AVAILABLE:
            self._setup_routes()

    def _setup_routes(self):
        """Setup API routes"""

        @self.app.get("/")
        async def root():
            return {"message": "HeartMAP API", "version": "1.0.0"}

        @self.app.get("/health")
        async def health_check():
            return {"status": "healthy"}

        @self.app.post("/analyze", response_model=AnalysisResponse)
        async def analyze_data(
            file: UploadFile = File(...),
            request: AnalysisRequest = AnalysisRequest()
        ):
            """Analyze single-cell data"""
            try:
                # Save uploaded file temporarily
                with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp_file:
                    content = await file.read()
                    tmp_file.write(content)
                    tmp_file_path = tmp_file.name

                # Update config with overrides
                if request.config_overrides:
                    # Apply config overrides (simplified)
                    pass

                # Create output directory
                with tempfile.TemporaryDirectory() as output_dir:
                    # Run analysis based on type
                    pipeline = self._get_pipeline(request.analysis_type)
                    results = pipeline.run(tmp_file_path, output_dir)

                    # Format response
                    response_data = self._format_results(results, request.output_format)

                    return AnalysisResponse(
                        status="success",
                        message="Analysis completed successfully",
                        results=response_data
                    )

            except Exception as e:
                raise HTTPException(status_code=500, detail=str(e))
            finally:
                # Clean up temporary file
                Path(tmp_file_path).unlink(missing_ok=True)

        @self.app.get("/models")
        async def list_models():
            """List available analysis models"""
            return {
                "models": [
                    "basic",
                    "advanced_communication",
                    "multi_chamber",
                    "comprehensive"
                ]
            }

        @self.app.get("/config")
        async def get_config():
            """Get current configuration"""
            return self.config.to_dict()

        @self.app.post("/config")
        async def update_config(new_config: Dict[str, Any]):
            """Update configuration"""
            try:
                self.config = Config.from_dict(new_config)
                return {"status": "success", "message": "Configuration updated"}
            except Exception as e:
                raise HTTPException(status_code=400, detail=str(e))

    def _get_pipeline(self, analysis_type: str):
        """Get analysis pipeline by type"""
        pipelines = {
            "basic": BasicPipeline(self.config),
            "advanced_communication": AdvancedCommunicationPipeline(self.config),
            "multi_chamber": MultiChamberPipeline(self.config),
            "comprehensive": ComprehensivePipeline(self.config)
        }

        if analysis_type not in pipelines:
            raise ValueError(f"Unknown analysis type: {analysis_type}")

        return pipelines[analysis_type]

    def _format_results(self, results: Dict[str, Any], output_format: str) -> Dict[str, Any]:
        """Format results for API response"""
        if output_format == "json":
            # Convert complex objects to serializable format
            formatted = {}

            if 'results' in results:
                res = results['results']
                formatted['summary'] = {
                    'n_cells': results.get('adata', {}).n_obs if 'adata' in results else 0,
                    'analysis_completed': True
                }

                # Extract key metrics
                if 'annotation' in res:
                    ann_res = res['annotation']
                    if 'metadata' in ann_res:
                        formatted['annotation_summary'] = ann_res['metadata']

                if 'communication' in res:
                    comm_res = res['communication']
                    if 'metadata' in comm_res:
                        formatted['communication_summary'] = comm_res['metadata']

                if 'multi_chamber' in res:
                    chamber_res = res['multi_chamber']
                    if 'metadata' in chamber_res:
                        formatted['multi_chamber_summary'] = chamber_res['metadata']

            return formatted

        else:
            raise ValueError(f"Unsupported output format: {output_format}")

    def run(self, host: str = "0.0.0.0", port: int = 8000, debug: bool = False):
        """Run the API server"""
        if not FASTAPI_AVAILABLE:
            raise ImportError("FastAPI not available. Install with: pip install fastapi uvicorn")

        if self.app is not None:
            uvicorn.run(self.app, host=host, port=port, log_level="debug" if debug else "info")


class CLIInterface:
    """Command line interface for HeartMAP"""

    def __init__(self):
        self.config = None

    def run_analysis(self,
                     data_path: str,
                     analysis_type: str = "comprehensive",
                     output_dir: str = "results",
                     config_path: Optional[str] = None):
        """Run analysis from command line"""

        # Load configuration
        self.config = load_config(config_path)

        # Update output directory
        self.config.update_paths(str(Path(output_dir).parent))
        self.config.create_directories()

        # Get pipeline
        pipelines = {
            "basic": BasicPipeline(self.config),
            "advanced": AdvancedCommunicationPipeline(self.config),
            "multi_chamber": MultiChamberPipeline(self.config),
            "comprehensive": ComprehensivePipeline(self.config)
        }

        if analysis_type not in pipelines:
            raise ValueError(f"Unknown analysis type: {analysis_type}")

        pipeline = pipelines[analysis_type]

        # Run analysis
        print(f"Starting {analysis_type} analysis...")
        results = pipeline.run(data_path, output_dir)

        print(f"Analysis completed! Results saved to: {output_dir}")
        return results


def create_api(config_path: Optional[str] = None) -> HeartMapAPI:
    """Create HeartMAP API instance"""
    return HeartMapAPI(config_path)


def run_cli():
    """Run command line interface"""
    import argparse

    parser = argparse.ArgumentParser(description="HeartMAP Analysis Platform")
    parser.add_argument("data_path", help="Path to input data file")
    parser.add_argument("--analysis-type", default="comprehensive",
                        choices=["basic", "advanced", "multi_chamber", "comprehensive"],
                        help="Type of analysis to run")
    parser.add_argument("--output-dir", default="results",
                        help="Output directory for results")
    parser.add_argument("--config", help="Path to configuration file")

    args = parser.parse_args()

    cli = CLIInterface()
    cli.run_analysis(
        data_path=args.data_path,
        analysis_type=args.analysis_type,
        output_dir=args.output_dir,
        config_path=args.config
    )


# Export API classes and functions
__all__ = [
    'AnalysisRequest',
    'AnalysisResponse',
    'HeartMapAPI',
    'CLIInterface',
    'create_api',
    'run_cli'
]
