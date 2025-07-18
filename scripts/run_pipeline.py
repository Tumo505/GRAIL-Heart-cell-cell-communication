import subprocess
import sys
from pathlib import Path
import time

def run_script(script_path, description):
    """Run a Python script and handle errors"""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Script: {script_path}")
    print(f"{'='*60}")
    
    try:
        start_time = time.time()
        result = subprocess.run([sys.executable, script_path], 
                              capture_output=True, text=True, check=True)
        end_time = time.time()
        
        print(f"✅ Completed in {end_time - start_time:.2f} seconds")
        if result.stdout:
            print("Output:", result.stdout)
        
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running {script_path}")
        print(f"Return code: {e.returncode}")
        print(f"Error output: {e.stderr}")
        return False

def main():
    """Run the complete cell communication analysis pipeline"""
    print("Starting Cell-Cell Communication Analysis Pipeline")
    print("=" * 60)
    
    # Define pipeline steps
    pipeline_steps = [
        ("scripts/00_download_data.py", "Data Download Instructions"),
        ("scripts/01_data_preprocessing.py", "Data Preprocessing"),
        ("scripts/02_quality_control.py", "Quality Control"),
        ("scripts/03_cell_annotation.py", "Cell Type Annotation"),
        ("scripts/04_communication_analysis.py", "Cell Communication Analysis"),
        ("scripts/05_visualization.py", "Visualization and Reporting"),
        ("scripts/06_advanced_communication.py", "Advanced Communication Analysis")
    ]
    
    # Check if data file exists
    data_file = Path("data/raw/healthy_human_4chamber_map_unnormalized_V3.h5ad")
    if not data_file.exists():
        print(f"⚠️  Data file not found: {data_file}")
        print("Please download the dataset first using:")
        print("python scripts/00_download_data.py")
        return
    
    # Run pipeline steps
    successful_steps = 0
    
    for script_path, description in pipeline_steps[1:]:  # Skip download step
        if run_script(script_path, description):
            successful_steps += 1
        else:
            print(f"\n❌ Pipeline failed at step: {description}")
            break
    
    # Summary
    print(f"\n{'='*60}")
    print("PIPELINE SUMMARY")
    print(f"{'='*60}")
    print(f"Successfully completed: {successful_steps}/{len(pipeline_steps)-1} steps")
    
    if successful_steps == len(pipeline_steps) - 1:
        print("🎉 Complete pipeline executed successfully!")
        print("\nResults are available in:")
        print("- results/figures/")
        print("- results/communication/")
        print("- results/advanced/")
        print("- results/final/")
    else:
        print("⚠️  Pipeline completed with errors. Check the output above.")

if __name__ == "__main__":
    main()