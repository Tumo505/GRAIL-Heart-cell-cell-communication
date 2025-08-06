#!/usr/bin/env python3
"""
Test script to verify all HeartMAP imports work correctly
This script is designed to be run in CI environments
"""

import sys
import traceback

def test_import(module_path, description):
    """Test importing a module and print result"""
    try:
        if '.' in module_path:
            # For classes, import the module and get the class
            parts = module_path.split('.')
            if len(parts) == 2:
                # Import like: from heartmap.config import Config
                module = __import__(parts[0], fromlist=[parts[1]])
                getattr(module, parts[1])  # Check if the class exists
            else:
                # Just import the module
                __import__(module_path)
        else:
            __import__(module_path)
        print(f"✅ {description}: {module_path}")
        return True
    except ImportError as e:
        print(f"❌ {description}: {module_path} - {e}")
        return False
    except AttributeError as e:
        print(f"❌ {description}: {module_path} - {e}")
        return False
    except Exception as e:
        print(f"⚠️  {description}: {module_path} - {e}")
        return False

def main():
    """Test all critical imports"""
    print("🧪 Testing HeartMAP package imports...")
    
    tests = [
        # Core package
        ("heartmap", "Core package"),
        
        # Configuration
        ("heartmap.config", "Config module"),
        ("heartmap.Config", "Config class"),
        
        # Models
        ("heartmap.models", "Models module"),
        ("heartmap.HeartMapModel", "HeartMapModel class"),
        
        # Pipelines
        ("heartmap.pipelines", "Pipelines module"),
        ("heartmap.BasicPipeline", "BasicPipeline class"),
        ("heartmap.AdvancedCommunicationPipeline", "AdvancedCommunicationPipeline class"),
        ("heartmap.MultiChamberPipeline", "MultiChamberPipeline class"),
        
        # Data processing
        ("heartmap.data", "Data module"),
        
        # Utilities
        ("heartmap.utils", "Utils module"),
        
        # API (may have warnings for optional dependencies)
        ("heartmap.api", "API module"),
    ]
    
    passed = 0
    failed = 0
    
    for module_path, description in tests:
        if test_import(module_path, description):
            passed += 1
        else:
            failed += 1
    
    print(f"\n📊 Results: {passed} passed, {failed} failed")
    
    if failed > 0:
        print("❌ Some imports failed!")
        sys.exit(1)
    else:
        print("✅ All imports successful!")
        sys.exit(0)

if __name__ == "__main__":
    main()
