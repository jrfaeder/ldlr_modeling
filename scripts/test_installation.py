#!/usr/bin/env python3
"""Test installation"""

def test_imports():
    """Test all required packages"""
    packages = {
        'bionetgen': 'bionetgen',
        'pandas': 'pandas',
        'numpy': 'numpy',
        'scipy': 'scipy.stats',
        'matplotlib': 'matplotlib.pyplot'
    }
    
    print("Testing package imports...")
    print("-" * 50)
    
    all_ok = True
    for name, import_path in packages.items():
        try:
            __import__(import_path)
            print(f"✓ {name:20s} OK")
        except ImportError as e:
            print(f"✗ {name:20s} FAILED: {e}")
            all_ok = False
    
    print("-" * 50)
    if all_ok:
        print("\n✅ All packages working!")
        return 0
    else:
        print("\n❌ Some packages failed")
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(test_imports())
