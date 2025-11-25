#!/bin/bash
# setup_project.sh
# Automated setup for LDLR modeling with PyBioNetGen

set -e

echo "=========================================="
echo "LDLR Modeling Setup (PyBioNetGen)"
echo "=========================================="
echo ""

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_ok() {
    echo -e "${GREEN}✓${NC} $1"
}

print_error() {
    echo -e "${RED}✗${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

# Check Python
echo "Checking prerequisites..."
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version | cut -d' ' -f2)
    print_ok "Python 3 found: $PYTHON_VERSION"
else
    print_error "Python 3 not found. Please install Python 3.8+"
    exit 1
fi

# Check Git
if command -v git &> /dev/null; then
    print_ok "Git found"
else
    print_warning "Git not found (optional)"
fi

echo ""
echo "=========================================="
echo "Creating Project Structure"
echo "=========================================="

# Create directories
echo "Creating directories..."
mkdir -p models
mkdir -p data
mkdir -p scripts
mkdir -p results/figures
mkdir -p results/data
mkdir -p docs

print_ok "Directories created"

# Create placeholders
touch results/figures/.gitkeep
touch results/data/.gitkeep

echo ""
echo "=========================================="
echo "Setting Up Python Environment"
echo "=========================================="

# Create virtual environment
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
    print_ok "Virtual environment created"
else
    print_warning "Virtual environment already exists"
fi

# Activate
echo "Activating virtual environment..."
source venv/bin/activate

# Create requirements
echo "Creating requirements.txt..."
cat > requirements.txt << EOL
bionetgen>=0.8.0
pandas>=1.3.0
numpy>=1.21.0
scipy>=1.7.0
matplotlib>=3.4.0
seaborn>=0.11.0
EOL

print_ok "requirements.txt created"

# Install packages
echo "Installing packages (this may take a minute)..."
pip install --upgrade pip --quiet
pip install -r requirements.txt --quiet

print_ok "Packages installed"

# Test PyBioNetGen
echo "Testing PyBioNetGen..."
python3 << EOPYTHON
import bionetgen.core.version as bng_version
print(f"PyBioNetGen version: {bng_version.get_version()}")
EOPYTHON

print_ok "PyBioNetGen working!"

echo ""
echo "=========================================="
echo "Creating Project Files"
echo "=========================================="

# .gitignore
echo "Creating .gitignore..."
cat > .gitignore << 'EOL'
# Python
__pycache__/
*.py[cod]
venv/
*.egg-info/

# BioNetGen outputs
*.net
*.gdat
*.cdat
*.rxn

# Results
results/data/*.gdat
results/data/*.cdat

# OS
.DS_Store
Thumbs.db

# IDE
.vscode/
.idea/
EOL

print_ok ".gitignore created"

# README
echo "Creating README.md..."
cat > README.md << 'EOL'
# LDLR Functional Landscape Modeling

Rule-based computational model of LDLR for understanding familial hypercholesterolemia.

## Quick Start

```bash
# Setup (first time)
./setup_project.sh

# Activate environment
source venv/bin/activate

# Run simulations
python scripts/run_variants.py

# Analyze
python scripts/analyze_results.py
```

## Technology

- **PyBioNetGen**: Rule-based modeling in Python
- **Python 3.8+**: Core language
- **Pandas/NumPy**: Data analysis
- **Matplotlib**: Visualization

## Project Structure

```
.
├── models/          # Model code
├── scripts/         # Analysis scripts
├── results/         # Figures and data
└── docs/            # Documentation
```

## Citation

Based on: Tabet et al., Science 2025 - "The functional landscape of coding 
variation in the familial hypercholesterolemia gene LDLR"

## License

MIT License
EOL

print_ok "README.md created"

# Example data
echo "Creating example data..."
cat > data/variants.csv << 'EOL'
variant,uptake_score,domain,clinvar
WT,1.00,WT,WT
C52Y,0.05,LA1,pathogenic
D147N,0.15,LA3,pathogenic
P526L,0.95,EGF-A,benign
T705I,0.92,beta_propeller,benign
EOL

print_ok "Example data created"

# Test script
echo "Creating test script..."
cat > scripts/test_installation.py << 'EOPYTHON'
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
EOPYTHON

chmod +x scripts/test_installation.py
print_ok "Test script created"

echo ""
echo "=========================================="
echo "Testing Installation"
echo "=========================================="

python scripts/test_installation.py

echo ""
echo "=========================================="
echo "Git Repository"
echo "=========================================="

if [ ! -d ".git" ]; then
    echo "Initializing git..."
    git init
    git add .
    git commit -m "Initial setup with PyBioNetGen"
    print_ok "Git initialized"
    
    echo ""
    echo "To connect to GitHub:"
    echo "  1. Create repo: https://github.com/new"
    echo "  2. Run:"
    echo "     git remote add origin https://github.com/YOUR_USERNAME/ldlr-modeling.git"
    echo "     git branch -M main"
    echo "     git push -u origin main"
else
    print_warning "Git repository already exists"
fi

echo ""
echo "=========================================="
echo "✅ Setup Complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo ""
echo "1. Follow Day 1 guide to create your first model"
echo "   See: docs/DAY1_GUIDE.md"
echo ""
echo "2. Activate environment:"
echo "   source venv/bin/activate"
echo ""
echo "3. Start coding! 🚀"
echo ""
echo "Helpful commands:"
echo "  python scripts/test_installation.py  # Test everything works"
echo "  bionetgen --version  # Check version"
echo ""
EOPYTHON

if __name__ == "__main__":
    import sys
    sys.exit(test_imports())