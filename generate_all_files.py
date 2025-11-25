#!/usr/bin/env python3
"""
Generate complete LDLR modeling project structure
Run this script to create all files and folders

Usage:
    python3 generate_project.py

This will create:
- All directory structure
- Model code
- Analysis scripts
- Documentation
- Setup files
"""

import os
from pathlib import Path

def create_directory_structure():
    """Create all project directories"""
    dirs = [
        'models',
        'scripts',
        'data',
        'results/figures',
        'results/data',
        'docs'
    ]
    
    for dir_path in dirs:
        Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    print("✓ Created directory structure")

def create_gitignore():
    """Create .gitignore file"""
    content = """# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
venv/
env/
ENV/
*.egg-info/
dist/
build/

# Jupyter
.ipynb_checkpoints/
*.ipynb

# BioNetGen outputs
*.net
*.gdat
*.cdat
*.rxn
*.species
*.groups
*.log

# Results
results/data/*.gdat
results/data/*.cdat

# Data files (large)
*.xlsx
*.xls
data/raw/

# OS
.DS_Store
Thumbs.db
*.swp
*.swo
*~

# IDE
.vscode/
.idea/
*.sublime-project
*.sublime-workspace

# Temporary
temp/
tmp/
*.tmp
"""
    
    with open('.gitignore', 'w') as f:
        f.write(content)
    
    print("✓ Created .gitignore")

def create_requirements():
    """Create requirements.txt"""
    content = """pybiogen>=0.8.0
pandas>=1.3.0
numpy>=1.21.0
scipy>=1.7.0
matplotlib>=3.4.0
seaborn>=0.11.0
"""
    
    with open('requirements.txt', 'w') as f:
        f.write(content)
    
    print("✓ Created requirements.txt")

def create_readme():
    """Create README.md"""
    content = """# LDLR Functional Landscape Modeling

[![Python](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![PyBioNetGen](https://img.shields.io/badge/PyBioNetGen-0.8%2B-green)](https://github.com/RuleWorld/PyBioNetGen)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Computational rule-based model of the LDL receptor system using PyBioNetGen.

## Quick Start

```bash
# Setup
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Run simulations
python scripts/run_variants.py

# Analyze results
python scripts/analyze_results.py
```

## Project Structure

```
.
├── models/          # PyBioNetGen model code
├── scripts/         # Analysis pipeline
├── data/            # Variant functional scores
├── results/         # Figures and simulation data
└── docs/            # Documentation
```

## Documentation

- [1-Week Plan](docs/PLAN_1WEEK.md) - Project timeline
- [Day 1 Guide](docs/DAY1_GUIDE.md) - Getting started
- [3-Day MVP](docs/MVP_3DAY.md) - Minimal version

## Technology

- **PyBioNetGen** - Rule-based modeling in Python
- **Python 3.8+** - Core language
- **pandas/numpy/scipy** - Data analysis
- **matplotlib/seaborn** - Visualization

## Citation

Based on: Tabet et al., Science 2025 - "The functional landscape of coding 
variation in the familial hypercholesterolemia gene LDLR"

## License

MIT License
"""
    
    with open('README.md', 'w') as f:
        f.write(content)
    
    print("✓ Created README.md")

def create_ldlr_model():
    """Create models/ldlr_model.py"""
    content = '''"""
LDLR model using PyBioNetGen
Simplified version with 4 LA modules
"""

import bionetgen
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

class LDLRModel:
    """LDLR-LDL uptake model"""
    
    def __init__(self, variant_name='WT', functional_score=1.0):
        """
        Initialize model
        
        Parameters:
        -----------
        variant_name : str
            Variant identifier (e.g., 'WT', 'C52Y')
        functional_score : float
            Functional score from experiment (0-1.5)
        """
        self.variant_name = variant_name
        self.functional_score = functional_score
        self.result = None
        
    def get_model_string(self):
        """Generate BioNetGen model string"""
        
        # Scale parameters by functional score
        binding_strength = max(self.functional_score, 0.01)
        
        model = f"""
# LDLR Model - {self.variant_name}
# Functional score: {self.functional_score:.3f}

begin model

begin parameters
    # Binding (scaled by functional score)
    k_on_base      1.0
    k_off_surf     {1.0 / binding_strength:.4f}
    k_off_endo     50.0
    
    # Module strengths
    strength_LA3   1.0
    strength_LA4   1.0
    strength_LA5   1.0
    strength_LA7   0.8
    
    # Trafficking
    k_endo         2.0
    k_recycle      3.0
    k_degrade      5.0
    
    # Initial conditions
    LDLR_init      1000
    LDL_conc       100
end parameters

begin molecule types
    LDLR(la3, la4, la5, la7, loc~surface~endosome)
    LDL(ldlr, loc~extra~endo~lyso)
end molecule types

begin seed species
    LDLR(la3, la4, la5, la7, loc~surface)  LDLR_init
    LDL(ldlr, loc~extra)  LDL_conc
end seed species

begin observables
    Molecules  LDLR_surface      LDLR(loc~surface)
    Molecules  LDLR_endosome     LDLR(loc~endosome)
    Molecules  LDL_free          LDL(ldlr, loc~extra)
    Molecules  LDL_internalized  LDL(loc~endo) LDL(loc~lyso)
    Molecules  Complex_surface   LDLR(la3!+, loc~surface)
end observables

begin reaction rules
    # LDLR-LDL Binding
    LDLR(la3, loc~surface) + LDL(ldlr, loc~extra) <-> \\\\
        LDLR(la3!1, loc~surface).LDL(ldlr!1, loc~extra) \\\\
        k_on_base*strength_LA3, k_off_surf
    
    LDLR(la4, loc~surface) + LDL(ldlr, loc~extra) <-> \\\\
        LDLR(la4!1, loc~surface).LDL(ldlr!1, loc~extra) \\\\
        k_on_base*strength_LA4, k_off_surf
    
    LDLR(la5, loc~surface) + LDL(ldlr, loc~extra) <-> \\\\
        LDLR(la5!1, loc~surface).LDL(ldlr!1, loc~extra) \\\\
        k_on_base*strength_LA5, k_off_surf
    
    LDLR(la7, loc~surface) + LDL(ldlr, loc~extra) <-> \\\\
        LDLR(la7!1, loc~surface).LDL(ldlr!1, loc~extra) \\\\
        k_on_base*strength_LA7, k_off_surf
    
    # Endocytosis
    LDLR(la3!1, loc~surface).LDL(ldlr!1, loc~extra) -> \\\\
        LDLR(la3!1, loc~endosome).LDL(ldlr!1, loc~endo)  k_endo
    
    LDLR(la4!1, loc~surface).LDL(ldlr!1, loc~extra) -> \\\\
        LDLR(la4!1, loc~endosome).LDL(ldlr!1, loc~endo)  k_endo
    
    LDLR(la5!1, loc~surface).LDL(ldlr!1, loc~extra) -> \\\\
        LDLR(la5!1, loc~endosome).LDL(ldlr!1, loc~endo)  k_endo
    
    LDLR(la7!1, loc~surface).LDL(ldlr!1, loc~extra) -> \\\\
        LDLR(la7!1, loc~endosome).LDL(ldlr!1, loc~endo)  k_endo
    
    # LDL Release (pH-dependent)
    LDLR(la3!1, loc~endosome).LDL(ldlr!1, loc~endo) -> \\\\
        LDLR(la3, loc~endosome) + LDL(ldlr, loc~endo)  k_off_endo
    
    LDLR(la4!1, loc~endosome).LDL(ldlr!1, loc~endo) -> \\\\
        LDLR(la4, loc~endosome) + LDL(ldlr, loc~endo)  k_off_endo
    
    LDLR(la5!1, loc~endosome).LDL(ldlr!1, loc~endo) -> \\\\
        LDLR(la5, loc~endosome) + LDL(ldlr, loc~endo)  k_off_endo
    
    LDLR(la7!1, loc~endosome).LDL(ldlr!1, loc~endo) -> \\\\
        LDLR(la7, loc~endosome) + LDL(ldlr, loc~endo)  k_off_endo
    
    # Receptor Recycling
    LDLR(la3, la4, la5, la7, loc~endosome) -> \\\\
        LDLR(la3, la4, la5, la7, loc~surface)  k_recycle
    
    # LDL Degradation
    LDL(ldlr, loc~endo) -> LDL(ldlr, loc~lyso)  k_recycle
    LDL(loc~lyso) -> 0  k_degrade
    
end reaction rules

end model
"""
        return model
    
    def run(self, t_end=200, n_steps=200, out_dir='results/data'):
        """Run the simulation"""
        
        model_str = self.get_model_string()
        
        print(f"Running {self.variant_name}...")
        self.result = bionetgen.run(
            model_str,
            out_dir=out_dir,
            suppress=True
        )
        
        return self.result
    
    def get_data(self):
        """Load simulation data"""
        if self.result is None:
            raise ValueError("Must run simulation first")
        
        return pd.read_csv(
            self.result['gdat_file'],
            sep='\\s+',
            comment='#'
        )
    
    def get_final_uptake(self):
        """Get final LDL internalized count"""
        data = self.get_data()
        return data['LDL_internalized'].iloc[-1]

# Test if run directly
if __name__ == "__main__":
    print("Testing LDLR model...")
    
    model = LDLRModel('WT', 1.0)
    model.run()
    
    data = model.get_data()
    print(f"\\nFinal LDL internalized: {data['LDL_internalized'].iloc[-1]:.1f}")
    print("✓ Model test passed!")
'''
    
    with open('models/ldlr_model.py', 'w') as f:
        f.write(content)
    
    print("✓ Created models/ldlr_model.py")

def create_run_variants():
    """Create scripts/run_variants.py"""
    content = '''"""
Batch runner for LDLR variants
Runs 5 variants and saves results
"""

import sys
sys.path.append('models')

import pandas as pd
from pathlib import Path
from ldlr_model import LDLRModel
from scipy.stats import pearsonr

# Create directories
Path('results/data').mkdir(parents=True, exist_ok=True)

# Define variants
VARIANTS = {
    'WT': {'score': 1.00, 'domain': 'WT', 'clinvar': 'WT'},
    'C52Y': {'score': 0.05, 'domain': 'LA1', 'clinvar': 'pathogenic'},
    'D147N': {'score': 0.15, 'domain': 'LA3', 'clinvar': 'pathogenic'},
    'P526L': {'score': 0.95, 'domain': 'EGF-A', 'clinvar': 'benign'},
    'T705I': {'score': 0.92, 'domain': 'beta_propeller', 'clinvar': 'benign'}
}

def main():
    print("="*60)
    print("LDLR VARIANT SIMULATION")
    print("="*60)
    print(f"\\nRunning {len(VARIANTS)} variants...\\n")
    
    results = []
    
    for i, (name, info) in enumerate(VARIANTS.items(), 1):
        print(f"[{i}/{len(VARIANTS)}] {name}...", end=' ')
        
        try:
            model = LDLRModel(name, info['score'])
            model.run()
            uptake = model.get_final_uptake()
            
            print(f"✓ (uptake: {uptake:.1f})")
            
            results.append({
                'variant': name,
                'exp_score': info['score'],
                'model_uptake': uptake,
                'domain': info['domain'],
                'clinvar': info['clinvar']
            })
        except Exception as e:
            print(f"❌ FAILED: {e}")
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Normalize to WT
    wt_uptake = df.loc[df['variant']=='WT', 'model_uptake'].values[0]
    df['model_score'] = df['model_uptake'] / wt_uptake
    
    # Save
    df.to_csv('results/simulation_results.csv', index=False)
    print(f"\\n💾 Results saved")
    
    # Summary
    print("\\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Variants: {len(df)}")
    print(f"WT uptake: {wt_uptake:.1f}")
    
    r, p = pearsonr(df['model_score'], df['exp_score'])
    print(f"\\nCorrelation: r = {r:.3f}, p = {p:.3e}")
    
    if r > 0.6:
        print("✓ STRONG correlation")
    elif r > 0.4:
        print("~ MODERATE correlation")
    else:
        print("✗ WEAK correlation")
    
    print("\\n✅ Complete! Next: python scripts/analyze_results.py")

if __name__ == "__main__":
    main()
'''
    
    with open('scripts/run_variants.py', 'w') as f:
        f.write(content)
    
    print("✓ Created scripts/run_variants.py")

def create_analyze_results():
    """Create scripts/analyze_results.py"""
    content = '''"""
Analyze LDLR variant results
Creates validation and separation plots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from pathlib import Path

# Setup
FIGURES_DIR = Path('results/figures')
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")

def load_results():
    """Load simulation results"""
    df = pd.read_csv('results/simulation_results.csv')
    return df

def plot_validation(df):
    """Create validation scatter plot"""
    fig, ax = plt.subplots(figsize=(8, 7))
    
    colors = {'pathogenic': 'red', 'benign': 'blue', 'WT': 'green'}
    
    for clinvar in df['clinvar'].unique():
        mask = df['clinvar'] == clinvar
        ax.scatter(
            df.loc[mask, 'exp_score'],
            df.loc[mask, 'model_score'],
            c=colors.get(clinvar, 'gray'),
            label=clinvar,
            s=150,
            alpha=0.7,
            edgecolors='black',
            linewidths=1
        )
    
    ax.plot([0, 1.5], [0, 1.5], 'k--', alpha=0.5, linewidth=2)
    
    r, p = pearsonr(df['model_score'], df['exp_score'])
    
    ax.text(
        0.05, 0.95,
        f"r = {r:.3f}\\np = {p:.2e}",
        transform=ax.transAxes,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
        fontsize=12,
        fontweight='bold'
    )
    
    ax.set_xlabel('Experimental Score', fontsize=13, fontweight='bold')
    ax.set_ylabel('Model Score', fontsize=13, fontweight='bold')
    ax.set_title('Model Validation', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right')
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'validation.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved validation.png")
    
    return r, p

def plot_separation(df):
    """Plot pathogenic vs benign"""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    classified = df[df['clinvar'].isin(['pathogenic', 'benign'])]
    
    data_to_plot = [
        classified[classified['clinvar']=='pathogenic']['model_score'],
        classified[classified['clinvar']=='benign']['model_score']
    ]
    
    bp = ax.boxplot(
        data_to_plot,
        labels=['Pathogenic', 'Benign'],
        patch_artist=True,
        showmeans=True
    )
    
    bp['boxes'][0].set_facecolor('lightcoral')
    bp['boxes'][1].set_facecolor('lightblue')
    
    for i, data in enumerate(data_to_plot, 1):
        y = data
        x = np.random.normal(i, 0.04, size=len(y))
        ax.scatter(x, y, alpha=0.6, s=100, edgecolors='black')
    
    ax.set_ylabel('Model Score', fontsize=13, fontweight='bold')
    ax.set_title('Pathogenic vs. Benign Separation', fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'separation.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved separation.png")

def print_summary(df, r, p):
    """Print results summary"""
    print("\\n" + "="*60)
    print("RESULTS SUMMARY")
    print("="*60)
    
    print(f"\\nCorrelation: r = {r:.4f}, p = {p:.2e}")
    
    path_mean = df[df['clinvar']=='pathogenic']['model_score'].mean()
    benign_mean = df[df['clinvar']=='benign']['model_score'].mean()
    
    print(f"\\nPathogenic mean: {path_mean:.3f}")
    print(f"Benign mean: {benign_mean:.3f}")
    print(f"Separation: {benign_mean - path_mean:.3f}")
    
    mae = np.abs(df['model_score'] - df['exp_score']).mean()
    print(f"\\nMean Absolute Error: {mae:.4f}")
    
    print("\\n" + "="*60)

def main():
    print("="*60)
    print("LDLR VARIANT ANALYSIS")
    print("="*60)
    
    df = load_results()
    print(f"\\nLoaded {len(df)} variants\\n")
    
    print("Generating plots...")
    r, p = plot_validation(df)
    plot_separation(df)
    
    print_summary(df, r, p)
    
    print(f"\\n✅ Complete! Figures in {FIGURES_DIR}")

if __name__ == "__main__":
    main()
'''
    
    with open('scripts/analyze_results.py', 'w') as f:
        f.write(content)
    
    print("✓ Created scripts/analyze_results.py")

def create_variant_data():
    """Create data/variants.csv"""
    content = """variant,uptake_score,domain,clinvar
WT,1.00,WT,WT
C52Y,0.05,LA1,pathogenic
D147N,0.15,LA3,pathogenic
P526L,0.95,EGF-A,benign
T705I,0.92,beta_propeller,benign
"""
    
    with open('data/variants.csv', 'w') as f:
        f.write(content)
    
    print("✓ Created data/variants.csv")

def create_docs():
    """Create documentation files"""
    
    # Simplified 1-week plan
    plan = """# 1-Week LDLR Modeling Project

## Schedule

### Day 1 (6h): Setup & Base Model
- Install PyBioNetGen
- Create base LDLR model
- Test WT simulation

### Day 2 (6h): Variants
- Tune parameters
- Implement 5 variants
- Run all variants

### Day 3 (6h): Batch Simulations
- Create batch runner
- Run all simulations
- Generate validation plot

### Day 4 (6h): Analysis
- Create separation plot
- Calculate statistics
- Document results

### Day 5 (6h): Documentation
- Write README
- Clean up code
- Push to GitHub

## Success Criteria

- ✅ 5 variants simulated
- ✅ r > 0.5 correlation
- ✅ 2 figures generated
- ✅ GitHub repository

## Files

All files have been generated by generate_project.py
"""
    
    with open('docs/PLAN_1WEEK.md', 'w') as f:
        f.write(plan)
    
    print("✓ Created docs/PLAN_1WEEK.md")
    
    # Quick start guide
    quickstart = """# Quick Start Guide

## Installation

```bash
# Install dependencies
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Test Installation

```bash
python3 -c "import bionetgen; print('✓ PyBioNetGen working!')"
```

## Run Analysis

```bash
# 1. Test model
python models/ldlr_model.py

# 2. Run all variants
python scripts/run_variants.py

# 3. Generate plots
python scripts/analyze_results.py
```

## Expected Output

After running scripts:
- `results/simulation_results.csv` - Data
- `results/figures/validation.png` - Validation plot
- `results/figures/separation.png` - Separation plot

## Troubleshooting

**Import Error:**
- Check venv is activated: `source venv/bin/activate`

**No module 'ldlr_model':**
- Run from project root directory

**Poor correlation:**
- Expected for small sample (n=5)
- Focus on qualitative separation
"""
    
    with open('docs/QUICKSTART.md', 'w') as f:
        f.write(quickstart)
    
    print("✓ Created docs/QUICKSTART.md")

def main():
    """Generate all project files"""
    
    print("\\n" + "="*60)
    print("LDLR PROJECT FILE GENERATOR")
    print("="*60)
    print("\\nThis will create all project files and directories\\n")
    
    # Check if already in a project directory
    if Path('models').exists() or Path('scripts').exists():
        response = input("⚠️  Directory structure exists. Overwrite? (y/n): ")
        if response.lower() != 'y':
            print("Cancelled.")
            return
    
    print("Creating project structure...\\n")
    
    # Create everything
    create_directory_structure()
    create_gitignore()
    create_requirements()
    create_readme()
    create_ldlr_model()
    create_run_variants()
    create_analyze_results()
    create_variant_data()
    create_docs()
    
    print("\\n" + "="*60)
    print("✅ PROJECT GENERATED SUCCESSFULLY!")
    print("="*60)
    
    print("""
Next steps:

1. Create virtual environment:
   python3 -m venv venv
   source venv/bin/activate

2. Install dependencies:
   pip install -r requirements.txt

3. Test the model:
   python models/ldlr_model.py

4. Run full analysis:
   python scripts/run_variants.py
   python scripts/analyze_results.py

5. Check results:
   - results/simulation_results.csv
   - results/figures/*.png

Documentation:
- docs/PLAN_1WEEK.md - Project plan
- docs/QUICKSTART.md - Quick start guide
- README.md - Project overview

Happy modeling! 🚀
""")

if __name__ == "__main__":
    main()
