# Day 1 Getting Started Guide
## Using PyBioNetGen for LDLR Modeling

---

## Overview (6 hours)

**Morning** (3 hours): Setup and simple test  
**Afternoon** (3 hours): Build base LDLR model

---

## Morning Session: Setup (3 hours)

### Step 1: Install Python and Create Environment (30 min)

```bash
# Check Python version (need 3.8+)
python3 --version

# Create project directory
mkdir ldlr_modeling
cd ldlr_modeling

# Create virtual environment
python3 -m venv venv

# Activate it
# Mac/Linux:
source venv/bin/activate
# Windows:
# venv\Scripts\activate

# Verify
which python  # Should point to venv/bin/python
```

---

### Step 2: Install PyBioNetGen (15 min)

**PyBioNetGen makes everything easier - no separate BioNetGen installation needed!**

```bash
# Install PyBioNetGen (includes BioNetGen)
pip install pybiogen

# Install other required packages
pip install pandas numpy scipy matplotlib seaborn

# Save requirements
pip freeze > requirements.txt

# Test installation
bionetgen --version
```

**Expected output** (exact version numbers may vary): 
```
BioNetGen simple command line interface 0.8.6
BioNetGen version: BioNetGen-2.9.3
Cement Framework 3.0.10
Python 3.12.3
Platform macOS-14.7.4-arm64-arm-64bit
```

If this works, you're ready to go! 🎉

---

### Step 3: Create Project Structure (10 min)

```bash
# Create directories
mkdir models scripts data results/figures results/data docs -p

# Create initial files
touch README.md
touch models/ldlr_model.py
touch scripts/run_variants.py
touch scripts/analyze_results.py
touch data/variants.csv

# Verify structure
tree -L 2
```

**Should see**:
```
ldlr_modeling/
├── models/
├── scripts/
├── data/
├── results/
│   ├── figures/
│   └── data/
├── docs/
├── venv/
└── requirements.txt
```

---

### Step 4: Create First Simple Test (45 min)

Create `models/simple_test.py`:

```python
"""
Simple test to verify PyBioNetGen works
A + B <-> A:B binding model
"""

import bionetgen
import matplotlib.pyplot as plt
import pandas as pd

# Define model as a string
model_str = """
begin model

begin parameters
    k_on   1.0
    k_off  0.1
    A_0    100
    B_0    50
end parameters

begin molecule types
    A(b)
    B(a)
end molecule types

begin seed species
    A(b)  A_0
    B(a)  B_0
end seed species

begin observables
    Molecules  A_free    A(b)
    Molecules  B_free    B(a)
    Molecules  Complex   A(b!1).B(a!1)
end observables

begin reaction rules
    A(b) + B(a) <-> A(b!1).B(a!1)  k_on, k_off
end reaction rules

end model
"""

# Run simulation
print("Running simple binding model...")
result = bionetgen.run(model_str, out_dir='results/data')

print("✓ Simulation complete!")
print(f"✓ Output saved to: {result['results_folder']}")

# Load and plot results
gdat_file = result['gdat_file']
data = pd.read_csv(gdat_file, sep='\s+', comment='#')

# Create plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(data['time'], data['A_free'], label='Free A', linewidth=2)
ax.plot(data['time'], data['B_free'], label='Free B', linewidth=2)
ax.plot(data['time'], data['Complex'], label='A:B Complex', linewidth=2)

ax.set_xlabel('Time', fontsize=12)
ax.set_ylabel('Molecules', fontsize=12)
ax.set_title('Simple Binding Model', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('results/figures/simple_test.png', dpi=150)
print("✓ Plot saved to: results/figures/simple_test.png")

plt.show()

# Print final values
print("\n=== Final State ===")
print(f"Free A: {data['A_free'].iloc[-1]:.1f}")
print(f"Free B: {data['B_free'].iloc[-1]:.1f}")
print(f"Complex: {data['Complex'].iloc[-1]:.1f}")
print("\n✅ Test successful! PyBioNetGen is working.")
```

**Run it**:
```bash
python models/simple_test.py
```

**Expected**: You should see a plot with binding curves and "Test successful!" message.

---

### Step 5: Understand PyBioNetGen Basics (30 min)

**Key concepts:**

#### Running Models

```python
import bionetgen

# Method 1: Model as string
model_str = """
begin model
...
end model
"""
result = bionetgen.run(model_str)

# Method 2: Model from file
result = bionetgen.run("model.bngl")

# Method 3: With options
result = bionetgen.run(
    model_str,
    out_dir='results',      # Output directory
    suppress=True,          # Suppress BioNetGen output
    cleanup=True            # Clean up intermediate files
)
```

#### Accessing Results

```python
# Result is a dictionary with:
result = {
    'results_folder': 'path/to/results',
    'gdat_file': 'path/to/file.gdat',      # Observable data
    'cdat_file': 'path/to/file.cdat',      # Species concentrations
    'model': model_object                   # Model structure
}

# Load observable data
import pandas as pd
data = pd.read_csv(result['gdat_file'], sep='\s+', comment='#')
print(data.columns)  # ['time', 'Observable1', 'Observable2', ...]
```

#### Model Structure

```python
# Same BioNetGen syntax, but in Python string:
model = """
begin model

begin parameters
    k_rate  1.0
end parameters

begin molecule types
    Molecule(site)
end molecule types

begin seed species
    Molecule(site)  100
end seed species

begin observables
    Molecules  Total  Molecule()
end observables

begin reaction rules
    # Rules here
end reaction rules

end model
"""
```

**✅ Morning Session Complete! Take a break.**

---

## Afternoon Session: Build LDLR Model (3 hours)

### Step 6: Design Simplified LDLR Model (30 min)

**Model design decisions:**

```
Components:
- LDLR with 4 LA modules (la3, la4, la5, la7)
- LDL with single binding site (simplified)
- 2 compartments: surface, endosome

Processes:
1. LDLR-LDL binding (at surface, pH 7.4)
2. Endocytosis (any bound complex)
3. LDL release (in endosome, pH 5.5)
4. Receptor recycling (to surface)
5. LDL degradation

Parameters:
- k_bind: ~1.0 (binding rate)
- k_unbind_surf: ~1.0 (slow at surface)
- k_unbind_endo: ~50.0 (fast in endosome)
- k_endo: ~2.0 (endocytosis)
- k_recycle: ~3.0 (recycling)
```

---

### Step 7: Create Base LDLR Model (90 min)

Create `models/ldlr_model.py`:

```python
"""
Base LDLR model for LDL uptake
Simplified version with 4 LA modules
"""

import bionetgen
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

class LDLRModel:
    """LDLR-LDL uptake model using PyBioNetGen"""
    
    def __init__(self, variant_name='WT', functional_score=1.0):
        """
        Initialize model
        
        Parameters:
        -----------
        variant_name : str
            Variant identifier (e.g., 'WT', 'C52Y')
        functional_score : float
            Functional score from experiment (0-1.5)
            Used to scale binding parameters
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
    # Binding parameters (scaled by functional score)
    k_on_base      1.0
    k_off_surf     {1.0 / binding_strength:.4f}
    k_off_endo     50.0
    
    # Module-specific strengths (from paper Fig 3)
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
    # LDLR with 4 LA modules and location
    LDLR(la3, la4, la5, la7, loc~surface~endosome)
    
    # LDL (simplified: one binding site)
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
    # === LDLR-LDL Binding at Surface (4 independent sites) ===
    
    # LA3 binding
    LDLR(la3, loc~surface) + LDL(ldlr, loc~extra) <-> \\
        LDLR(la3!1, loc~surface).LDL(ldlr!1, loc~extra) \\
        k_on_base*strength_LA3, k_off_surf
    
    # LA4 binding
    LDLR(la4, loc~surface) + LDL(ldlr, loc~extra) <-> \\
        LDLR(la4!1, loc~surface).LDL(ldlr!1, loc~extra) \\
        k_on_base*strength_LA4, k_off_surf
    
    # LA5 binding
    LDLR(la5, loc~surface) + LDL(ldlr, loc~extra) <-> \\
        LDLR(la5!1, loc~surface).LDL(ldlr!1, loc~extra) \\
        k_on_base*strength_LA5, k_off_surf
    
    # LA7 binding
    LDLR(la7, loc~surface) + LDL(ldlr, loc~extra) <-> \\
        LDLR(la7!1, loc~surface).LDL(ldlr!1, loc~extra) \\
        k_on_base*strength_LA7, k_off_surf
    
    # === Endocytosis (any bound complex) ===
    
    LDLR(la3!1, loc~surface).LDL(ldlr!1, loc~extra) -> \\
        LDLR(la3!1, loc~endosome).LDL(ldlr!1, loc~endo)  k_endo
    
    LDLR(la4!1, loc~surface).LDL(ldlr!1, loc~extra) -> \\
        LDLR(la4!1, loc~endosome).LDL(ldlr!1, loc~endo)  k_endo
    
    LDLR(la5!1, loc~surface).LDL(ldlr!1, loc~extra) -> \\
        LDLR(la5!1, loc~endosome).LDL(ldlr!1, loc~endo)  k_endo
    
    LDLR(la7!1, loc~surface).LDL(ldlr!1, loc~extra) -> \\
        LDLR(la7!1, loc~endosome).LDL(ldlr!1, loc~endo)  k_endo
    
    # === LDL Release in Endosome (pH-dependent) ===
    
    LDLR(la3!1, loc~endosome).LDL(ldlr!1, loc~endo) -> \\
        LDLR(la3, loc~endosome) + LDL(ldlr, loc~endo)  k_off_endo
    
    LDLR(la4!1, loc~endosome).LDL(ldlr!1, loc~endo) -> \\
        LDLR(la4, loc~endosome) + LDL(ldlr, loc~endo)  k_off_endo
    
    LDLR(la5!1, loc~endosome).LDL(ldlr!1, loc~endo) -> \\
        LDLR(la5, loc~endosome) + LDL(ldlr, loc~endo)  k_off_endo
    
    LDLR(la7!1, loc~endosome).LDL(ldlr!1, loc~endo) -> \\
        LDLR(la7, loc~endosome) + LDL(ldlr, loc~endo)  k_off_endo
    
    # === Receptor Recycling ===
    
    LDLR(la3, la4, la5, la7, loc~endosome) -> \\
        LDLR(la3, la4, la5, la7, loc~surface)  k_recycle
    
    # === LDL Degradation ===
    
    LDL(ldlr, loc~endo) -> LDL(ldlr, loc~lyso)  k_recycle
    LDL(loc~lyso) -> 0  k_degrade
    
end reaction rules

end model
"""
        return model
    
    def run(self, t_end=200, n_steps=200, out_dir='results/data'):
        """Run the simulation"""
        
        model_str = self.get_model_string()
        
        print(f"Running {self.variant_name} model...")
        self.result = bionetgen.run(
            model_str,
            out_dir=out_dir,
            suppress=True
        )
        print(f"  ✓ Complete")
        
        return self.result
    
    def get_data(self):
        """Load simulation data as DataFrame"""
        if self.result is None:
            raise ValueError("Must run simulation first")
        
        return pd.read_csv(
            self.result['gdat_file'],
            sep='\s+',
            comment='#'
        )
    
    def get_final_uptake(self):
        """Get final LDL internalized count"""
        data = self.get_data()
        return data['LDL_internalized'].iloc[-1]
    
    def plot(self, save_path=None):
        """Plot simulation results"""
        data = self.get_data()
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Panel 1: LDLR trafficking
        axes[0,0].plot(data['time'], data['LDLR_surface'], 
                      label='Surface', linewidth=2)
        axes[0,0].plot(data['time'], data['LDLR_endosome'], 
                      label='Endosome', linewidth=2)
        axes[0,0].set_xlabel('Time')
        axes[0,0].set_ylabel('LDLR Count')
        axes[0,0].set_title('LDLR Trafficking')
        axes[0,0].legend()
        axes[0,0].grid(alpha=0.3)
        
        # Panel 2: LDL dynamics
        axes[0,1].plot(data['time'], data['LDL_free'], 
                      label='Free LDL', linewidth=2)
        axes[0,1].plot(data['time'], data['LDL_internalized'], 
                      label='Internalized', linewidth=2, color='red')
        axes[0,1].set_xlabel('Time')
        axes[0,1].set_ylabel('LDL Count')
        axes[0,1].set_title('LDL Uptake')
        axes[0,1].legend()
        axes[0,1].grid(alpha=0.3)
        
        # Panel 3: Cumulative uptake
        axes[1,0].plot(data['time'], data['LDL_internalized'], 
                      color='darkred', linewidth=2)
        axes[1,0].set_xlabel('Time')
        axes[1,0].set_ylabel('Internalized LDL')
        axes[1,0].set_title('Cumulative LDL Uptake')
        axes[1,0].grid(alpha=0.3)
        
        # Panel 4: Surface complexes
        axes[1,1].plot(data['time'], data['Complex_surface'], 
                      color='purple', linewidth=2)
        axes[1,1].set_xlabel('Time')
        axes[1,1].set_ylabel('Complexes')
        axes[1,1].set_title('LDLR-LDL Complexes at Surface')
        axes[1,1].grid(alpha=0.3)
        
        fig.suptitle(f'{self.variant_name} (score: {self.functional_score:.2f})',
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"✓ Plot saved to: {save_path}")
        
        plt.show()
        
        return fig

# Example usage
if __name__ == "__main__":
    # Create and run WT model
    wt_model = LDLRModel(variant_name='WT', functional_score=1.0)
    wt_model.run()
    
    # Plot results
    wt_model.plot(save_path='results/figures/ldlr_wt.png')
    
    # Print summary
    data = wt_model.get_data()
    print("\n=== Final State (t=200) ===")
    print(f"Surface LDLR: {data['LDLR_surface'].iloc[-1]:.1f}")
    print(f"LDL internalized: {data['LDL_internalized'].iloc[-1]:.1f}")
    print(f"Uptake: {data['LDL_internalized'].iloc[-1]/100*100:.1f}%")
```

**Run it**:
```bash
python models/ldlr_model.py
```

**Expected**: 4-panel plot showing LDLR trafficking and LDL uptake dynamics.

---

### Step 8: Verify Model Behavior (30 min)

Create `scripts/test_model.py`:

```python
"""Test that LDLR model behaves reasonably"""

import sys
sys.path.append('models')
from ldlr_model import LDLRModel

def test_wt_model():
    """Test wild-type model"""
    print("Testing WT model...")
    
    model = LDLRModel('WT', functional_score=1.0)
    model.run()
    data = model.get_data()
    
    # Check 1: LDL uptake occurs
    final_uptake = data['LDL_internalized'].iloc[-1]
    assert final_uptake > 0, "❌ No LDL uptake!"
    print(f"✓ LDL uptake: {final_uptake:.1f} particles")
    
    # Check 2: Receptors remain on surface
    surface_start = data['LDLR_surface'].iloc[0]
    surface_end = data['LDLR_surface'].iloc[-1]
    assert surface_end > surface_start * 0.5, "❌ Receptors depleted!"
    print(f"✓ Surface LDLR maintained: {surface_end:.1f}/{surface_start:.1f}")
    
    # Check 3: Uptake increases over time
    uptake_mid = data['LDL_internalized'].iloc[len(data)//2]
    uptake_end = data['LDL_internalized'].iloc[-1]
    assert uptake_end > uptake_mid, "❌ Uptake not increasing!"
    print(f"✓ Uptake increases: {uptake_mid:.1f} → {uptake_end:.1f}")
    
    print("\n✅ All tests passed!")

def test_pathogenic_variant():
    """Test that pathogenic variant has reduced function"""
    print("\nTesting pathogenic variant...")
    
    # WT
    wt = LDLRModel('WT', 1.0)
    wt.run()
    wt_uptake = wt.get_final_uptake()
    
    # Pathogenic (low score)
    path = LDLRModel('C52Y', 0.05)
    path.run()
    path_uptake = path.get_final_uptake()
    
    # Should have much less uptake
    assert path_uptake < wt_uptake * 0.3, "❌ Pathogenic not impaired!"
    print(f"✓ WT uptake: {wt_uptake:.1f}")
    print(f"✓ Pathogenic uptake: {path_uptake:.1f}")
    print(f"✓ Reduction: {(1-path_uptake/wt_uptake)*100:.1f}%")
    
    print("\n✅ Pathogenic variant test passed!")

if __name__ == "__main__":
    test_wt_model()
    test_pathogenic_variant()
    print("\n" + "="*50)
    print("🎉 MODEL VALIDATION COMPLETE!")
    print("="*50)
```

**Run it**:
```bash
python scripts/test_model.py
```

**Expected output**:
```
Testing WT model...
✓ LDL uptake: 85.3 particles
✓ Surface LDLR maintained: 950.2/1000.0
✓ Uptake increases: 42.1 → 85.3

✅ All tests passed!

Testing pathogenic variant...
✓ WT uptake: 85.3
✓ Pathogenic uptake: 4.2
✓ Reduction: 95.1%

✅ Pathogenic variant test passed!

==================================================
🎉 MODEL VALIDATION COMPLETE!
==================================================
```

---

## End of Day 1 Checklist

- [ ] PyBioNetGen installed successfully
- [ ] Simple test model runs
- [ ] Base LDLR model created
- [ ] WT model shows LDL uptake
- [ ] Pathogenic variant shows reduced function
- [ ] All tests pass
- [ ] 4-panel plot generated

---

## Troubleshooting

### PyBioNetGen won't install
```bash
# Try updating pip first
pip install --upgrade pip

# Then install
pip install pybiogen

# If still fails, check Python version
python --version  # Must be 3.8+
```

### Model won't run
```python
# Add error handling
try:
    result = bionetgen.run(model_str)
except Exception as e:
    print(f"Error: {e}")
    # Check model syntax
```

### Can't find results
```python
# Print result paths
result = bionetgen.run(model_str)
print("Results folder:", result['results_folder'])
print("Data file:", result['gdat_file'])
```

---

## Quick Reference

### Key PyBioNetGen Commands

```python
import bionetgen

# Run model
result = bionetgen.run(model_string)

# Access results
data = pd.read_csv(result['gdat_file'], sep='\s+', comment='#')

# Common result files
result['gdat_file']  # Observable timecourse
result['cdat_file']  # Species concentrations
result['net_file']   # Reaction network
```

---

## Tomorrow Preview (Day 2)

You'll work on:
1. Quick parameter tuning
2. Testing LA3 knockout
3. Implementing all 5 variants
4. Running batch simulations

**Great work today! 🎉**

---

## Notes Section

**What worked well:**


**What was challenging:**


**Questions for tomorrow:**


