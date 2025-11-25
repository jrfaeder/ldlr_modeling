# Quick Start Guide

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
