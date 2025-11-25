# LDLR Functional Landscape Modeling

Rule-based computational model of LDLR for understanding familial hypercholesterolemia.

## Quick Start

```bash
# Setup (first time)
./setup_project_pybngl.sh

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
