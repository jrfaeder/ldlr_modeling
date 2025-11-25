# LDLR Functional Landscape Modeling

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
