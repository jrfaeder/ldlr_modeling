# LDLR Modeling Project: 1-Week Sprint

## Project Goals (Streamlined)

### Primary Objective
Build a functional LDLR model using PyBioNetGen that demonstrates:
- LDL uptake by wild-type LDLR
- Functional differences between pathogenic and benign variants
- Basic validation against experimental data

### Success Criteria
1. ✅ Working model with realistic uptake kinetics
2. ✅ 5 variants simulated (WT + 2 pathogenic + 2 benign)
3. ✅ Correlation r > 0.5 with experimental scores
4. ✅ Clear separation between pathogenic and benign
5. ✅ 2 publication-quality figures
6. ✅ Documented GitHub repository

---

## Daily Schedule (5 days × 6 hours = 30 hours)

### **Day 1: Setup & Basic Model** (Monday)

#### Afternoon Session (3 hours)
- [ ] Build base LDLR model (simplified)
  - 4 LA modules (LA3, LA4, LA5, LA7)
  - Basic LDL binding
  - Simple trafficking
- [ ] Run simulation
- [ ] Plot basic dynamics

**Deliverable**: Working WT LDLR model showing LDL uptake

---

### **Day 2: Variants & Calibration** (Tuesday)

#### Morning Session (3 hours)
- [ ] Quick parameter tuning
  - Adjust uptake timescale
  - Reasonable receptor numbers
- [ ] Test LA3 knockout (confirm it's critical)

#### Afternoon Session (3 hours)
- [ ] Define 5 variants with functional scores
- [ ] Create variant generator function
- [ ] Run all 5 variants manually
- [ ] Verify outputs look reasonable

**Deliverable**: 5 variants simulated successfully

---

### **Day 3: Batch Automation** (Wednesday)

#### Morning Session (3 hours)
- [ ] Create batch simulation script
- [ ] Parse all outputs into DataFrame
- [ ] Calculate model scores (normalized to WT)

#### Afternoon Session (3 hours)
- [ ] Create validation scatter plot
- [ ] Calculate correlation with experimental data
- [ ] Identify any major outliers

**Deliverable**: Validation plot with r value

---

### **Day 4: Analysis & Visualization** (Thursday)

#### Morning Session (3 hours)
- [ ] Create pathogenic vs benign comparison plot
- [ ] Calculate separation statistics
- [ ] Test if they're distinguishable

#### Afternoon Session (3 hours)
- [ ] Create combined summary figure
- [ ] Write results summary
- [ ] Document key findings

**Deliverable**: 2 publication-ready figures + results text

---

### **Day 5: Documentation & Wrap-Up** (Friday)

#### Morning Session (3 hours)
- [ ] Write comprehensive README
  - Project overview
  - Results summary
  - Key findings
  - Limitations
- [ ] Clean up code
- [ ] Add comments to key sections

#### Afternoon Session (3 hours)
- [ ] Create requirements.txt
- [ ] Set up GitHub repository
- [ ] Push all code and results
- [ ] Write future directions

**Deliverable**: Complete shareable GitHub repository

---

## Simplified Model Details

### Model Components

```
LDLR molecule:
- 4 LA modules: la3, la4, la5, la7
- Location: surface or endosome
- No cooperative binding factors

LDL molecule:
- Single binding site (simplified)
- Location: extracellular or internalized

Processes:
1. LDLR-LDL binding (4 independent sites)
2. Endocytosis
3. Receptor recycling
4. LDL degradation
```

### 5 Variants to Test

| Variant | Domain | Uptake Score | ClinVar |
|---------|--------|--------------|---------|
| WT | - | 1.00 | WT |
| C52Y | LA1 | 0.05 | Pathogenic |
| D147N | LA3 | 0.15 | Pathogenic |
| P526L | EGF-A | 0.95 | Benign |
| T705I | β-propeller | 0.92 | Benign |

---

## Project Structure

```
ldlr_modeling/
├── README.md
├── requirements.txt
├── setup.py (optional)
├── models/
│   └── ldlr_model.py          # Model as Python class
├── scripts/
│   ├── run_variants.py        # Batch simulation
│   └── analyze_results.py     # Plotting & stats
├── data/
│   └── variants.csv           # Variant functional scores
├── results/
│   ├── figures/
│   │   ├── validation.png
│   │   └── separation.png
│   └── simulation_data.csv
└── docs/
    ├── PLAN_1WEEK.md
    └── RESULTS.md
```

---

## Key Files to Create

### Day 1
1. `models/ldlr_model.py` - Base model class
2. `scripts/test_model.py` - Test script

### Day 2
3. `scripts/run_variants.py` - Batch runner
4. `data/variants.csv` - Variant data

### Day 3-4
5. `scripts/analyze_results.py` - Analysis
6. `results/figures/*.png` - Plots

### Day 5
7. `README.md` - Documentation
8. `requirements.txt` - Dependencies

---

## Technology Stack

### Required
- Python 3.8+
- PyBioNetGen (`pip install pybiogen`)
- pandas, numpy, scipy
- matplotlib, seaborn

### Optional
- Jupyter (for interactive development)
- pytest (for testing)

---

## Success Metrics (Reduced Bar)

### Minimum Viable
- ✅ r ≥ 0.4 (weak but positive correlation)
- ✅ Pathogenic mean < 0.3
- ✅ Benign mean > 0.7
- ✅ Working code on GitHub

### Target
- ✅ r ≥ 0.5 (moderate correlation)
- ✅ Clear separation (benign - pathogenic > 0.5)
- ✅ Model explains >50% of variance
- ✅ Complete documentation

### Stretch
- ✅ r ≥ 0.6 (good correlation)
- ✅ All pathogenic < 0.4, all benign > 0.8
- ✅ Mechanistic insight documented
- ✅ Tutorial notebook

---

## 3-Day MVP Backup Plan

If you're behind schedule, here's the absolute minimum:

### Day 1: Model + WT only
- Setup PyBioNetGen
- Create working WT model
- One plot showing uptake

### Day 2: Add 3 variants
- WT + C52Y (pathogenic) + P526L (benign)
- Run all 3
- Simple comparison

### Day 3: Wrap up
- One validation plot
- Calculate r value
- Basic README
- GitHub push

**MVP deliverable**: Proof that model can distinguish pathogenic from benign (n=3)

---

## Daily Time Commitment

| Day | Hours | Can't Miss | Optional |
|-----|-------|------------|----------|
| Mon | 6 | Setup, base model | - |
| Tue | 6 | Variants | LA3 knockout |
| Wed | 6 | Batch script, validation plot | - |
| Thu | 6 | Separation plot | Combined figure |
| Fri | 6 | README | GitHub setup |

**Total: 30 hours**

---

## Risk Mitigation

### If PyBioNetGen installation fails
- **Backup**: Use BioNetGen CLI with subprocess calls
- **Time**: Add 2 hours for troubleshooting

### If model doesn't work by Day 2
- **Backup**: Use even simpler model (2 LA modules)
- **Cut**: Module knockout analysis

### If correlation is poor (r < 0.4)
- **Backup**: Focus on qualitative separation
- **Document**: Explain why (model limitations)

### If running behind
- **Cut first**: Mechanistic analysis, fancy plots
- **Cut second**: Documentation details
- **Never cut**: Basic validation, GitHub repo

---

## Post-Week Extensions

After completing the 1-week sprint, you can:

**Week 2** (if time):
- Add remaining 5-10 variants
- Implement cooperative binding
- Add module-specific analysis
- Improve parameterization

**Week 3+** (future work):
- VLDL binding
- Clinical risk prediction
- Spatial modeling
- Publication-quality manuscript

---

## Checklist

### Pre-Week Preparation
- [ ] Python 3.8+ installed
- [ ] Project folder created
- [ ] GitHub account ready
- [ ] Paper PDF downloaded

### End-of-Week Deliverables
- [ ] Working LDLR model (code)
- [ ] 5 variants simulated (data)
- [ ] 2 figures (validation + separation)
- [ ] Results summary (text)
- [ ] GitHub repository (public)
- [ ] README with findings

### Bonus (if time)
- [ ] Module knockout analysis
- [ ] Parameter documentation
- [ ] Tutorial notebook
- [ ] Future work section

---

## Daily Progress Tracking

### Monday: ___/___/___
**Completed**: 
**Blockers**:
**Tomorrow**:

### Tuesday: ___/___/___
**Completed**:
**Blockers**:
**Tomorrow**:

### Wednesday: ___/___/___
**Completed**:
**Blockers**:
**Tomorrow**:

### Thursday: ___/___/___
**Completed**:
**Blockers**:
**Tomorrow**:

### Friday: ___/___/___
**Completed**:
**Blockers**:
**Next steps**:

---

## Final Notes

This 1-week version is:
- ✅ **Achievable** with focused effort
- ✅ **Publishable** as proof-of-concept
- ✅ **Expandable** to full 2-week version
- ✅ **Defendable** with clear limitations documented

**Good luck! 🚀**
