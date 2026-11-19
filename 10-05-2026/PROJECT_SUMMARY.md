# Project summary: modeling LDLR variant effects

**Tianyi Zhang** · Faeder Lab, Department of Computational and Systems Biology, University of Pittsburgh
**Advisor:** James Faeder · **Collaboration:** Roth Lab (Fritz Roth, Cindy Zhang, Warren)

## Goal

Understand how missense variants in the LDL receptor (LDLR) disrupt LDL uptake, by linking deep mutational scanning data (Tabet et al., *Science* 2025: cell-surface abundance **A** and LDL uptake **F** for ~14,600 missense variants) to the steps of a mechanistic model of the LDLR cycle. The guiding question is which minimal parameter changes turn the wild-type phenotype into each mutant phenotype.

## What was done

### 1. Mechanistic model of LDL uptake
- Built a rule-based **BioNetGen** (ODE) model of the LDLR cycle: LDL binding at the surface through four ligand-binding modules, endocytosis, pH-dependent endosomal release, receptor recycling, and lysosomal LDL degradation.
- Variant scores map onto model parameters (A sets receptor supply; F sets uptake efficiency), and simulated surface LDLR and uptake are compared with measured A and F. With A = F = 1 the model reproduces the wild-type steady state.
- **Status:** working model (`ldlr_model_a9.py`).

### 2. Phenotype structure and parameter-space exploration
- **k-means clustering** (k = 4) of A–F space into WT-like, severe/null, abundance-defect and function-defect variants; **PCA** of the phenotype space.
- **Latin Hypercube Sampling** of parameter space (1,000 parameter sets) to map which phenotypes the model can reach; **random forest surrogate models** of model outputs.
- **Global sensitivity analysis** (log–log regression on LHS samples): receptor supply (`LDLR_init`) and endocytosis (`k_endo`) dominate; endosomal release (`k_off_endo`) is the least sensitive parameter.
- **Inverse fitting** with **PyPESTO** (multi-start optimization): many phenotypes can be fit, but very low abundance with near-normal function is unreachable, and different parameter sets fit equally well, so A and F alone cannot identify the mechanism.
- **Status:** complete.

### 3. Grounding the model in variant biology
- Analyzed phenotype classes by mutation chemistry and by LDLR region. Region matters more than mutation type: LA3–LA5 are strongly enriched for function defects.
- Derived three domain-to-mechanism hypotheses: LA3–5 → surface binding (`k_off_surf`), cytoplasmic NPxY motif → endocytosis (`k_endo`), β-propeller → endosomal release (`k_off_endo`).
- **Status:** hypotheses defined; a domain-aware model that applies them (`ldlr_model_a13.py`) is in development.

### 4. Bayesian parameter inference
- Fit a hierarchical Bayesian model (**PyMC**, NUTS sampling) mapping variant features to model parameters through a linear forward model.
- **Finding:** parameters learned under the linear approximation were about 9× too large to reproduce the observed phenotypes in the nonlinear mass-action model, so linear-surrogate inference does not transfer to the ODE model.
- **Status:** complete; this result motivated the feature-based feasibility test below.

### 5. Feature-based prediction of phenotype components
- Defined two targets: **S = A** (surface abundance) and **B = F − LOWESS(F | A)** (function beyond what abundance explains).
- Built a 26-feature design: substitution chemistry (BLOSUM62, hydrophobicity, volume, polarity, charge, Cys/Pro/Gly), domain identity, and structural context from **PDB 9BDE** (LDLR–ApoB100 complex) and **PDB 1N7D** (LDLR at endosomal pH). MutateX ΔΔG is part of the design but awaits values from the Roth lab (25 features used).
- **Ridge regression** with domain × feature interactions and nested 5-fold cross-validation, on 14,617 missense variants.
- **Results:**
  - Global out-of-fold r = 0.42 for B and 0.49 for S; domain interactions improve both (from 0.35 and 0.40).
  - Prediction is domain-specific: B is best predicted in LA3–LA5 (r up to 0.57), S in EGF-A, EGF-C, LA7 and the signal peptide (r up to 0.65).
  - Local predictability tracks within-domain variance (r = 0.77 for B, 0.74 for S across 17 domains).
  - The **β-propeller** is a genuine failure region: real variation but r = 0.14 (B) and 0.29 (S); 28% of variants but 46% of the worst-predicted 5%. Adding its structure from 1N7D did not help, which points to a missing mechanism (pH-dependent release) rather than missing structure.
- **Status:** complete; presented at the joint Faeder–Roth lab meeting (October 2026).

### 6. Clinical classification check (initial)
- Compared assay scores and model predictions with ClinGen VCEP and ClinVar benign/pathogenic labels (Tabet Data S4–S5).
- **Findings:** the measured uptake score classifies best (AUC 0.87), ahead of feature-based predictions (about 0.67). 55 of 93 VCEP-pathogenic missense variants look wild-type-like in the assay, concentrated in the β-propeller and LA2/LA6. These are candidates for mechanisms the cell assay does not capture.
- **Status:** initial analysis done; to be continued with the PCSK9 model.

### 7. Data quality and reproducibility
- Removed an outcome-derived label that leaked information into the features, and 774 stop-gain variants hidden in the missense set.
- Verified every model input against the Tabet supplementary tables and the PDB structures (structural features computed from the LDLR chain only, with residue-identity checks).
- Released a documented, reproducible code base: mechanistic models, B-score computation, structural feature extraction and the regression pipeline (see `README.md`).

## Current status

| Component | Status |
|---|---|
| Mechanistic model (A9) | Working |
| Parameter exploration, sensitivity, inverse fitting | Complete |
| Bayesian parameter inference | Complete (negative result informed next steps) |
| Feature-based prediction | Complete, presented October 2026 |
| Domain-aware model (A13) | In development |
| Clinical classification check | Initial analysis |
| PCSK9 extension | Next step |

## Next steps

Work will resume after the PhD application period.

1. **Add PCSK9 to the model** (agreed with the Faeder and Roth labs). PCSK9 binds LDLR's EGF-A domain and diverts the receptor from recycling to degradation. This requires adding receptor synthesis and degradation, PCSK9 binding at the surface and in the endosome, and PCSK9-driven degradation, then mapping EGF-A variants to PCSK9 affinity.
2. **Test whether the PCSK9-aware model explains the "silent" pathogenic variants** that look normal in the assay.
3. **Test the β-propeller release mechanism** (`k_off_endo`, alone or with endocytosis and recycling) in the model.
4. **Complete the structural features:** MutateX ΔΔG (requested from the Roth lab) and propeller–LA4/5 distances from 1N7D.

## Methods and tools

BioNetGen (rule-based ODE modeling), PyPESTO, PyMC, Latin Hypercube Sampling, global sensitivity analysis, k-means, PCA, random forest, ridge regression with cross-validation, LOWESS, Biopython, DSSP, PyMOL; Python (numpy, pandas, scikit-learn, statsmodels, matplotlib).

## Materials

- `README.md`: every file in this repository and how to reproduce the results
- `LDLR_presentation_Oct2026.pdf`: slides from the joint Faeder–Roth lab meeting
