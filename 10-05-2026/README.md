# LDLR variant modeling

Mechanistic and statistical modeling of LDLR missense variant effects, built on the deep mutational scanning data of Tabet et al. (*Science* 2025), which measured LDL uptake (function, **F**) and cell-surface abundance (**A**) for nearly all possible LDLR amino-acid substitutions.

The work has three parts:

1. **Mechanistic model.** A rule-based BioNetGen model of the LDLR uptake cycle (surface binding, endocytosis, endosomal release, receptor recycling, LDL degradation), used to explore which parameter changes can reproduce mutant A/F phenotypes.
2. **Domain rules.** Phenotype clusters, mutation type and LDLR region are combined to link groups of variants to specific model parameters.
3. **Feature-based prediction.** A ridge regression tests whether variant features can predict two phenotype components directly:
   **S = A** (surface abundance) and **B = F − LOWESS(F | A)** (function beyond what abundance explains).

---

## Files

| File | Purpose |
|---|---|
| `ldlr_model_a9.py` | BioNetGen model of LDL uptake; simulates one variant from its A and F scores |
| `functional_abundance_clusters_k4.csv` | Input for A9: A and F scores with phenotype cluster |
| `ldlr_model_a13.py` | Development version of the model with domain-specific parameter mapping |
| `ldlr_variant_scores_v5.csv` | Main variant table: scores, B, annotations and structural features |
| `ldlr_compute_B.py` | Computes B from Tabet Data S1 and S2 |
| `ldlr_structure_features.py` | Computes structural features from PDB 9BDE and PDB 1N7D |
| `ldlr_ridge_v5.py` | Ridge regression predicting B and S; produces all regression results |
| `science_ady7186_data_s1.csv` | Tabet et al. Data S1: LDL uptake scores (F) |
| `science_ady7186_data_s2.csv` | Tabet et al. Data S2: cell-surface abundance scores (A) |
| `science_ady7186_data_s4.csv` | Tabet et al. Data S4: ClinGen FH VCEP variant classifications |

---

## Scripts

### `ldlr_model_a9.py`: mechanistic model
Builds and simulates (ODE) the BNGL model of the LDLR cycle for one variant. The variant's scores set the parameters:

- A scales receptor supply: `LDLR_init`, `k_recycle_ldlr`
- F scales uptake efficiency: `k_endo`, `k_off_surf = 1/F`
- All other rates stay at wild type

With A = F = 1 the model reproduces the wild-type steady state. Simulated surface LDLR and internalized LDL are reported relative to wild type and compared with the measured A and F. Any rate constant can be set directly from the command line, which is how the Latin Hypercube Sampling, sensitivity analysis and inverse fitting were run.

```bash
python ldlr_model_a9.py p.Pro3Ile --csv functional_abundance_clusters_k4.csv
python ldlr_model_a9.py p.Pro3Ile --k_endo 1.2 --k_off_surf 0.8   # direct parameter values
```

### `ldlr_model_a13.py`: domain-aware model (in development)
Uses the same reaction network as A9. Each variant changes only the parameter(s) linked to its phenotype cluster and LDLR region (the `route` column of v5). This version is not finished: its outputs have not been validated and are not used in the current results. It is the starting point for future steps, including the β-propeller endosomal-release (`k_off_endo`) mechanism.

```bash
python ldlr_model_a13.py p.Asp168Lys --csv ldlr_variant_scores_v5.csv
python ldlr_model_a13.py --batch --subset C3_LA3_5 --csv ldlr_variant_scores_v5.csv --max 20
```

### `ldlr_compute_B.py`: B score
Merges Data S1 (F) and S2 (A), fits LOWESS of F on A using the 14,617 true missense variants (`frac = 0.3`, 3 robustifying iterations), and evaluates the curve for every variant by linear interpolation. B is the residual: B < 0 means less uptake than expected for the variant's surface level. Reproduces `uptake_loess_predicted` and `B_score` in v5 exactly.

```bash
python ldlr_compute_B.py --s1 science_ady7186_data_s1.csv --s2 science_ady7186_data_s2.csv \
    --out B_scores.csv --check ldlr_variant_scores_v5.csv
```

### `ldlr_structure_features.py`: structural features
Computes, per LDLR residue:

- `dist_apob100`: minimum Cα–Cα distance to ApoB100, from **PDB 9BDE** (LDLR bound to ApoB100; LDLR chain R, residues 66–354)
- `is_in_helix`, `is_in_strand`: DSSP secondary structure, from 9BDE where the residue is resolved, otherwise from **PDB 1N7D** (LDLR extracellular domain at endosomal pH; chain A, residues 65–714)

Only LDLR and ApoB100 chains are used. Every resolved residue is checked against the variant table's reference amino acid before anything is written. Download the two structures from [rcsb.org](https://www.rcsb.org) (9BDE, 1N7D); DSSP (`mkdssp`) must be installed.

```bash
python ldlr_structure_features.py --pdb_9bde 9BDE.cif --pdb_1n7d 1N7D.cif \
    --csv ldlr_variant_scores_v5.csv --out structure_by_position.csv \
    --out_variants ldlr_variant_scores_v5.csv
```

### `ldlr_ridge_v5.py`: feature-based prediction
Predicts B and S for the 14,617 true missense variants from the 26-feature design:

- **Substitution (8):** BLOSUM62, Δhydrophobicity, Δvolume, Δpolarity, Δcharge, Cys / Pro / Gly flags
- **Domain (14):** domain one-hot (all 17 LDLR region labels) + VLDL blind-spot flag (LA2, LA6)
- **Structural (4):** distance to ApoB100, helix, strand, ΔΔG

MutateX ΔΔG was used by Tabet et al. but per-variant values were not published, so it is not yet included (25 features are used). Missing structural values are set to 0 with missing-data indicators.

Two models are fit, flat and with domain × feature interactions, using ridge regression with 5-fold cross-validation (scaling and penalty chosen inside each training fold). The script reports global and per-domain Pearson r and MAE, the variance-vs-predictability relationship, the failure analysis of the worst-predicted variants, and position-window scans (β-propeller, LA1–LA7), with tables and figures written to `--out_dir`.

```bash
python ldlr_ridge_v5.py --csv ldlr_variant_scores_v5.csv --out_dir ridge_results
```

---

## Data files

### `ldlr_variant_scores_v5.csv` (16,330 variants × 21 columns)
All variants scored in both Data S1 and S2: 14,617 missense, plus 774 stop-gain, 685 synonymous and 254 in-frame deletions. Analyses filter to missense variants in the scripts.

| Column | Description |
|---|---|
| `variant`, `position`, `ref_aa3`, `alt_aa3` | HGVS protein variant, residue number, reference / alternate amino acid (empty for synonymous and deletions) |
| `functional_score` | F, LDL uptake (Data S1); 0 = nonsense-like, 1 = synonymous-like |
| `abundance_score` | A, cell-surface abundance (Data S2), same scaling |
| `uptake_loess_predicted`, `B_score` | LOWESS fit F̂(A) and B = F − F̂(A) (`ldlr_compute_B.py`) |
| `domain` | LDLR region (17 labels: signal, linker_pre_LA1, LA1–LA7, EGF-A, EGF-B, beta-prop, EGF-C, O-sugar, linker_pre_TM, TM, NPxY) |
| `cluster`, `cluster_label` | Phenotype cluster in A–F space: 0 WT-like, 1 severe/null, 2 abundance_defect, 3 functional_defect |
| `mutation_type` | Substitution class (conservative, charge_gain_or_loss, charge_reversal, Cys_change, Gly_change, Pro_change, other_nonconservative; `unknown` for non-missense) |
| `delta_charge` | Charge change; Arg/Lys/His = +1, Asp/Glu = −1 |
| `is_vldl_blind_spot` | True for LA2 and LA6 |
| `route` | Cluster × region group used by A13 (C3_LA3_5, C2_beta_propeller, C3_beta_propeller, C3_cyto_tail, default) |
| `clinvar_vcep`, `clinvar_class` | ClinGen FH VCEP classification (Data S4); filled only for the ~120 expert-classified variants |
| `dist_apob100` | Distance to ApoB100 (Å), from 9BDE; filled for residues 66–354 only |
| `is_in_helix`, `is_in_strand` | Secondary structure (1/0); filled for residues 65–714 |
| `ss_source` | Structure providing the secondary structure: `9BDE` or `1N7D` (label only, not a model feature) |

Blank cells mean no value exists for that variant (not classified by the VCEP, or residue not resolved in either structure).

### `functional_abundance_clusters_k4.csv` (16,330 variants)
Input for A9: `variant`, `functional_score`, `abundance_score`, `cluster` (k = 4 phenotype clusters, same as in v5).

### Tabet et al. supplementary data
Only Data S1, S2 and S4 are used. The other supplementary files (S3 LDL uptake with excess VLDL, S5 evidence calibration, S6 oligonucleotides, S7/S8 PyMOL sessions) are not needed.

---

## Reproducing the results

```bash
# 1. B score (already in v5; optional check)
python ldlr_compute_B.py --check ldlr_variant_scores_v5.csv
# 2. Structural features (already in v5; optional rebuild)
python ldlr_structure_features.py --pdb_9bde 9BDE.cif --pdb_1n7d 1N7D.cif \
    --csv ldlr_variant_scores_v5.csv --out_variants ldlr_variant_scores_v5.csv
# 3. Ridge regression
python ldlr_ridge_v5.py --csv ldlr_variant_scores_v5.csv --out_dir ridge_results
# 4. Mechanistic model
python ldlr_model_a9.py p.Pro3Ile --csv functional_abundance_clusters_k4.csv
```

## Requirements

Python 3.10+, numpy, pandas, scipy, scikit-learn, statsmodels, matplotlib, biopython, [PyBioNetGen](https://pybionetgen.github.io/PyBioNetGen/) (`pip install bionetgen`), and DSSP (`mkdssp`) for the structure script only.

## References

- Tabet DR et al. The functional landscape of coding variation in the familial hypercholesterolemia gene *LDLR*. *Science* (2025). doi:10.1126/science.ady7186
- Reimund M et al. Structure of apolipoprotein B100 bound to the low-density lipoprotein receptor. *Nature* 638, 829–835 (2025). PDB 9BDE
- Rudenko G et al. Structure of the LDL receptor extracellular domain at endosomal pH. *Science* 298, 2353–2358 (2002). PDB 1N7D
