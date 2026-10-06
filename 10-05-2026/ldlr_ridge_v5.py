"""
LDLR ridge regression — can variant features predict B and S?
=================================================================

PURPOSE
  Tests whether a fixed set of variant features (the 26-feature design) carries enough
  information to predict the two phenotype components of LDLR missense
  variants (Tabet et al., Science 2025), without the mechanistic model:
    S = abundance_score                       (surface abundance, = A)
    B = B_score = F - LOESS(F | A)            (function beyond abundance)
  This script produces every ridge-regression number used in the
  presentation: global r, domain-level r and MAE, the variance-vs-r
  relationship, the β-propeller failure analysis and the position-window
  scans.

INPUT
  ldlr_variant_scores_v5.csv (16,330 rows). Only true missense variants
  are analysed: rows with no ref/alt amino acid (synonymous, in-frame
  deletions) and stop-gain rows (alt = Ter) are removed -> 14,617 variants.
  The structural columns (dist_apob100, is_in_helix, is_in_strand) are
  computed by ldlr_structure_features.py from PDB 9BDE and PDB 1N7D; the
  column ss_source records which structure the secondary structure came
  from and is not used by the model.

THE 26-FEATURE DESIGN (25 used in this version)
  The design has 26 features: 8 substitution + 14 domain + 4 structural.
  The fourth structural feature, MutateX ΔΔG, is not included: Tabet et
  al. used MutateX ΔΔG in their analysis but did not publish per-variant
  values. It will be added once the values are available; until then 25
  features are used. Domain identity is encoded with all
  17 region labels in the CSV, so the model receives 8 + 17 + 1 + 3 = 29
  feature columns, plus 2 missing-data indicators (31 columns, flat model).

  A. Substitution properties (8), computed from ref_aa3 / alt_aa3:
       blosum62               BLOSUM62[ref, alt]
       delta_hydrophobicity   Kyte-Doolittle(alt) - Kyte-Doolittle(ref)
       delta_volume           Grantham volume(alt) - volume(ref)
       delta_polarity         Grantham polarity(alt) - polarity(ref)
       delta_charge           charge(alt) - charge(ref); R/K/H = +1, D/E = -1
       is_cys_change          1 if ref or alt is Cys
       is_pro_intro           1 if alt is Pro
       is_gly_change          1 if ref or alt is Gly
  B. Domain identity (17 one-hot + 1):
       domain                 one-hot over the 17 region labels in the CSV:
                              signal, linker_pre_LA1, LA1-LA7, EGF-A,
                              EGF-B, beta-prop, EGF-C, O-sugar,
                              linker_pre_TM, TM, NPxY
       is_vldl_blind_spot     1 if domain is LA2 or LA6
  C. Structural (3 used):
       dist_apob100           min Cα–Cα distance to ApoB100
                              (PDB 9BDE, LDLR residues 66-354)
       is_in_helix            DSSP helix (H, G, I)
       is_in_strand           DSSP β-strand (E)
                              (PDB 9BDE where resolved, otherwise PDB 1N7D;
                              together LDLR residues 65-714)
       (MutateX ΔΔG           not included, see above)

MISSING STRUCTURAL DATA
  dist_apob100 exists for ~34% of missense variants (9BDE covers LA2 to
  EGF-A) and secondary structure for ~74% (9BDE + 1N7D cover LA2 to
  EGF-C). Signal, LA1, O-sugar, TM and the cytoplasmic tail are in neither
  structure. Missing values are set to 0, and two indicator columns
  (structure_missing for distance, ss_missing for secondary structure) let
  the model tell "missing" from "value 0". The indicators are bookkeeping
  for missing data, not additional biological features.

MODELS
  Flat          : all feature columns above (+ missing-data indicators).
  + interactions: flat model plus domain × feature terms for the 11
                  non-domain features (A and C), so each feature can have
                  a different slope in each domain.
  Both: standardization -> ridge regression, penalty chosen by internal
  cross-validation. Performance = Pearson r between observed values and
  5-fold out-of-fold predictions (KFold, shuffled, random_state = 0).
  Standardization and penalty selection are refit inside every training
  fold, so no information from the held-out fold is used.

OUTPUTS (written to --out_dir, default ridge_results/)
  summary.txt                 all headline numbers
  oof_predictions.csv         per-variant observed and predicted B and S
  domain_metrics.csv          per-domain n, r (95% bootstrap CI), MAE, SD
  window_beta_propeller.csv   20-residue window scan, β-propeller
  window_LA_repeats.csv       10-residue window scan, LA1-LA7
  fig_domain_r.png            local r by domain (B and S)
  fig_variance_vs_r.png       within-domain SD vs local r (B and S)
  fig_window_beta_propeller.png, fig_window_LA_repeats.png

USAGE
  python ldlr_ridge_v5.py --csv ldlr_variant_scores_v5.csv --out_dir ridge_results

REQUIREMENTS
  numpy, pandas, scipy, scikit-learn, matplotlib, biopython
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.Align import substitution_matrices
from scipy.stats import pearsonr
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


# =========================================================================
# Settings
# =========================================================================
N_FOLDS      = 5
RANDOM_STATE = 0
ALPHAS       = np.logspace(-2, 4, 25)   # ridge penalty grid
N_BOOT       = 1000                     # bootstrap resamples for CIs
WORST_FRAC   = 0.05                     # failure analysis: worst / best 5%

DOMAIN_ORDER = ["signal", "linker_pre_LA1", "LA1", "LA2", "LA3", "LA4", "LA5",
                "LA6", "LA7", "EGF-A", "EGF-B", "beta-prop", "EGF-C",
                "O-sugar", "linker_pre_TM", "TM", "NPxY"]
LA_DOMAINS = ["LA1", "LA2", "LA3", "LA4", "LA5", "LA6", "LA7"]

# Amino-acid property scales
AA3_TO_1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
            "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
            "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
            "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
KYTE_DOOLITTLE = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
                  "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
                  "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
                  "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
GRANTHAM_VOLUME = {"A": 31, "R": 124, "N": 56, "D": 54, "C": 55, "Q": 85,
                   "E": 83, "G": 3, "H": 96, "I": 111, "L": 111, "K": 119,
                   "M": 105, "F": 132, "P": 32.5, "S": 32, "T": 61,
                   "W": 170, "Y": 136, "V": 84}
GRANTHAM_POLARITY = {"A": 8.1, "R": 10.5, "N": 11.6, "D": 13.0, "C": 5.5,
                     "Q": 10.5, "E": 12.3, "G": 9.0, "H": 10.4, "I": 5.2,
                     "L": 4.9, "K": 11.3, "M": 5.7, "F": 5.2, "P": 8.0,
                     "S": 9.2, "T": 8.6, "W": 5.4, "Y": 6.2, "V": 5.9}

SUBSTITUTION_FEATS = ["blosum62", "delta_hydrophobicity", "delta_volume",
                      "delta_polarity", "delta_charge", "is_cys_change",
                      "is_pro_intro", "is_gly_change"]
STRUCTURAL_FEATS   = ["dist_apob100", "is_in_helix", "is_in_strand"]
MISSING_FLAGS      = ["structure_missing", "ss_missing"]
INTERACT_FEATS     = SUBSTITUTION_FEATS + STRUCTURAL_FEATS   # 11 features


# =========================================================================
# 1. Data
# =========================================================================

def load_missense(csv_path: str) -> pd.DataFrame:
    """Read the variant table and keep true missense variants only."""
    df = pd.read_csv(csv_path, low_memory=False)
    n_all = len(df)
    df = df.dropna(subset=["ref_aa3", "alt_aa3"])        # synonymous, deletions
    df = df[df["alt_aa3"] != "Ter"]                      # stop-gain
    df = df[df["ref_aa3"] != df["alt_aa3"]]
    df = df.dropna(subset=["B_score", "abundance_score"]).reset_index(drop=True)
    print(f"Variants: {n_all} rows in CSV -> {len(df)} true missense variants")
    return df


def build_features(df: pd.DataFrame) -> pd.DataFrame:
    """Compute the features (+ missing-data indicators)."""
    df = df.copy()
    ref = df["ref_aa3"].map(AA3_TO_1)
    alt = df["alt_aa3"].map(AA3_TO_1)
    if ref.isna().any() or alt.isna().any():
        raise ValueError("Unrecognised amino-acid code in ref_aa3 / alt_aa3.")

    # A. Substitution properties
    blosum = substitution_matrices.load("BLOSUM62")
    df["blosum62"] = [float(blosum[r, a]) for r, a in zip(ref, alt)]
    df["delta_hydrophobicity"] = alt.map(KYTE_DOOLITTLE) - ref.map(KYTE_DOOLITTLE)
    df["delta_volume"] = alt.map(GRANTHAM_VOLUME) - ref.map(GRANTHAM_VOLUME)
    df["delta_polarity"] = alt.map(GRANTHAM_POLARITY) - ref.map(GRANTHAM_POLARITY)
    # delta_charge is stored in the CSV (R/K/H = +1, D/E = -1)
    df["is_cys_change"] = ((ref == "C") | (alt == "C")).astype(int)
    df["is_pro_intro"] = ((alt == "P") & (ref != "P")).astype(int)
    df["is_gly_change"] = ((ref == "G") | (alt == "G")).astype(int)

    # B. Domain identity
    df["is_vldl_blind_spot"] = df["is_vldl_blind_spot"].astype(int)

    # C. Structural (9BDE / 1N7D), with missing-data indicators
    df["structure_missing"] = df["dist_apob100"].isna().astype(int)
    df["ss_missing"] = df["is_in_helix"].isna().astype(int)
    for col in STRUCTURAL_FEATS:
        df[col] = df[col].fillna(0.0)
    return df


def design_matrices(df: pd.DataFrame):
    """Return (X_flat, X_interactions) as DataFrames."""
    dom = pd.get_dummies(df["domain"], prefix="dom").astype(float)
    base = pd.concat([df[SUBSTITUTION_FEATS + STRUCTURAL_FEATS + MISSING_FLAGS
                         + ["is_vldl_blind_spot"]].astype(float), dom], axis=1)
    inter = {f"{d}__x__{f}": dom[d].values * df[f].values
             for d in dom.columns for f in INTERACT_FEATS}
    X_int = pd.concat([base, pd.DataFrame(inter, index=df.index)], axis=1)
    return base, X_int


# =========================================================================
# 2. Model
# =========================================================================

def oof_predict(X: pd.DataFrame, y: np.ndarray) -> np.ndarray:
    """5-fold out-of-fold predictions; scaling and penalty fit per fold."""
    model = make_pipeline(StandardScaler(), RidgeCV(alphas=ALPHAS))
    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    return cross_val_predict(model, X.values, y, cv=kf)


def r_and_ci(obs, pred, n_boot=N_BOOT, seed=0):
    obs, pred = np.asarray(obs, float), np.asarray(pred, float)
    if len(obs) < 8 or obs.std() == 0 or pred.std() == 0:
        return np.nan, np.nan, np.nan
    r = pearsonr(obs, pred)[0]
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, len(obs), len(obs))
        if obs[i].std() > 0 and pred[i].std() > 0:
            boots.append(pearsonr(obs[i], pred[i])[0])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return r, lo, hi


# =========================================================================
# 3. Analyses
# =========================================================================

def domain_metrics(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dom in DOMAIN_ORDER:
        g = df[df["domain"] == dom]
        if g.empty:
            continue
        rB, rB_lo, rB_hi = r_and_ci(g["B_score"], g["pred_B"])
        rS, rS_lo, rS_hi = r_and_ci(g["abundance_score"], g["pred_S"])
        rows.append(dict(
            domain=dom, n=len(g),
            r_B=rB, r_B_lo=rB_lo, r_B_hi=rB_hi,
            r_S=rS, r_S_lo=rS_lo, r_S_hi=rS_hi,
            MAE_B=(g["B_score"] - g["pred_B"]).abs().mean(),
            MAE_S=(g["abundance_score"] - g["pred_S"]).abs().mean(),
            SD_B=g["B_score"].std(), SD_S=g["abundance_score"].std(),
        ))
    return pd.DataFrame(rows)


def window_scan(df, domain, width, starts=None):
    """Local r in consecutive position windows inside one domain."""
    g = df[df["domain"] == domain]
    lo, hi = int(g["position"].min()), int(g["position"].max())
    if starts is None:
        starts = list(range(lo, hi + 1, width))
    rows = []
    for k, s in enumerate(starts):
        e = s + width - 1
        if k == len(starts) - 1:
            e = max(e, hi)                       # last window takes the remainder
        w = g[(g["position"] >= s) & (g["position"] <= e)]
        rB = pearsonr(w["B_score"], w["pred_B"])[0] if len(w) > 8 else np.nan
        rS = pearsonr(w["abundance_score"], w["pred_S"])[0] if len(w) > 8 else np.nan
        rows.append(dict(domain=domain, start=s, end=min(e, hi), n=len(w), r_B=rB, r_S=rS))
    return pd.DataFrame(rows)


def failure_analysis(df: pd.DataFrame) -> dict:
    """Which variants are worst / best predicted for B?"""
    err = (df["B_score"] - df["pred_B"]).abs()
    k = int(round(WORST_FRAC * len(df)))
    worst = df.loc[err.nlargest(k).index]
    best = df.loc[err.nsmallest(k).index]
    return dict(
        bp_share_all=(df["domain"] == "beta-prop").mean(),
        bp_share_worst=(worst["domain"] == "beta-prop").mean(),
        cluster_all=df["cluster_label"].value_counts(normalize=True),
        cluster_worst=worst["cluster_label"].value_counts(normalize=True),
        wt_like_best=(best["cluster_label"] == "WT-like").mean(),
    )


# =========================================================================
# 4. Figures
# =========================================================================

def plot_domain_r(dm, path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, col, title, color in [(axes[0], "r_B", "B_score", "#028090"),
                                  (axes[1], "r_S", "abundance_score", "#D2691E")]:
        d = dm.sort_values(col)
        ax.barh(d["domain"], d[col], color=color)
        for y, v in enumerate(d[col]):
            ax.text(v + (0.01 if v >= 0 else -0.01), y, f"{v:.3f}",
                    va="center", ha="left" if v >= 0 else "right", fontsize=8)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_xlim(-0.35, 0.9)
        ax.set_title(f"{title}\nlocal r by domain")
        ax.set_xlabel("r")
    plt.tight_layout(); plt.savefig(path, dpi=150); plt.close(fig)


def plot_variance_vs_r(dm, path):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, sd, r, title, color in [(axes[0], "SD_B", "r_B", "Binding (B)", "#028090"),
                                    (axes[1], "SD_S", "r_S", "Abundance (S)", "#D2691E")]:
        ax.scatter(dm[sd], dm[r], color=color)
        for _, row in dm.iterrows():
            ax.annotate(row["domain"], (row[sd], row[r]), fontsize=7,
                        xytext=(3, 3), textcoords="offset points")
        cc = pearsonr(dm[sd], dm[r])[0]
        ax.set_title(f"{title}: variance vs. predictability (r = {cc:.2f})")
        ax.set_xlabel("Within-domain SD of true score")
        ax.set_ylabel("Local r")
        ax.axhline(0, color="black", lw=0.6)
    plt.tight_layout(); plt.savefig(path, dpi=150); plt.close(fig)


def plot_windows(win, path, by_domain):
    if not by_domain:
        fig, ax = plt.subplots(figsize=(10, 4))
        labels = [f"{s}-{e}" for s, e in zip(win["start"], win["end"])]
        ax.plot(labels, win["r_B"], "o-", color="#028090", label="r_B")
        ax.plot(labels, win["r_S"], "o-", color="#D2691E", label="r_S")
        ax.axhline(0, color="black", lw=0.6)
        ax.set_ylabel("local r"); ax.legend()
        plt.xticks(rotation=45, ha="right")
    else:
        fig, axes = plt.subplots(2, 4, figsize=(14, 6), sharey=True)
        for ax, dom in zip(axes.flat, LA_DOMAINS):
            w = win[win["domain"] == dom]
            ax.plot(w["start"], w["r_B"], "o-", color="#028090", label="r_B")
            ax.plot(w["start"], w["r_S"], "o-", color="#D2691E", label="r_S")
            ax.axhline(0, color="black", lw=0.6)
            ax.set_title(dom); ax.set_xticks(w["start"])
        axes.flat[-1].axis("off")
        axes.flat[0].legend()
    plt.tight_layout(); plt.savefig(path, dpi=150); plt.close(fig)


# =========================================================================
# Main
# =========================================================================

def main():
    ap = argparse.ArgumentParser(description="LDLR ridge regression (26-feature design)")
    ap.add_argument("--csv", default="ldlr_variant_scores_v5.csv")
    ap.add_argument("--out_dir", default="ridge_results")
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    df = build_features(load_missense(args.csv))
    X_flat, X_int = design_matrices(df)
    y_B, y_S = df["B_score"].values, df["abundance_score"].values
    print(f"Design matrix: flat {X_flat.shape[1]} columns, "
          f"with interactions {X_int.shape[1]} columns")

    # Global performance
    res = {}
    for name, X in [("flat", X_flat), ("interactions", X_int)]:
        pB, pS = oof_predict(X, y_B), oof_predict(X, y_S)
        res[name] = (pearsonr(y_B, pB)[0], pearsonr(y_S, pS)[0], pB, pS)
    df["pred_B"], df["pred_S"] = res["interactions"][2], res["interactions"][3]

    # Downstream analyses use the interaction model
    dm = domain_metrics(df)
    var_r_B = pearsonr(dm["SD_B"], dm["r_B"])[0]
    var_r_S = pearsonr(dm["SD_S"], dm["r_S"])[0]
    fail = failure_analysis(df)
    win_bp = window_scan(df, "beta-prop", 20)
    win_la = pd.concat([window_scan(df, d, 10,
                        starts=[int(df.loc[df.domain == d, "position"].min()) + 10 * k
                                for k in range(4)]) for d in LA_DOMAINS])

    # Save tables
    df[["variant", "position", "domain", "cluster_label", "B_score", "pred_B",
        "abundance_score", "pred_S"]].to_csv(out / "oof_predictions.csv", index=False)
    dm.round(4).to_csv(out / "domain_metrics.csv", index=False)
    win_bp.round(4).to_csv(out / "window_beta_propeller.csv", index=False)
    win_la.round(4).to_csv(out / "window_LA_repeats.csv", index=False)

    # Figures
    plot_domain_r(dm, out / "fig_domain_r.png")
    plot_variance_vs_r(dm, out / "fig_variance_vs_r.png")
    plot_windows(win_bp, out / "fig_window_beta_propeller.png", by_domain=False)
    plot_windows(win_la, out / "fig_window_LA_repeats.png", by_domain=True)

    # Summary
    lines = [
        "LDLR ridge regression (26-feature design, 25 used) — summary",
        f"n = {len(df)} true missense variants; {N_FOLDS}-fold CV, random_state = {RANDOM_STATE}",
        "",
        "GLOBAL PERFORMANCE (Pearson r, out-of-fold)",
        f"  Flat            r_B = {res['flat'][0]:.3f}   r_S = {res['flat'][1]:.3f}",
        f"  + interactions  r_B = {res['interactions'][0]:.3f}   r_S = {res['interactions'][1]:.3f}",
        "",
        "LOCAL r BY DOMAIN (interaction model)",
        dm[["domain", "n", "r_B", "r_S", "MAE_B", "MAE_S"]].round(3).to_string(index=False),
        "",
        "VARIANCE vs PREDICTABILITY across domains",
        f"  corr(within-domain SD, local r):  B = {var_r_B:.2f}   S = {var_r_S:.2f}",
        "",
        "FAILURE ANALYSIS (B)",
        f"  β-propeller share of all variants        = {fail['bp_share_all']:.1%}",
        f"  β-propeller share of worst-predicted 5%  = {fail['bp_share_worst']:.1%}",
        "  Cluster composition, all variants:",
        fail["cluster_all"].map("{:.1%}".format).to_string(),
        "  Cluster composition, worst-predicted 5%:",
        fail["cluster_worst"].map("{:.1%}".format).to_string(),
        f"  WT-like share of best-predicted 5%       = {fail['wt_like_best']:.1%}",
    ]
    text = "\n".join(lines)
    (out / "summary.txt").write_text(text)
    print("\n" + text + f"\n\nAll outputs written to {out.resolve()}")


if __name__ == "__main__":
    main()
