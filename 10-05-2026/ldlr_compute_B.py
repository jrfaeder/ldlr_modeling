"""
B score: LDL-uptake function beyond what surface abundance explains
====================================================================

PURPOSE
  Computes, for every LDLR variant measured in both Tabet et al. (Science
  2025) maps,
    A  = abundance_score   cell-surface abundance       (Data S2, "score")
    F  = functional_score  LDL uptake                   (Data S1, "score")
    F_hat(A)               LOWESS fit of F on A
    B  = F - F_hat(A)      function not explained by abundance
  B > 0: more uptake than expected for the receptor's surface level;
  B < 0: less uptake than expected (a per-receptor function defect).

METHOD
  1. Keep variants scored in both S1 and S2 (16,330 variants).
  2. Fit LOWESS of F on A using true missense variants only (14,617;
     synonymous, in-frame deletions and stop-gain variants excluded),
     with frac = 0.3 and 3 robustifying iterations (statsmodels defaults
     otherwise).
  3. Evaluate the fitted curve at every variant's A by linear
     interpolation, so all 16,330 variants receive F_hat and B.

  These settings reproduce the columns uptake_loess_predicted and
  B_score in the project variant table exactly (max difference 0).

INPUT
  science_ady7186_data_s1.csv   LDL uptake scores (Tabet Data S1)
  science_ady7186_data_s2.csv   surface abundance scores (Tabet Data S2)

OUTPUT
  --out  CSV with: variant, position, ref_aa3, alt_aa3, is_missense,
         abundance_score, functional_score, uptake_loess_predicted, B_score

USAGE
  python ldlr_compute_B.py --s1 science_ady7186_data_s1.csv \\
      --s2 science_ady7186_data_s2.csv --out B_scores.csv
  # optional: confirm against an existing variant table
  python ldlr_compute_B.py ... --check ldlr_variant_scores_v5.csv

REQUIREMENTS
  numpy, pandas, statsmodels
"""

import argparse
import re

import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

LOWESS_FRAC = 0.3
LOWESS_ITER = 3
HGVS = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|del)$")


def parse_hgvs(variant: str):
    """p.Gly2Leu -> ('Gly', 2, 'Leu'); synonymous / deletion -> (None, pos, None)."""
    m = HGVS.match(variant)
    if not m:
        return None, np.nan, None
    ref, pos, alt = m.groups()
    if alt in {"=", "del"}:
        return None, int(pos), None
    return ref, int(pos), alt


def build_table(s1_path: str, s2_path: str) -> pd.DataFrame:
    f = pd.read_csv(s1_path)[["hgvs_pro", "score"]].rename(columns={"score": "functional_score"})
    a = pd.read_csv(s2_path)[["hgvs_pro", "score"]].rename(columns={"score": "abundance_score"})
    df = f.merge(a, on="hgvs_pro", how="inner").rename(columns={"hgvs_pro": "variant"})
    df = df.dropna(subset=["functional_score", "abundance_score"]).reset_index(drop=True)

    parsed = df["variant"].map(parse_hgvs)
    df["ref_aa3"] = [p[0] for p in parsed]
    df["position"] = [p[1] for p in parsed]
    df["alt_aa3"] = [p[2] for p in parsed]
    df["is_missense"] = (df["ref_aa3"].notna() & df["alt_aa3"].notna()
                         & (df["alt_aa3"] != "Ter") & (df["ref_aa3"] != df["alt_aa3"]))
    return df


def compute_B(df: pd.DataFrame) -> pd.DataFrame:
    fit_rows = df[df["is_missense"]]
    curve = lowess(fit_rows["functional_score"].values,
                   fit_rows["abundance_score"].values,
                   frac=LOWESS_FRAC, it=LOWESS_ITER, return_sorted=True)
    df = df.copy()
    df["uptake_loess_predicted"] = np.interp(df["abundance_score"].values,
                                             curve[:, 0], curve[:, 1])
    df["B_score"] = df["functional_score"] - df["uptake_loess_predicted"]
    print(f"Variants scored in S1 and S2: {len(df)}; "
          f"LOWESS fit on {len(fit_rows)} missense variants "
          f"(frac={LOWESS_FRAC}, it={LOWESS_ITER})")
    return df


def main():
    ap = argparse.ArgumentParser(description="Compute B = F - LOWESS(F | A)")
    ap.add_argument("--s1", default="science_ady7186_data_s1.csv")
    ap.add_argument("--s2", default="science_ady7186_data_s2.csv")
    ap.add_argument("--out", default="B_scores.csv")
    ap.add_argument("--check", default=None,
                    help="optional variant table to compare B_score against")
    args = ap.parse_args()

    df = compute_B(build_table(args.s1, args.s2))
    cols = ["variant", "position", "ref_aa3", "alt_aa3", "is_missense",
            "abundance_score", "functional_score", "uptake_loess_predicted", "B_score"]
    df[cols].to_csv(args.out, index=False)
    print(f"Wrote {args.out}")

    if args.check:
        ref = pd.read_csv(args.check, low_memory=False)[["variant", "B_score"]]
        m = df.merge(ref, on="variant", suffixes=("", "_ref"))
        diff = (m["B_score"] - m["B_score_ref"]).abs().max()
        print(f"Check vs {args.check}: {len(m)} variants matched, "
              f"max |B difference| = {diff:.2e}")


if __name__ == "__main__":
    main()
