"""
LDLR model A13
=====================================================================

WHAT A13 DOES (in one file):
  1. Look up a variant (or a whole route) directly from the CSV.
  2. Compute scaling factors for that variant directly from its own
     abundance_score / B_score 
  3. Build and simulate the BNGL model with PyBioNetGen, printing
     progress at each stage (variant lookup -> BNGL write -> network
     generation + ODE simulation -> results parsed).
  4. Report every scaled parameter explicitly, with its route, before
     and after the simulation.

ROUTING TABLE — route -> factor used -> parameter(s) scaled
  ┌────────────────────┬──────────────────┬───────────────────────────┐
  │ Route               │ factor used       │ parameter(s) scaled       │
  ├────────────────────┼──────────────────┼───────────────────────────┤
  │ C3_LA3_5            │ binding factor     │ k_on_base                 │
  │ C2_beta_propeller   │ stability factor   │ LDLR_init, k_recycle_ldlr │
  │ C3_cyto_tail        │ endo factor        │ k_endo                    │
  │ C3_beta_propeller   │ (none for now   )  │ all held at WT base       │
  │ default             │ all three          │ LDLR_init, k_on_base,     │
  │                     │                    │ k_recycle_ldlr, k_endo    │
  └────────────────────┴──────────────────┴───────────────────────────┘

HOW EACH FACTOR IS COMPUTED — directly from the variant's own CSV row:
  binding factor    = max(0, 1 + B_score_raw / B_REF)
  stability factor  = clipped abundance_score
  endo factor       = clipped B_score

MODES:
  Single variant:
    python ldlr_model_a13.py p.Asp168Lys --csv ldlr_variant_scores_v5.csv

  Batch run over a route:
    python ldlr_model_a13.py --batch --subset C3_LA3_5 \\
      --csv ldlr_variant_scores_v5.csv --max 20

Run from anywhere — all paths are given as arguments; nothing is
hardcoded to a particular working directory.
"""

import argparse
import os
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import bionetgen


# =========================================================================
# Input table (overridable via CLI arg)
# =========================================================================
VARIANT_TABLE_CSV = "ldlr_variant_scores_v5.csv"

# =========================================================================
# Score handling
# =========================================================================
# B_score (residualized F) and abundance_score (A) can be negative in the
# raw data; clip to [0, 1.5] and floor at EPS for numerical validity.
F_MIN, F_MAX = 0.0, 1.5
A_MIN, A_MAX = 0.0, 1.5
EPS = 0.01

# Minimum value allowed for any scaling factor before it is used as a BNGL
# rate-law multiplier. BioNetGen errors on a zero/negative rate constant,
# so every factor is clamped to this floor in get_scaled_parameters()
# before it touches k_on_base, LDLR_init, or k_endo.
FACTOR_FLOOR = 1e-4

# =========================================================================
# WT base parameters
# =========================================================================
LDLR_INIT_BASE      = 1000
K_ON_BASE           = 1.0     # base LDL on-rate; scaled by binding factor in C3_LA3_5
B_REF               = 1.0     # reference B_score deviation for full binding loss
                               # (binding factor = max(0, 1 + B_raw/B_REF))
K_ENDO_BASE         = 2.0
K_OFF_ENDO_BASE     = 50.0    # hardcoded; direction TBD (see C3_beta_propeller)
K_RECYCLE_LDLR_BASE = 3.0
K_LYSO_LDL_BASE     = 3.0
K_DEGRADE_LDL_BASE  = 5.0

# =========================================================================
# WT steady-state references (for ratio-to-WT reporting)
# =========================================================================
WT_SURFACE_SS = 587.775
WT_UPTAKE_SS  = 621.906

# =========================================================================
# Allowed direct parameter overrides (sanity checks / LHS workflows)
# =========================================================================
ALLOWED_OVERRIDE_KEYS = {
    "LDLR_init",
    "k_recycle_ldlr",
    "k_lyso_ldl",
    "k_degrade_ldl",
    "k_endo",
    "k_off_surf",
    "k_off_endo",
    "k_on_base",
}


# =========================================================================
# Small utilities
# =========================================================================

@contextmanager
def pushd(new_dir: Path):
    prev = Path.cwd()
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(prev)


def normalize_variant_name(user_input: str) -> str:
    s = (user_input or "").strip().replace(" ", "")
    if not s:
        return s
    if not s.lower().startswith("p."):
        s = "p." + s
    return "p." + s[2:]


def clip_and_floor(x: float, lo: float, hi: float, eps: float = EPS) -> float:
    if pd.isna(x):
        raise ValueError("Score is NaN in table.")
    x = float(x)
    x = max(lo, min(hi, x))
    return max(x, eps)


def clip_scores(A_raw: float, F_raw: float, B_raw: float) -> tuple:
    """
    Clip and floor scores for numerical validity.
    Returns (A_used, F_used, B_used, note).
    """
    A_used = clip_and_floor(A_raw, lo=A_MIN, hi=A_MAX, eps=EPS)
    F_used = clip_and_floor(F_raw, lo=F_MIN, hi=F_MAX, eps=EPS)
    B_used = clip_and_floor(B_raw, lo=F_MIN, hi=F_MAX, eps=EPS)
    note = (
        f"A_raw={A_raw:.4f} -> A_used={A_used:.4f}; "
        f"F_raw={F_raw:.4f} -> F_used={F_used:.4f} (raw); "
        f"B_raw={B_raw:.4f} -> B_used={B_used:.4f} (residualized)"
    )
    return A_used, F_used, B_used, note


def lookup_variant(variant_name: str, csv_path: str = VARIANT_TABLE_CSV) -> dict:
    """
    Look up a single variant directly from the CSV (no other script needed).

    Returns:
        variant        : matched variant string (p.XxxNYyy)
        abundance_score: raw A score
        functional_score: raw F score
        B_score        : LOESS-residualized F
        cluster        : int cluster assignment
        cluster_label  : human-readable cluster name
        domain         : domain name (e.g. 'LA4', 'beta-prop', 'NPxY')
        mutation_type  : e.g. 'Cys_change', 'charge_gain_or_loss', etc.
        route          : pre-computed routing string from the CSV
    """
    df = pd.read_csv(csv_path)
    key = normalize_variant_name(variant_name)

    hit = df.loc[df["variant"] == key]
    if hit.empty:
        hit = df.loc[df["variant"].str.lower() == key.lower()]
    if hit.empty:
        raise KeyError(f"Variant not found: {variant_name!r} (normalized to {key!r})")

    row = hit.iloc[0]
    return {
        "variant":          str(row["variant"]),
        "abundance_score":  float(row["abundance_score"]),
        "functional_score": float(row["functional_score"]),
        "B_score":          float(row["B_score"]),
        "cluster":          int(row["cluster"]),
        "cluster_label":    str(row["cluster_label"]),
        "domain":           str(row["domain"]),
        "mutation_type":    str(row["mutation_type"]),
        "route":            str(row["route"]),
    }


def lookup_route_variants(route: str, csv_path: str) -> pd.DataFrame:
    """Return all variants for a given route, with the same filtering used
    everywhere else in this file (drop NaN scores, drop nonsense/Ter, drop
    VLDL blind spots)."""
    df = pd.read_csv(csv_path)
    sub = df[df["route"] == route].copy()
    sub = sub.dropna(subset=["abundance_score", "B_score"])
    sub = sub[sub["alt_aa3"] != "Ter"]
    sub = sub[~sub["is_vldl_blind_spot"]].reset_index(drop=True)
    return sub


# =========================================================================
# LDLR model — standalone: BNGL string, routing, ODE run, plotting
# =========================================================================

class LDLRModel:
    """
    Full LDLR receptor-trafficking model (4 LA-module binding, endocytosis,
    endosomal release, recycling, lysosomal degradation), simulated with
    PyBioNetGen. Self-contained — nothing else needs importing or running
    first.
    """

    def __init__(
        self,
        variant_name="WT",
        functional_score=1.0,     # clipped raw functional_score
        B_score=1.0,                # clipped B_score (residualized)
        B_score_raw=0.0,            # raw (unclipped) B_score — for binding factor
        abundance_score=1.0,
        cluster=0,
        cluster_label="WT-like",
        domain="unknown",
        mutation_type="unknown",
        route="default",
        LDLR_init_base=LDLR_INIT_BASE,
        wt_surface_ss=WT_SURFACE_SS,
        wt_uptake_ss=WT_UPTAKE_SS,
        param_overrides=None,
    ):
        self.variant_name     = variant_name
        self.functional_score = float(functional_score)
        self.B_score          = float(B_score)
        self.B_score_raw      = float(B_score_raw)
        self.abundance_score  = float(abundance_score)
        self.cluster          = int(cluster)
        self.cluster_label    = str(cluster_label)
        self.domain           = str(domain)
        self.mutation_type    = str(mutation_type)
        self.routing          = str(route)
        self.LDLR_init_base   = int(LDLR_init_base)
        self.wt_surface_ss    = float(wt_surface_ss)
        self.wt_uptake_ss     = float(wt_uptake_ss)
        self.param_overrides  = param_overrides or {}

        self.result            = None
        self.expected_gdat_key = None
        self.run_dir           = None

    # -- Scaling factors -----------------------------------------------------

    def get_scaling_factors(self) -> dict:
        """
        Compute the three scaling factors directly from this variant's own
        scores — the only source of scaling in A13:

          binding_factor   = max(0, 1 + B_score_raw / B_REF)
          stability_factor = clipped abundance_score
          endo_factor      = clipped B_score
        """
        binding_factor   = max(0.0, 1.0 + self.B_score_raw / B_REF)
        stability_factor = max(self.abundance_score, EPS)
        endo_factor      = max(self.B_score, EPS)
        return {
            "binding_factor":   binding_factor,
            "stability_factor": stability_factor,
            "endo_factor":      endo_factor,
        }

    # -- Parameter scaling ---------------------------------------------------

    def get_scaled_parameters(self) -> dict:
        """
        Route -> factor -> BNGL parameter:

          C3_LA3_5           binding_factor   -> k_on_base
                              (pure binding-affinity defect; k_off_surf stays WT)
          C2_beta_propeller  stability_factor -> LDLR_init, k_recycle_ldlr
                              (recycling / receptor-abundance defect)
          C3_cyto_tail       endo_factor      -> k_endo
                              (endocytosis-rate defect; NPxY/cytoplasmic tail)
          C3_beta_propeller  (none)           -> all params held at WT base
                              (endosomal-release defect; direction unresolved,
                              TODO?: scale k_off_endo once resolved)
          default            all three        -> LDLR_init, k_on_base,
                                                   k_recycle_ldlr, k_endo
                              (C0, C1, and any unmatched C2/C3 variant)
        """
        factors = self.get_scaling_factors()
        route   = self.routing

        fb = max(factors["binding_factor"],   FACTOR_FLOOR)
        fs = max(factors["stability_factor"], FACTOR_FLOOR)
        fe = max(factors["endo_factor"],      FACTOR_FLOOR)

        # Start every parameter at WT base
        LDLR_init      = self.LDLR_init_base
        k_on_base      = K_ON_BASE
        k_recycle_ldlr = K_RECYCLE_LDLR_BASE
        k_endo         = K_ENDO_BASE
        k_off_surf     = 1.0
        k_off_endo     = K_OFF_ENDO_BASE
        k_lyso_ldl     = K_LYSO_LDL_BASE
        k_degrade_ldl  = K_DEGRADE_LDL_BASE

        if route == "C3_LA3_5":
            k_on_base  = K_ON_BASE * fb
            k_off_surf = 1.0

        elif route == "C3_beta_propeller":
            # Intentional no-op: all params stay at WT base for this route.
            # not done yet: scale k_off_endo once scaling direction is resolved
            # from vector geometry.
            pass

        elif route == "C2_beta_propeller":
            LDLR_init      = int(round(self.LDLR_init_base * fs))
            k_recycle_ldlr = K_RECYCLE_LDLR_BASE * fs

        elif route == "C3_cyto_tail":
            k_endo = K_ENDO_BASE * fe

        else:  # default fallback (C0, C1, unmatched C2/C3)
            LDLR_init      = int(round(self.LDLR_init_base * fs))
            k_recycle_ldlr = K_RECYCLE_LDLR_BASE * fs
            k_on_base      = K_ON_BASE * fb
            k_endo         = K_ENDO_BASE * fe

        params = {
            "route":            route,
            "binding_factor":   fb,
            "stability_factor": fs,
            "endo_factor":      fe,
            "A_eff_model":      max(self.abundance_score, EPS),
            "F_used":           max(self.functional_score, EPS),
            "B_used":           max(self.B_score, EPS),
            "LDLR_init":        LDLR_init,
            "k_on_base":        k_on_base,
            "k_recycle_ldlr":   k_recycle_ldlr,
            "k_lyso_ldl":       k_lyso_ldl,
            "k_degrade_ldl":    k_degrade_ldl,
            "k_endo":           k_endo,
            "k_off_surf":       k_off_surf,
            "k_off_endo":       k_off_endo,
        }

        # Direct overrides (sanity checks / LHS workflows) applied last
        for key, value in self.param_overrides.items():
            if key not in ALLOWED_OVERRIDE_KEYS:
                raise KeyError(
                    f"Unknown parameter override: {key!r}. "
                    f"Allowed: {sorted(ALLOWED_OVERRIDE_KEYS)}"
                )
            params[key] = float(value)

        # Enforce valid ranges
        params["LDLR_init"] = max(int(round(params["LDLR_init"])), 1)
        for key in ["k_on_base", "k_recycle_ldlr", "k_lyso_ldl", "k_degrade_ldl",
                    "k_endo", "k_off_surf", "k_off_endo"]:
            params[key] = max(float(params[key]), 1e-9)

        return params

    # -- BNGL model string ---------------------------------------------------

    def get_model_string(self) -> str:
        """
        The full BNGL model: 4 LA-module surface binding, endocytosis,
        endosomal release, receptor recycling, and lysosomal LDL
        degradation. Reaction structure is fixed; only the parameter
        VALUES (computed above) differ per variant.
        """
        p = self.get_scaled_parameters()
        return f"""
# LDLR Model A13 (standalone) - {self.variant_name}
# Cluster:          {self.cluster} ({self.cluster_label})
# Domain:           {self.domain}
# Mutation type:    {self.mutation_type}
# Route:            {p["route"]}
# binding_factor:   {p["binding_factor"]:.4f}
# stability_factor: {p["stability_factor"]:.4f}
# endo_factor:      {p["endo_factor"]:.4f}
# Direct overrides: {self.param_overrides if self.param_overrides else "None"}

begin model

begin parameters
    # Binding / uptake-side
    k_on_base          {p["k_on_base"]:.6f}
    k_off_surf         {p["k_off_surf"]:.6f}
    k_off_endo         {p["k_off_endo"]:.6f}
    k_endo             {p["k_endo"]:.6f}

    # Module strengths (LA4 most critical; LA4 > LA3 ~ LA5 > LA7)
    strength_LA3       1.0
    strength_LA4       1.0
    strength_LA5       1.0
    strength_LA7       0.8

    # Trafficking / degradation parameters
    k_recycle_ldlr     {p["k_recycle_ldlr"]:.6f}
    k_lyso_ldl         {p["k_lyso_ldl"]:.6f}
    k_degrade_ldl      {p["k_degrade_ldl"]:.6f}

    # Initial supply
    LDLR_init          {p["LDLR_init"]}
    LDL_conc           100
end parameters

begin functions
  k_on_LA3() = k_on_base*strength_LA3
  k_on_LA4() = k_on_base*strength_LA4
  k_on_LA5() = k_on_base*strength_LA5
  k_on_LA7() = k_on_base*strength_LA7
end functions

begin molecule types
    LDLR(la3,la4,la5,la7,loc~surface~endosome)
    LDL(ldlr,loc~extra~endo~lyso)
end molecule types

begin seed species
    LDLR(la3,la4,la5,la7,loc~surface)  LDLR_init
    LDL(ldlr,loc~extra)  LDL_conc
end seed species

begin observables
    Molecules  LDLR_surface      LDLR(loc~surface)
    Molecules  LDLR_endosome     LDLR(loc~endosome)
    Molecules  LDL_free          LDL(ldlr,loc~extra)
    Molecules  LDL_endo          LDL(ldlr,loc~endo)
    Molecules  LDL_lyso          LDL(ldlr,loc~lyso)
    Molecules  Surf_LA3          LDLR(la3!+,loc~surface)
    Molecules  Surf_LA4          LDLR(la4!+,loc~surface)
    Molecules  Surf_LA5          LDLR(la5!+,loc~surface)
    Molecules  Surf_LA7          LDLR(la7!+,loc~surface)
end observables

begin reaction rules
    # Surface binding / unbinding
    LDLR(la3,la4,la5,la7,loc~surface) + LDL(ldlr,loc~extra) -> \\
        LDLR(la3!1,la4,la5,la7,loc~surface).LDL(ldlr!1,loc~extra) + LDL(ldlr,loc~extra) \\
        k_on_LA3()
    LDLR(la3!1,la4,la5,la7,loc~surface).LDL(ldlr!1,loc~extra) -> \\
        LDLR(la3,la4,la5,la7,loc~surface) \\
        k_off_surf

    LDLR(la3,la4,la5,la7,loc~surface) + LDL(ldlr,loc~extra) -> \\
        LDLR(la3,la4!1,la5,la7,loc~surface).LDL(ldlr!1,loc~extra) + LDL(ldlr,loc~extra) \\
        k_on_LA4()
    LDLR(la3,la4!1,la5,la7,loc~surface).LDL(ldlr!1,loc~extra) -> \\
        LDLR(la3,la4,la5,la7,loc~surface) \\
        k_off_surf

    LDLR(la3,la4,la5,la7,loc~surface) + LDL(ldlr,loc~extra) -> \\
        LDLR(la3,la4,la5!1,la7,loc~surface).LDL(ldlr!1,loc~extra) + LDL(ldlr,loc~extra) \\
        k_on_LA5()
    LDLR(la3,la4,la5!1,la7,loc~surface).LDL(ldlr!1,loc~extra) -> \\
        LDLR(la3,la4,la5,la7,loc~surface) \\
        k_off_surf

    LDLR(la3,la4,la5,la7,loc~surface) + LDL(ldlr,loc~extra) -> \\
        LDLR(la3,la4,la5,la7!1,loc~surface).LDL(ldlr!1,loc~extra) + LDL(ldlr,loc~extra) \\
        k_on_LA7()
    LDLR(la3,la4,la5,la7!1,loc~surface).LDL(ldlr!1,loc~extra) -> \\
        LDLR(la3,la4,la5,la7,loc~surface) \\
        k_off_surf

    # Endocytosis
    LDLR(la3!1,la4,la5,la7,loc~surface).LDL(ldlr!1,loc~extra) -> \\
        LDLR(la3!1,la4,la5,la7,loc~endosome).LDL(ldlr!1,loc~endo)  k_endo
    LDLR(la3,la4!1,la5,la7,loc~surface).LDL(ldlr!1,loc~extra) -> \\
        LDLR(la3,la4!1,la5,la7,loc~endosome).LDL(ldlr!1,loc~endo)  k_endo
    LDLR(la3,la4,la5!1,la7,loc~surface).LDL(ldlr!1,loc~extra) -> \\
        LDLR(la3,la4,la5!1,la7,loc~endosome).LDL(ldlr!1,loc~endo)  k_endo
    LDLR(la3,la4,la5,la7!1,loc~surface).LDL(ldlr!1,loc~extra) -> \\
        LDLR(la3,la4,la5,la7!1,loc~endosome).LDL(ldlr!1,loc~endo)  k_endo

    # Endosomal release
    LDLR(la3!1,la4,la5,la7,loc~endosome).LDL(ldlr!1,loc~endo) -> \\
        LDLR(la3,la4,la5,la7,loc~endosome) + LDL(ldlr,loc~endo)  k_off_endo
    LDLR(la3,la4!1,la5,la7,loc~endosome).LDL(ldlr!1,loc~endo) -> \\
        LDLR(la3,la4,la5,la7,loc~endosome) + LDL(ldlr,loc~endo)  k_off_endo
    LDLR(la3,la4,la5!1,la7,loc~endosome).LDL(ldlr!1,loc~endo) -> \\
        LDLR(la3,la4,la5,la7,loc~endosome) + LDL(ldlr,loc~endo)  k_off_endo
    LDLR(la3,la4,la5,la7!1,loc~endosome).LDL(ldlr!1,loc~endo) -> \\
        LDLR(la3,la4,la5,la7,loc~endosome) + LDL(ldlr,loc~endo)  k_off_endo

    # Receptor recycling: endosome -> surface
    LDLR(la3,la4,la5,la7,loc~endosome) -> \\
        LDLR(la3,la4,la5,la7,loc~surface)  k_recycle_ldlr

    # LDL cargo: lysosomal routing and degradation
    LDL(ldlr,loc~endo) -> LDL(ldlr,loc~lyso)  k_lyso_ldl
    LDL(ldlr,loc~lyso) -> 0                    k_degrade_ldl

end reaction rules

end model
"""

    # -- Run / results ---------------------------------------------------

    def run(self, t_end=200, n_steps=1000, out_root="results/data/results",
            quiet=False, verbose=True):
        """
        Writes the BNGL file, then runs BioNetGen (network generation +
        ODE simulation). Progress is printed at each stage; set quiet=True
        to also suppress BioNetGen's own console output (network-generation
        and integrator messages), or verbose=False to suppress this
        script's own step announcements.
        """
        out_root = Path(out_root)
        out_root.mkdir(parents=True, exist_ok=True)

        safe_name = self.variant_name.replace("/", "_")
        run_dir   = out_root / safe_name
        run_dir.mkdir(parents=True, exist_ok=True)
        self.run_dir = run_dir

        temp_model_file        = run_dir / f"{safe_name}_temp.bngl"
        self.expected_gdat_key = temp_model_file.stem

        if verbose:
            print(f"[A13] Step 1/4: variant={self.variant_name}  route={self.routing}")

        with open(temp_model_file, "w") as f:
            f.write(self.get_model_string())
            f.write("\ngenerate_network({overwrite=>1})\n")
            f.write(f"simulate({{method=>\"ode\", t_end=>{t_end}, n_steps=>{n_steps}}})\n")

        if verbose:
            print(f"[A13] Step 2/4: BNGL written to {temp_model_file}")
            print(f"[A13] Step 3/4: running BioNetGen "
                  f"(network generation + ODE integration, t_end={t_end}, n_steps={n_steps}) ...")

        with pushd(run_dir):
            self.result = bionetgen.run(temp_model_file.name, out=".", suppress=quiet)

        if verbose:
            print(f"[A13] Step 4/4: simulation complete for {self.variant_name}")

        return self.result

    def _pick_correct_gdat_key(self) -> str:
        keys = list(self.result.gdats.keys())
        if self.expected_gdat_key in self.result.gdats:
            return self.expected_gdat_key
        for k in keys:
            if self.variant_name in k:
                return k
        return keys[0]

    def get_data(self):
        if self.result is None:
            raise ValueError("Must run simulation first.")

        model_name = self._pick_correct_gdat_key()
        df = pd.DataFrame(self.result.gdats[model_name])

        endo_cols = [c for c in df.columns if c.lower() == "ldl_endo"]
        lyso_cols = [c for c in df.columns if c.lower() == "ldl_lyso"]
        if not endo_cols or not lyso_cols:
            raise KeyError("No endosome or lysosome columns found in GDAT!")

        df["LDL_internalized"]     = df[endo_cols].sum(axis=1) + df[lyso_cols].sum(axis=1)
        df["LDLR_surface_complex"] = df[["Surf_LA3", "Surf_LA4", "Surf_LA5", "Surf_LA7"]].sum(axis=1)

        final_surface = float(df["LDLR_surface"].iloc[-1])
        final_uptake  = float(df["LDL_internalized"].iloc[-1])

        df["surface_ratio_vs_WTss"] = final_surface / self.wt_surface_ss
        df["uptake_ratio_vs_WTss"]  = final_uptake  / self.wt_uptake_ss

        return df

    def get_summary_metrics(self) -> dict:
        data = self.get_data()

        final_surface = float(data["LDLR_surface"].iloc[-1])
        final_uptake  = float(data["LDL_internalized"].iloc[-1])

        return {
            "final_surface":         final_surface,
            "final_uptake":          final_uptake,
            "surface_ratio_vs_WTss": final_surface / self.wt_surface_ss,
            "uptake_ratio_vs_WTss":  final_uptake  / self.wt_uptake_ss,
            "target_A":              self.abundance_score,
            "target_F":              self.functional_score,
            "target_B":              self.B_score,
        }

    def plot(self, save_path=None):
        data    = self.get_data()
        metrics = self.get_summary_metrics()
        params  = self.get_scaled_parameters()

        fig, axes = plt.subplots(2, 2, figsize=(13, 10))

        axes[0, 0].plot(data["time"], data["LDLR_surface"],  label="Surface LDLR",  linewidth=2)
        axes[0, 0].plot(data["time"], data["LDLR_endosome"], label="Endosome LDLR", linewidth=2)
        axes[0, 0].axhline(self.wt_surface_ss, linestyle="--", linewidth=1,
                           label=f"WT ss surface={self.wt_surface_ss:.1f}")
        axes[0, 0].set_title(
            f"Surface LDLR\nfinal/WTss={metrics['surface_ratio_vs_WTss']:.3f}  "
            f"target A={metrics['target_A']:.3f}")
        axes[0, 0].legend(); axes[0, 0].grid(alpha=0.3)

        axes[0, 1].plot(data["time"], data["LDL_free"],         label="Free LDL",         linewidth=2)
        axes[0, 1].plot(data["time"], data["LDL_internalized"], label="Internalized LDL",  linewidth=2, color="red")
        axes[0, 1].set_title("LDL Dynamics")
        axes[0, 1].legend(); axes[0, 1].grid(alpha=0.3)

        axes[1, 0].plot(data["time"], data["LDL_internalized"], linewidth=2, color="darkred")
        axes[1, 0].axhline(self.wt_uptake_ss, linestyle="--", linewidth=1,
                           label=f"WT ss uptake={self.wt_uptake_ss:.1f}")
        axes[1, 0].set_title(
            f"Steady-state LDL Uptake\nfinal/WTss={metrics['uptake_ratio_vs_WTss']:.3f}  "
            f"F(raw)={metrics['target_F']:.3f}  B(resid)={metrics['target_B']:.3f}")
        axes[1, 0].legend(); axes[1, 0].grid(alpha=0.3)

        axes[1, 1].plot(data["time"], data["LDLR_surface_complex"], linewidth=2, color="purple")
        axes[1, 1].set_title("LDLR-LDL Complexes at Surface")
        axes[1, 1].grid(alpha=0.3)

        override_text = ("none" if not self.param_overrides else
                         ", ".join(f"{k}={v:.4g}" if isinstance(v, (int, float)) else f"{k}={v}"
                                   for k, v in self.param_overrides.items()))

        fig.suptitle(
            f"{self.variant_name}  |  cluster={self.cluster} ({self.cluster_label})  "
            f"domain={self.domain}  mutation_type={self.mutation_type}\n"
            f"Route: {self.routing}  |  "
            f"factors=[binding={params['binding_factor']:.3f}, "
            f"stability={params['stability_factor']:.3f}, endo={params['endo_factor']:.3f}]\n"
            f"Overrides: {override_text}\n"
            f"LDLR_init={params['LDLR_init']}  k_on_base={params['k_on_base']:.3f}  "
            f"k_recycle_ldlr={params['k_recycle_ldlr']:.3f}  k_endo={params['k_endo']:.3f}  "
            f"k_off_surf={params['k_off_surf']:.3f}  k_off_endo={params['k_off_endo']:.3f} [hardcoded]",
            fontsize=10, fontweight="bold",
        )
        plt.tight_layout()

        if save_path:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"[A13] Plot saved to: {save_path}")

        plt.close(fig)
        return fig


# =========================================================================
# Helper to build a model instance from a CSV row
# =========================================================================

def build_model(info: dict, param_overrides: dict | None = None) -> LDLRModel:
    """Given a variant's CSV row (from lookup_variant), build an LDLRModel.
    Scaling factors are computed inside the model directly from the
    variant's own scores — nothing else needs to be loaded."""
    A_used, F_used, B_used, _ = clip_scores(
        info["abundance_score"], info["functional_score"], info["B_score"]
    )

    return LDLRModel(
        variant_name=info["variant"],
        functional_score=F_used,
        B_score=B_used,
        B_score_raw=info["B_score"],
        abundance_score=A_used,
        cluster=info["cluster"],
        cluster_label=info["cluster_label"],
        domain=info["domain"],
        mutation_type=info["mutation_type"],
        route=info["route"],
        param_overrides=param_overrides,
    )


# =========================================================================
# Single-variant mode
# =========================================================================

def run_single(variant_name: str, csv_path: str, t_end: float, n_steps: int,
                out_dir: str, overrides: dict | None = None, quiet: bool = False):

    info = lookup_variant(variant_name, csv_path=csv_path)

    print("=== Variant identification ===")
    print(f"  Matched variant : {info['variant']}")
    print(f"  Cluster         : {info['cluster']} ({info['cluster_label']})")
    print(f"  Domain          : {info['domain']}")
    print(f"  Mutation type   : {info['mutation_type']}")
    print(f"  Route           : {info['route']}")

    model = build_model(info, param_overrides=overrides)
    params = model.get_scaled_parameters()

    print("\n=== Scaling factors (computed directly from this variant's scores) ===")
    print(f"  binding_factor    = {params['binding_factor']:.4f}")
    print(f"  stability_factor  = {params['stability_factor']:.4f}")
    print(f"  endo_factor       = {params['endo_factor']:.4f}")

    print("\n=== Scaled BNGL parameters ===")
    print(f"  Route             = {params['route']}")
    print(f"  LDLR_init         = {params['LDLR_init']}")
    print(f"  k_on_base         = {params['k_on_base']:.6f}")
    print(f"  k_recycle_ldlr    = {params['k_recycle_ldlr']:.6f}")
    print(f"  k_lyso_ldl        = {params['k_lyso_ldl']:.6f}")
    print(f"  k_degrade_ldl     = {params['k_degrade_ldl']:.6f}")
    print(f"  k_endo            = {params['k_endo']:.6f}")
    print(f"  k_off_surf        = {params['k_off_surf']:.6f}")
    print(f"  k_off_endo        = {params['k_off_endo']:.6f}  [hardcoded; TODO A14]")
    if overrides:
        print(f"  Direct overrides  = {overrides}")

    print()
    model.run(t_end=t_end, n_steps=n_steps, out_root=out_dir, quiet=quiet)
    metrics = model.get_summary_metrics()

    print("\n=== Steady-state comparison to WT references ===")
    print(f"  Surface ratio vs WT ss = {metrics['surface_ratio_vs_WTss']:.3f}  "
          f"(observed A={metrics['target_A']:.3f})")
    print(f"  Uptake ratio vs WT ss  = {metrics['uptake_ratio_vs_WTss']:.3f}  "
          f"(observed F={metrics['target_F']:.3f}, B={metrics['target_B']:.3f})")

    save_fig = str(Path(out_dir) / info["variant"] / f"{info['variant']}_a13.png")
    model.plot(save_path=save_fig)


# =========================================================================
# Batch mode
# =========================================================================

def run_batch(subset: str, csv_path: str, out_dir: Path, t_end=200, n_steps=500,
              max_variants=None, quiet=True):

    print(f"\n{'=' * 65}")
    print(f"A13 Batch ODE Run — {subset}")
    print(f"{'=' * 65}")

    sub = lookup_route_variants(subset, csv_path)
    if max_variants:
        sub = sub.head(max_variants)
    print(f"Variants: {len(sub)}")

    results = []
    for i, row in sub.iterrows():
        info = {
            "variant":          row["variant"],
            "abundance_score":  row["abundance_score"],
            "functional_score": row["functional_score"],
            "B_score":          row["B_score"],
            "cluster":          int(row["cluster"]),
            "cluster_label":    str(row["cluster_label"]),
            "domain":           str(row["domain"]),
            "mutation_type":    str(row["mutation_type"]),
            "route":            subset,
        }

        model = build_model(info)
        S_pred = B_pred = np.nan
        try:
            model.run(t_end=t_end, n_steps=n_steps,
                      out_root=str(out_dir / "ode_runs"), quiet=quiet, verbose=False)
            m = model.get_summary_metrics()
            S_pred = m["surface_ratio_vs_WTss"]
            B_pred = m["uptake_ratio_vs_WTss"] - S_pred
        except Exception as e:
            print(f"  ODE failed for {info['variant']}: {e}")

        results.append({
            "variant": info["variant"],
            "domain":  info["domain"],
            "S_obs":   info["abundance_score"],
            "B_obs":   info["B_score"],
            "S_pred":  S_pred,
            "B_pred":  B_pred,
        })

        if len(results) % 5 == 0:
            print(f"  {len(results)}/{len(sub)} done")

    result_df = pd.DataFrame(results)

    def r(x, y):
        m = np.isfinite(x) & np.isfinite(y)
        return float(np.corrcoef(x[m], y[m])[0, 1]) if m.sum() > 2 else np.nan

    S_obs = result_df["S_obs"].values
    B_obs = result_df["B_obs"].values

    print(f"\n{'=' * 65}")
    print(f"RESULTS — {subset}")
    print(f"{'=' * 65}")
    print(f"r_S = {r(S_obs, result_df['S_pred'].values):.3f}")
    print(f"r_B = {r(B_obs, result_df['B_pred'].values):.3f}")

    out_csv = out_dir / f"a13_batch_{subset}.csv"
    result_df.to_csv(out_csv, index=False)
    print(f"\nSaved: {out_csv}")

    fig, ax = plt.subplots(figsize=(6, 5))
    mask = np.isfinite(result_df["B_pred"].values)
    ax.scatter(B_obs[mask], result_df["B_pred"].values[mask], s=10, alpha=0.5,
               color="#028090", edgecolors="none")
    ax.plot([-1.6, 0.5], [-1.6, 0.5], "k--", lw=0.8, alpha=0.5)
    ax.set_xlabel("Observed B (Tabet)")
    ax.set_ylabel("Predicted B (ODE)")
    ax.set_title(f"A13 — {subset}")
    ax.grid(alpha=0.2)
    plot_path = out_dir / f"a13_batch_{subset}.png"
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved: {plot_path}")

    return result_df


# =========================================================================
# Main
# =========================================================================

def main():
    ap = argparse.ArgumentParser(
        description="LDLR model A13 — standalone, cluster+domain-aware routing")
    ap.add_argument("variant", nargs="?", default=None,
                    help="Single variant, e.g. p.Asp168Lys or Asp168Lys")
    ap.add_argument("--batch",      action="store_true")
    ap.add_argument("--subset",     default="C3_LA3_5")
    ap.add_argument("--csv",        default=VARIANT_TABLE_CSV)
    ap.add_argument("--t_end",      type=float, default=200.0)
    ap.add_argument("--n_steps",    type=int,   default=1000)
    ap.add_argument("--max",        type=int,   default=None)
    ap.add_argument("--out_dir",    default="results/a13")
    ap.add_argument("--quiet",      action="store_true",
                    help="Suppress BioNetGen's own console progress output.")

    for key in sorted(ALLOWED_OVERRIDE_KEYS):
        ap.add_argument(f"--{key}", type=float, default=None)

    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.batch:
        run_batch(args.subset, args.csv, out_dir, args.t_end, args.n_steps,
                  args.max, quiet=args.quiet)
    elif args.variant:
        overrides = {k: getattr(args, k) for k in ALLOWED_OVERRIDE_KEYS
                     if getattr(args, k, None) is not None}
        run_single(args.variant, args.csv, args.t_end, args.n_steps, str(out_dir),
                   overrides=overrides or None, quiet=args.quiet)
    else:
        print("Specify a variant or use --batch")
        print("Examples:")
        print("  python ldlr_model_a13.py p.Asp168Lys --csv ldlr_variant_scores_v5.csv")
        print("  python ldlr_model_a13.py --batch --subset C3_LA3_5 --max 20")
        ap.print_help()


if __name__ == "__main__":
    main()
