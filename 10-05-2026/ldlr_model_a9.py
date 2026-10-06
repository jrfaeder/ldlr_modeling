"""
LDLR model A9 — BioNetGen model of LDL uptake driven by variant A and F scores
==============================================================================

PURPOSE
  Simulates the LDLR uptake cycle for one LDLR missense variant with
  PyBioNetGen and compares the simulated steady state with the variant's
  measured scores (Tabet et al., Science 2025):
    A (abundance_score)  : surface LDLR abundance relative to WT
    F (functional_score) : LDL uptake relative to WT

  The same script is the simulation engine for the parameter-space
  exploration: Latin Hypercube Sampling, sensitivity analysis and inverse
  fitting all call it with direct parameter overrides.

MODEL STRUCTURE (BNGL)
  - LDLR has four LDL-binding modules (LA3, LA4, LA5, LA7) and two
    locations (surface, endosome). LDL has three locations
    (extracellular, endosome, lysosome).
  - Reactions: surface binding / unbinding (k_on, k_off_surf),
    endocytosis of LDLR–LDL complexes (k_endo), endosomal LDL release
    (k_off_endo), receptor recycling endosome -> surface (k_recycle_ldlr),
    LDL routing to the lysosome (k_lyso_ldl) and LDL degradation
    (k_degrade_ldl).
  - Receptors are recycled, not degraded, so the receptor pool is set by
    LDLR_init. Extracellular LDL is held constant (LDL_conc).
  - Receptor recycling and LDL cargo trafficking have separate rate
    constants, so receptor fate and LDL fate can change independently.

HOW A VARIANT'S SCORES SET THE PARAMETERS
  Scores are clipped to [0, 1.5] and floored at 0.01 so every rate is
  positive (some measured scores are slightly negative). Then:
    A  -> receptor supply      LDLR_init      = 1000 * A
                               k_recycle_ldlr = 3.0  * A
    F  -> uptake efficiency    k_endo         = 2.0  * F
                               k_off_surf     = 1.0  / F
    k_off_endo, k_lyso_ldl and k_degrade_ldl stay at their WT values.
  The phenotype cluster is read from the CSV and reported, but it does not
  change the scaling: every variant is treated with the same rule.
  With A = F = 1 the model reproduces the WT steady state.

OUTPUTS
  Simulated steady-state surface LDLR and internalized LDL (endosome +
  lysosome), each divided by the WT steady state (WT_SURFACE_SS,
  WT_UPTAKE_SS), are compared with the target A and F. A 4-panel time-
  course figure is saved for each run.

DIRECT PARAMETER OVERRIDES
  Any rate constant (or LDLR_init) can be set directly from the command
  line. Overrides are applied after score-based scaling, which is how the
  LHS, sensitivity and inverse-fitting workflows explore parameter space.

USAGE
  python ldlr_model_a9.py p.Pro3Ile
  python ldlr_model_a9.py p.Pro3Ile --csv functional_abundance_clusters_k4.csv
  python ldlr_model_a9.py p.Pro3Ile --t_end 200 --n_steps 1000
  python ldlr_model_a9.py p.Pro3Ile --k_endo 1.2 --k_off_surf 0.8

  Input CSV columns: variant, functional_score, abundance_score, cluster.
"""

import argparse
import os
from contextlib import contextmanager
from pathlib import Path

import bionetgen
import pandas as pd
import matplotlib.pyplot as plt


# =========================
# Input table
# =========================
VARIANT_TABLE_CSV = "functional_abundance_clusters_k4.csv"

# =========================
# Score handling
# =========================
# Scores are clipped to [MIN, MAX] and floored at EPS so that every
# scaled rate constant stays positive.
F_MIN, F_MAX = 0.0, 1.5
A_MIN, A_MAX = 0.0, 1.5
EPS = 0.01

# =========================
# WT base parameters
# =========================
LDLR_INIT_BASE = 1000
K_ENDO_BASE = 2.0
K_OFF_ENDO_BASE = 50.0

# =========================
# Trafficking / degradation (WT values)
# =========================
K_RECYCLE_LDLR_BASE = 3.0   # receptor recycling
K_LYSO_LDL_BASE = 3.0       # LDL cargo routing to lysosome
K_DEGRADE_LDL_BASE = 5.0    # LDL degradation in lysosome

# =========================
# WT steady-state references
# (surface LDLR and internalized LDL at A = F = 1; used as ratio denominators)
# =========================
WT_SURFACE_SS = 587.775
WT_UPTAKE_SS = 621.906

# =========================
# Parameters that can be set directly (LHS / sensitivity / fitting)
# =========================
ALLOWED_OVERRIDE_KEYS = {
    "LDLR_init",
    "k_recycle_ldlr",
    "k_lyso_ldl",
    "k_degrade_ldl",
    "k_endo",
    "k_off_surf",
    "k_off_endo",
}


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


def lookup_variant_scores(variant_name: str, csv_path: str = VARIANT_TABLE_CSV) -> dict:
    df = pd.read_csv(csv_path)
    key = normalize_variant_name(variant_name)

    hit = df.loc[df["variant"] == key]
    if hit.empty:
        hit = df.loc[df["variant"].str.lower() == key.lower()]

    if hit.empty:
        raise KeyError(f"Variant not found: {variant_name} (normalized to '{key}')")

    row = hit.iloc[0]
    return {
        "variant": str(row["variant"]),
        "functional_score": float(row["functional_score"]),
        "abundance_score": float(row["abundance_score"]),
        "cluster": int(row["cluster"]),
    }


def choose_scores_by_cluster(cluster: int, F_raw: float, A_raw: float):
    """
    Prepare a variant's scores for the model.

    F and A are clipped to [0, 1.5] and floored at EPS so that the scaled
    rate constants stay positive. The same rule applies to every variant;
    the cluster label is passed through for reporting only.

    Returns:
        F_used, A_used, note, skip_sim, cluster
        (skip_sim is always False; it is kept so callers can use one
        return signature.)
    """
    F_used = clip_and_floor(F_raw, lo=F_MIN, hi=F_MAX, eps=EPS)
    A_used = clip_and_floor(A_raw, lo=A_MIN, hi=A_MAX, eps=EPS)
    skip_sim = False
    note = (
        f"Cluster {cluster} (reported only; not used in scaling). "
        f"F_raw={float(F_raw):.4f} -> F_used={F_used:.4f}; "
        f"A_raw={float(A_raw):.4f} -> A_used={A_used:.4f}."
    )
    return F_used, A_used, note, skip_sim, cluster


class LDLRModel:
    def __init__(
        self,
        variant_name="WT",
        functional_score=1.0,
        abundance_score=1.0,
        cluster=0,
        LDLR_init_base=LDLR_INIT_BASE,
        wt_surface_ss=WT_SURFACE_SS,
        wt_uptake_ss=WT_UPTAKE_SS,
        param_overrides=None,
    ):
        self.variant_name = variant_name
        self.functional_score = float(functional_score)
        self.abundance_score = float(abundance_score)
        self.cluster = int(cluster)
        self.LDLR_init_base = int(LDLR_init_base)

        self.wt_surface_ss = float(wt_surface_ss)
        self.wt_uptake_ss = float(wt_uptake_ss)

        self.param_overrides = param_overrides or {}

        self.result = None
        self.expected_gdat_key = None
        self.run_dir = None

    def get_effective_abundance(self) -> float:
        """Abundance used for scaling: A_used, floored at EPS."""
        return max(self.abundance_score, EPS)

    def get_scaled_parameters(self):
        """
        Map the variant's scores onto BNGL parameters.

        A -> receptor supply:
             LDLR_init      = LDLR_init_base * A
             k_recycle_ldlr = K_RECYCLE_LDLR_BASE * A
        F -> uptake efficiency:
             k_endo         = K_ENDO_BASE * F
             k_off_surf     = 1 / F   (WT value 1.0)
        Held at WT:
             k_off_endo, k_lyso_ldl, k_degrade_ldl

        param_overrides, if given, replace any of these values afterwards.
        """
        F = max(self.functional_score, EPS)
        A_eff_model = self.get_effective_abundance()

        # Receptor-side abundance parameters
        LDLR_init_scaled = int(round(self.LDLR_init_base * A_eff_model))
        k_recycle_ldlr_scaled = K_RECYCLE_LDLR_BASE * A_eff_model

        # LDL cargo parameters kept constant
        k_lyso_ldl_scaled = K_LYSO_LDL_BASE
        k_degrade_ldl_scaled = K_DEGRADE_LDL_BASE

        # Uptake/function-side parameters
        k_endo_scaled = K_ENDO_BASE * F
        k_off_surf_scaled = 1.0 / F

        params = {
            "A_eff_model": A_eff_model,
            "LDLR_init": LDLR_init_scaled,
            "k_recycle_ldlr": k_recycle_ldlr_scaled,
            "k_lyso_ldl": k_lyso_ldl_scaled,
            "k_degrade_ldl": k_degrade_ldl_scaled,
            "k_endo": k_endo_scaled,
            "k_off_surf": k_off_surf_scaled,
            "k_off_endo": K_OFF_ENDO_BASE,
        }

        # Apply direct overrides
        for key, value in self.param_overrides.items():
            if key not in ALLOWED_OVERRIDE_KEYS:
                raise KeyError(
                    f"Unknown parameter override: {key}. "
                    f"Allowed keys: {sorted(ALLOWED_OVERRIDE_KEYS)}"
                )
            params[key] = float(value)

        # Keep LDLR_init as valid positive integer
        params["LDLR_init"] = int(round(params["LDLR_init"]))
        params["LDLR_init"] = max(params["LDLR_init"], 1)

        # Keep positive rates
        for key in ["k_recycle_ldlr", "k_lyso_ldl", "k_degrade_ldl", "k_endo", "k_off_surf", "k_off_endo"]:
            params[key] = max(float(params[key]), 1e-9)

        return params

    def get_model_string(self):
        p = self.get_scaled_parameters()

        return f"""
# LDLR model A9 - {self.variant_name}
# Cluster: {self.cluster}
# F_used: {self.functional_score:.3f}
# A_used_raw: {self.abundance_score:.3f}
# A_used_model: {p["A_eff_model"]:.3f}
# Direct overrides: {self.param_overrides if self.param_overrides else "None"}

begin model

begin parameters
    # Binding / uptake-side
    k_on_base          1.0
    k_off_surf         {p["k_off_surf"]:.6f}
    k_off_endo         {p["k_off_endo"]:.6f}
    k_endo             {p["k_endo"]:.6f}

    # Module strengths
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
        LDLR(la3,la4,la5,la7!1,loc~surface).LDL(ldlr!1,loc~extra) + LDL(ldlr,loc~extra)\\
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

    # Receptor-side recycling
    LDLR(la3,la4,la5,la7,loc~endosome) -> \\
        LDLR(la3,la4,la5,la7,loc~surface)  k_recycle_ldlr

    # LDL cargo-side lysosomal routing / degradation
    LDL(ldlr,loc~endo) -> LDL(ldlr,loc~lyso)  k_lyso_ldl
    LDL(ldlr,loc~lyso) -> 0  k_degrade_ldl

end reaction rules

end model
"""

    def run(self, t_end=200, n_steps=1000, out_root="results/data/results"):
        out_root = Path(out_root)
        out_root.mkdir(parents=True, exist_ok=True)

        safe_variant_name = self.variant_name.replace("/", "_")
        run_dir = out_root / safe_variant_name
        run_dir.mkdir(parents=True, exist_ok=True)
        self.run_dir = run_dir

        temp_model_file = run_dir / f"{safe_variant_name}_temp.bngl"
        self.expected_gdat_key = temp_model_file.stem

        with open(temp_model_file, "w") as f:
            f.write(self.get_model_string())
            f.write("\ngenerate_network({overwrite=>1})\n")
            f.write(f"simulate({{method=>\"ode\", t_end=>{t_end}, n_steps=>{n_steps}}})\n")

        print(f"BNGL file written: {temp_model_file}")

        with pushd(run_dir):
            self.result = bionetgen.run(temp_model_file.name, out=".", suppress=True)

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
            raise ValueError("Must run simulation first")

        model_name = self._pick_correct_gdat_key()
        gdat_array = self.result.gdats[model_name]
        df = pd.DataFrame(gdat_array)

        endo_cols = [c for c in df.columns if c.lower() == "ldl_endo"]
        lyso_cols = [c for c in df.columns if c.lower() == "ldl_lyso"]
        if not endo_cols or not lyso_cols:
            raise KeyError("No endosome or lysosome columns found in GDAT!")

        df["LDL_internalized"] = df[endo_cols].sum(axis=1) + df[lyso_cols].sum(axis=1)
        df["LDLR_surface_complex"] = df[["Surf_LA3", "Surf_LA4", "Surf_LA5", "Surf_LA7"]].sum(axis=1)

        final_surface = float(df["LDLR_surface"].iloc[-1])
        final_uptake = float(df["LDL_internalized"].iloc[-1])

        df["surface_ratio_vs_WTss"] = final_surface / self.wt_surface_ss
        df["uptake_ratio_vs_WTss"] = final_uptake / self.wt_uptake_ss

        print(f"Reading simulation output: {model_name}")
        return df

    def get_summary_metrics(self):
        data = self.get_data()

        final_surface = float(data["LDLR_surface"].iloc[-1])
        final_uptake = float(data["LDL_internalized"].iloc[-1])

        surface_ratio = final_surface / self.wt_surface_ss
        uptake_ratio = final_uptake / self.wt_uptake_ss

        return {
            "final_surface": final_surface,
            "final_uptake": final_uptake,
            "surface_ratio_vs_WTss": surface_ratio,
            "uptake_ratio_vs_WTss": uptake_ratio,
            "target_A": self.abundance_score,
            "target_F": self.functional_score,
        }

    def plot(self, save_path=None):
        data = self.get_data()
        metrics = self.get_summary_metrics()
        params = self.get_scaled_parameters()

        fig, axes = plt.subplots(2, 2, figsize=(13, 10))

        axes[0, 0].plot(data["time"], data["LDLR_surface"], label="Surface LDLR", linewidth=2)
        axes[0, 0].plot(data["time"], data["LDLR_endosome"], label="Endosome LDLR", linewidth=2)
        axes[0, 0].axhline(
            self.wt_surface_ss,
            linestyle="--",
            linewidth=1,
            label=f"WT ss surface={self.wt_surface_ss:.3f}"
        )
        axes[0, 0].set_title(
            f"Surface LDLR\nfinal/WTss = {metrics['surface_ratio_vs_WTss']:.3f}, "
            f"target A = {metrics['target_A']:.3f}"
        )
        axes[0, 0].legend()
        axes[0, 0].grid(alpha=0.3)

        axes[0, 1].plot(data["time"], data["LDL_free"], label="Free LDL", linewidth=2)
        axes[0, 1].plot(data["time"], data["LDL_internalized"], label="Internalized LDL", linewidth=2, color="red")
        axes[0, 1].set_title("LDL Dynamics")
        axes[0, 1].legend()
        axes[0, 1].grid(alpha=0.3)

        axes[1, 0].plot(data["time"], data["LDL_internalized"], linewidth=2, color="darkred")
        axes[1, 0].axhline(
            self.wt_uptake_ss,
            linestyle="--",
            linewidth=1,
            label=f"WT ss uptake={self.wt_uptake_ss:.3f}"
        )
        axes[1, 0].set_title(
            f"Steady-state LDL Uptake\nfinal/WTss = {metrics['uptake_ratio_vs_WTss']:.3f}, "
            f"target F = {metrics['target_F']:.3f}"
        )
        axes[1, 0].legend()
        axes[1, 0].grid(alpha=0.3)

        axes[1, 1].plot(data["time"], data["LDLR_surface_complex"], linewidth=2, color="purple")
        axes[1, 1].set_title("LDLR-LDL Complexes at Surface")
        axes[1, 1].grid(alpha=0.3)

        override_text = "none" if not self.param_overrides else ", ".join(
            [f"{k}={v:.4g}" if isinstance(v, (int, float)) else f"{k}={v}" for k, v in self.param_overrides.items()]
        )

        fig.suptitle(
            f"{self.variant_name} (cluster={self.cluster}, F_used={self.functional_score:.2f}, "
            f"A_used={self.abundance_score:.2f})\n"
            f"Overrides: {override_text}\n"
            f"LDLR_init={params['LDLR_init']}, k_recycle_ldlr={params['k_recycle_ldlr']:.3f}, "
            f"k_endo={params['k_endo']:.3f}, k_off_surf={params['k_off_surf']:.3f}",
            fontsize=12,
            fontweight="bold",
        )
        plt.tight_layout()

        if save_path:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"✓ Plot saved to: {save_path}")

        plt.show()
        return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("variant", help="Variant name, e.g. Gly2Leu or p.Gly2Leu")
    ap.add_argument("--csv", default=VARIANT_TABLE_CSV)
    ap.add_argument("--t_end", type=float, default=200.0)
    ap.add_argument("--n_steps", type=int, default=1000)
    ap.add_argument(
        "--out_dir",
        default="results/data/results",
        help="Root output directory; a subfolder per variant will be created.",
    )
    ap.add_argument("--save_fig", default=None)
    ap.add_argument("--LDLR_init_base", type=int, default=1000)
    ap.add_argument("--wt_surface_ss", type=float, default=WT_SURFACE_SS)
    ap.add_argument("--wt_uptake_ss", type=float, default=WT_UPTAKE_SS)

    # Optional direct overrides
    ap.add_argument("--LDLR_init", type=float, default=None)
    ap.add_argument("--k_recycle_ldlr", type=float, default=None)
    ap.add_argument("--k_lyso_ldl", type=float, default=None)
    ap.add_argument("--k_degrade_ldl", type=float, default=None)
    ap.add_argument("--k_endo", type=float, default=None)
    ap.add_argument("--k_off_surf", type=float, default=None)
    ap.add_argument("--k_off_endo", type=float, default=None)

    args = ap.parse_args()

    info = lookup_variant_scores(args.variant, csv_path=args.csv)
    F_used, A_used, note, skip_sim, cluster = choose_scores_by_cluster(
        info["cluster"], info["functional_score"], info["abundance_score"]
    )

    print("=== Variant treatment ===")
    print(f"Matched variant: {info['variant']}")
    print(f"Cluster:         {cluster}")
    print(f"F_raw={info['functional_score']:.4f} -> F_used={F_used:.4f}")
    print(f"A_raw={info['abundance_score']:.4f} -> A_used={A_used:.4f}")
    print(note)

    if skip_sim:
        print("=== Skipped simulation ===")
        return

    overrides = {}
    for key in ALLOWED_OVERRIDE_KEYS:
        value = getattr(args, key, None)
        if value is not None:
            overrides[key] = value

    if overrides:
        print("=== Direct parameter overrides requested ===")
        for k, v in overrides.items():
            print(f"{k:18s} = {v}")

    model = LDLRModel(
        variant_name=info["variant"],
        functional_score=F_used,
        abundance_score=A_used,
        cluster=cluster,
        LDLR_init_base=args.LDLR_init_base,
        wt_surface_ss=args.wt_surface_ss,
        wt_uptake_ss=args.wt_uptake_ss,
        param_overrides=overrides if overrides else None,
    )

    model.run(t_end=args.t_end, n_steps=args.n_steps, out_root=args.out_dir)

    if args.save_fig is None:
        args.save_fig = str(Path(args.out_dir) / info["variant"] / f"{info['variant']}.png")

    params = model.get_scaled_parameters()
    print("=== Final parameters used in simulation ===")
    print(f"A_eff_model       = {params['A_eff_model']:.4f}")
    print(f"LDLR_init         = {params['LDLR_init']}")
    print(f"k_recycle_ldlr    = {params['k_recycle_ldlr']:.6f}")
    print(f"k_lyso_ldl        = {params['k_lyso_ldl']:.6f}")
    print(f"k_degrade_ldl     = {params['k_degrade_ldl']:.6f}")
    print(f"k_endo            = {params['k_endo']:.6f}")
    print(f"k_off_surf        = {params['k_off_surf']:.6f}")
    print(f"k_off_endo        = {params['k_off_endo']:.6f}")

    metrics = model.get_summary_metrics()
    print("=== Steady-state comparison to WT references ===")
    print(f"Final surface LDLR          = {metrics['final_surface']:.3f}")
    print(f"Surface ratio vs WT ss      = {metrics['surface_ratio_vs_WTss']:.3f}")
    print(f"Target A                    = {metrics['target_A']:.3f}")
    print(f"Final LDL uptake            = {metrics['final_uptake']:.3f}")
    print(f"Uptake ratio vs WT ss       = {metrics['uptake_ratio_vs_WTss']:.3f}")
    print(f"Target F                    = {metrics['target_F']:.3f}")

    model.plot(save_path=args.save_fig)


if __name__ == "__main__":
    main()