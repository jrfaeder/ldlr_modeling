"""
LDLR model using PyBioNetGen
Simplified version with 4 LA modules
Cluster-aware mapping of Functional score (F) and Abundance score (A)

Method 3B:
- No cluster-2 abundance rescue
- Cluster 2 uses both A and F directly
- For cluster 2 only, strengthen uptake rescue by:
    k_endo = K_ENDO_BASE * F^alpha
  where alpha > 1

Biological idea:
- A keeps surface LDLR low (matches abundance phenotype)
- F boosts per-receptor uptake efficiency (rescues function without inflating abundance)

Run from ldlr_modeling/:
  python models/ldlr_model_a6.py p.Leu8Gln
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
F_MIN, F_MAX = 0.0, 1.5
A_MIN, A_MAX = 0.0, 1.5
EPS = 0.01

F_SNAP_LO, F_SNAP_HI = 0.90, 1.10
SKIP_SIM_FOR_CLUSTER1 = False
NULL_THRESHOLD = 0.02

# =========================
# WT base parameters
# =========================
LDLR_INIT_BASE = 1000
K_ENDO_BASE = 2.0
K_OFF_ENDO_BASE = 50.0

# separated trafficking / degradation bases
K_RECYCLE_LDLR_BASE = 3.0
K_LYSO_LDL_BASE = 3.0
K_DEGRADE_LDL_BASE = 5.0

# WT steady-state references
WT_SURFACE_SS = 587.858
WT_UPTAKE_SS = 621.993

# =========================
# METHOD 3B KEY PARAMETER
# Cluster 2 only:
# k_endo = K_ENDO_BASE * F^alpha
# =========================
CLUSTER2_F_ENDO_ALPHA = 1.5


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
    Returns:
        F_used, A_used, note, skip_sim, cluster
    """
    F_eff = clip_and_floor(F_raw, lo=F_MIN, hi=F_MAX, eps=EPS)
    A_eff = clip_and_floor(A_raw, lo=A_MIN, hi=A_MAX, eps=EPS)
    skip_sim = False

    # Cluster 0: scaled WT
    if cluster == 0:
        A_used = A_eff
        if F_SNAP_LO <= F_eff <= F_SNAP_HI:
            F_used = 1.0
            note = f"Cluster0 (scaled WT): A_used={A_used:.3f}, F snapped to 1.0."
        else:
            F_used = F_eff
            note = f"Cluster0 (scaled WT): A_used={A_used:.3f}, F_used={F_used:.3f}."
        return F_used, A_used, note, skip_sim, cluster

    # Cluster 2: METHOD 3B
    if cluster == 2:
        F_used = F_eff
        A_used = A_eff
        note = (
            f"Cluster2 (method 3B): A_used={A_used:.3f}, F_used={F_used:.3f}. "
            f"No abundance rescue; cluster-2 k_endo scales as F^{CLUSTER2_F_ENDO_ALPHA:.2f}."
        )
        return F_used, A_used, note, skip_sim, cluster

    # Cluster 3: functional-defect
    if cluster == 3:
        F_used = F_eff
        A_used = 1.0
        note = (
            f"Cluster3 (functional-defect): F_used={F_used:.3f} "
            f"scales uptake-side parameters; A_used=1.0."
        )
        return F_used, A_used, note, skip_sim, cluster

    # Cluster 1: near-null
    if cluster == 1:
        F_used = EPS
        A_used = EPS
        note = f"Cluster1 (near-null): set F_used=A_used=eps={EPS}."
        if SKIP_SIM_FOR_CLUSTER1 and (F_eff <= NULL_THRESHOLD and A_eff <= NULL_THRESHOLD):
            skip_sim = True
            note += f" Skipping simulation because F_eff and A_eff <= {NULL_THRESHOLD}."
        return F_used, A_used, note, skip_sim, cluster

    return F_eff, A_eff, f"Unknown cluster {cluster}: using clipped scores.", skip_sim, cluster


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
    ):
        self.variant_name = variant_name
        self.functional_score = float(functional_score)
        self.abundance_score = float(abundance_score)
        self.cluster = int(cluster)
        self.LDLR_init_base = int(LDLR_init_base)

        self.wt_surface_ss = float(wt_surface_ss)
        self.wt_uptake_ss = float(wt_uptake_ss)

        self.result = None
        self.expected_gdat_key = None
        self.run_dir = None

    def get_effective_abundance(self) -> float:
        # No abundance rescue in method 3B
        return max(self.abundance_score, EPS)

    def get_cluster2_kendo_scale(self, F: float) -> float:
        return F ** CLUSTER2_F_ENDO_ALPHA

    def get_scaled_parameters(self):
        """
        A -> receptor-side abundance / surface availability:
             - LDLR_init
             - k_recycle_ldlr

        F -> uptake-side performance:
             - k_endo
             - k_off_surf

        Cluster 2 METHOD 3B:
             k_endo = K_ENDO_BASE * F^alpha
             instead of just K_ENDO_BASE * F
        """
        F = max(self.functional_score, EPS)
        A_eff_model = self.get_effective_abundance()

        # receptor-side abundance parameters
        LDLR_init_scaled = int(round(self.LDLR_init_base * A_eff_model))
        k_recycle_ldlr_scaled = K_RECYCLE_LDLR_BASE * A_eff_model

        # LDL cargo parameters kept constant
        k_lyso_ldl_scaled = K_LYSO_LDL_BASE
        k_degrade_ldl_scaled = K_DEGRADE_LDL_BASE

        # uptake/function-side parameters
        if self.cluster == 2:
            k_endo_scaled = K_ENDO_BASE * self.get_cluster2_kendo_scale(F)
        else:
            k_endo_scaled = K_ENDO_BASE * F

        k_off_surf_scaled = 1.0 / F

        return {
            "A_eff_model": A_eff_model,
            "LDLR_init": LDLR_init_scaled,
            "k_recycle_ldlr": k_recycle_ldlr_scaled,
            "k_lyso_ldl": k_lyso_ldl_scaled,
            "k_degrade_ldl": k_degrade_ldl_scaled,
            "k_endo": k_endo_scaled,
            "k_off_surf": k_off_surf_scaled,
            "k_off_endo": K_OFF_ENDO_BASE,
        }

    def get_model_string(self):
        p = self.get_scaled_parameters()

        return f"""
# LDLR Model - {self.variant_name}
# Cluster: {self.cluster}
# F_used: {self.functional_score:.3f}
# A_used_raw: {self.abundance_score:.3f}
# A_used_model: {p["A_eff_model"]:.3f}

begin model

begin parameters
    # Binding / uptake-side
    k_on_base          1.0
    k_off_surf         {p["k_off_surf"]:.4f}
    k_off_endo         {p["k_off_endo"]:.4f}
    k_endo             {p["k_endo"]:.4f}

    # Module strengths
    strength_LA3       1.0
    strength_LA4       1.0
    strength_LA5       1.0
    strength_LA7       0.8

    # Separated trafficking / degradation parameters
    k_recycle_ldlr     {p["k_recycle_ldlr"]:.4f}
    k_lyso_ldl         {p["k_lyso_ldl"]:.4f}
    k_degrade_ldl      {p["k_degrade_ldl"]:.4f}

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

        run_dir = out_root / self.variant_name
        run_dir.mkdir(parents=True, exist_ok=True)
        self.run_dir = run_dir

        temp_model_file = run_dir / f"{self.variant_name}_temp.bngl"
        self.expected_gdat_key = temp_model_file.stem

        with open(temp_model_file, "w") as f:
            f.write(self.get_model_string())
            f.write("\ngenerate_network({overwrite=>1})\n")
            f.write(f"simulate({{method=>\"ode\", t_end=>{t_end}, n_steps=>{n_steps}}})\n")

        print(f"[DEBUG] BNGL written to: {temp_model_file}")

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

        fig, axes = plt.subplots(2, 2, figsize=(13, 10))

        axes[0, 0].plot(data["time"], data["LDLR_surface"], label="Surface LDLR", linewidth=2)
        axes[0, 0].plot(data["time"], data["LDLR_endosome"], label="Endosome LDLR", linewidth=2)
        axes[0, 0].axhline(self.wt_surface_ss, linestyle="--", linewidth=1, label=f"WT ss surface={self.wt_surface_ss:.3f}")
        axes[0, 0].set_title(
            f"Surface LDLR\nfinal/WTss = {metrics['surface_ratio_vs_WTss']:.3f}, target A = {metrics['target_A']:.3f}"
        )
        axes[0, 0].legend()
        axes[0, 0].grid(alpha=0.3)

        axes[0, 1].plot(data["time"], data["LDL_free"], label="Free LDL", linewidth=2)
        axes[0, 1].plot(data["time"], data["LDL_internalized"], label="Internalized LDL", linewidth=2, color="red")
        axes[0, 1].set_title("LDL Dynamics")
        axes[0, 1].legend()
        axes[0, 1].grid(alpha=0.3)

        axes[1, 0].plot(data["time"], data["LDL_internalized"], linewidth=2, color="darkred")
        axes[1, 0].axhline(self.wt_uptake_ss, linestyle="--", linewidth=1, label=f"WT ss uptake={self.wt_uptake_ss:.3f}")
        axes[1, 0].set_title(
            f"Steady-state LDL Uptake\nfinal/WTss = {metrics['uptake_ratio_vs_WTss']:.3f}, target F = {metrics['target_F']:.3f}"
        )
        axes[1, 0].legend()
        axes[1, 0].grid(alpha=0.3)

        axes[1, 1].plot(data["time"], data["LDLR_surface_complex"], linewidth=2, color="purple")
        axes[1, 1].set_title("LDLR-LDL Complexes at Surface")
        axes[1, 1].grid(alpha=0.3)

        fig.suptitle(
            f"{self.variant_name} (cluster={self.cluster}, F_used={self.functional_score:.2f}, A_used={self.abundance_score:.2f})",
            fontsize=14,
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
    ap.add_argument("--out_dir", default="results/data/results",
                    help="Root output directory; a subfolder per variant will be created.")
    ap.add_argument("--save_fig", default=None)
    ap.add_argument("--LDLR_init_base", type=int, default=1000)
    ap.add_argument("--wt_surface_ss", type=float, default=WT_SURFACE_SS)
    ap.add_argument("--wt_uptake_ss", type=float, default=WT_UPTAKE_SS)
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

    model = LDLRModel(
        variant_name=info["variant"],
        functional_score=F_used,
        abundance_score=A_used,
        cluster=cluster,
        LDLR_init_base=args.LDLR_init_base,
        wt_surface_ss=args.wt_surface_ss,
        wt_uptake_ss=args.wt_uptake_ss,
    )

    model.run(t_end=args.t_end, n_steps=args.n_steps, out_root=args.out_dir)

    if args.save_fig is None:
        args.save_fig = str(Path(args.out_dir) / info["variant"] / f"{info['variant']}.png")

    params = model.get_scaled_parameters()
    print("=== Scaled parameters used ===")
    print(f"A_eff_model       = {params['A_eff_model']:.4f}")
    print(f"LDLR_init         = {params['LDLR_init']}")
    print(f"k_endo            = {params['k_endo']:.4f}")
    print(f"k_off_surf        = {params['k_off_surf']:.4f}")
    print(f"k_recycle_ldlr    = {params['k_recycle_ldlr']:.4f}")
    print(f"k_lyso_ldl        = {params['k_lyso_ldl']:.4f}")
    print(f"k_degrade_ldl     = {params['k_degrade_ldl']:.4f}")

    if cluster == 2:
        print(f"cluster2 k_endo exponent alpha = {CLUSTER2_F_ENDO_ALPHA:.2f}")

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