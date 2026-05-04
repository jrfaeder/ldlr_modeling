"""
LDLR Model for PyPESTO Parameter Fitting
A8-based version: wraps models/ldlr_model_a8.py.

Purpose
-------
This module provides a PyPESTO-compatible objective for fitting LDLR model
parameters to target abundance and functional scores.

Conceptual mapping
------------------
    k -> LDLRModel(A8) -> (A_pred, F_pred)

where:
    A_pred = surface_ratio_vs_WTss
    F_pred = uptake_ratio_vs_WTss

Important
---------
This wrapper treats A and F as OUTPUTS, not inputs. The model is initialized
with WT-like scores (F=1, A=1), and candidate parameters are applied through
A8 direct parameter overrides.

Expected location in your repo
------------------------------
    models/ldlr_pypesto_model2.py

and A8 should remain here:
    models/ldlr_model_a8.py
"""

from pathlib import Path
import numpy as np

# ---------------------------------------------------------------------
# Import A8 model and constants
# ---------------------------------------------------------------------
# This works when running notebooks/scripts from project root.
# The fallback works if this file is executed/imported from inside models/.
try:
    from models.ldlr_model_a8 import (
        LDLRModel,
        LDLR_INIT_BASE,
        K_ENDO_BASE,
        K_OFF_ENDO_BASE,
        K_RECYCLE_LDLR_BASE,
        K_LYSO_LDL_BASE,
        K_DEGRADE_LDL_BASE,
        WT_SURFACE_SS,
        WT_UPTAKE_SS,
    )
except ImportError:
    from ldlr_model_a8 import (
        LDLRModel,
        LDLR_INIT_BASE,
        K_ENDO_BASE,
        K_OFF_ENDO_BASE,
        K_RECYCLE_LDLR_BASE,
        K_LYSO_LDL_BASE,
        K_DEGRADE_LDL_BASE,
        WT_SURFACE_SS,
        WT_UPTAKE_SS,
    )


EPS = 0.01

# Parameter order expected by the PyPESTO notebook.
PARAM_NAMES = [
    "LDLR_init_base",
    "k_endo_base",
    "k_off_endo_base",
    "k_recycle_ldlr_base",
    "k_lyso_ldl_base",
    "k_degrade_ldl_base",
    "k_off_surf",
]


class LDLRModelForFitting:
    """
    PyPESTO-compatible wrapper around ldlr_model_a8.LDLRModel.

    The fitting parameters are interpreted as candidate model parameters k.
    For each candidate k, the wrapper:
      1. builds the A8 model with WT scores: F=1, A=1
      2. applies k using A8 param_overrides
      3. runs the BioNetGen ODE simulation
      4. returns A_pred and F_pred ratios

    This keeps the PyPESTO fitting workflow consistent with the recent A8-based
    reachability, inverse-search, and surrogate analyses.
    """

    def __init__(
        self,
        variant_name="WT",
        cluster=0,
        work_dir=None,
        t_end=200,
        n_steps=1000,
    ):
        self.variant_name = str(variant_name)
        self.cluster = int(cluster)

        if work_dir is None:
            work_dir = Path("results/data/pypesto_fitting_a8")
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

        self.t_end = t_end
        self.n_steps = n_steps

        self.last_simulation_data = None
        self.last_metrics = None
        self.last_params = None
        self.last_model = None
        self.eval_count = 0

    def get_model_parameters(self, params_dict):
        """
        Convert PyPESTO fitting parameter names into A8 direct overrides.

        Input parameter names are kept compatible with the original notebook:
            LDLR_init_base
            k_endo_base
            k_off_endo_base
            k_recycle_ldlr_base
            k_lyso_ldl_base
            k_degrade_ldl_base
            k_off_surf

        A8 override names are:
            LDLR_init
            k_endo
            k_off_endo
            k_recycle_ldlr
            k_lyso_ldl
            k_degrade_ldl
            k_off_surf
        """
        ldlr_init = params_dict.get("LDLR_init_base", LDLR_INIT_BASE)
        k_endo = params_dict.get("k_endo_base", K_ENDO_BASE)
        k_off_endo = params_dict.get("k_off_endo_base", K_OFF_ENDO_BASE)
        k_recycle_ldlr = params_dict.get("k_recycle_ldlr_base", K_RECYCLE_LDLR_BASE)
        k_lyso_ldl = params_dict.get("k_lyso_ldl_base", K_LYSO_LDL_BASE)
        k_degrade_ldl = params_dict.get("k_degrade_ldl_base", K_DEGRADE_LDL_BASE)
        k_off_surf = params_dict.get("k_off_surf", 1.0)

        return {
            "LDLR_init": max(1, int(round(float(ldlr_init)))),
            "k_endo": max(float(k_endo), 1e-9),
            "k_off_endo": max(float(k_off_endo), 1e-9),
            "k_recycle_ldlr": max(float(k_recycle_ldlr), 1e-9),
            "k_lyso_ldl": max(float(k_lyso_ldl), 1e-9),
            "k_degrade_ldl": max(float(k_degrade_ldl), 1e-9),
            "k_off_surf": max(float(k_off_surf), 1e-9),
        }

    def simulate(self, params_dict, t_end=None, n_steps=None):
        """
        Run A8 model with direct parameter overrides.

        Returns
        -------
        dict
            final_surface, final_uptake, surface_ratio, uptake_ratio
        """
        if t_end is None:
            t_end = self.t_end
        if n_steps is None:
            n_steps = self.n_steps

        overrides = self.get_model_parameters(params_dict)
        self.eval_count += 1

        # IMPORTANT:
        # We use WT-like scores so A/F are not inputs. Candidate k is applied
        # via param_overrides, and A/F are read from the simulation outputs.
        model = LDLRModel(
            variant_name=self.variant_name,
            functional_score=1.0,
            abundance_score=1.0,
            cluster=0,
            LDLR_init_base=LDLR_INIT_BASE,
            wt_surface_ss=WT_SURFACE_SS,
            wt_uptake_ss=WT_UPTAKE_SS,
            param_overrides=overrides,
        )

        safe_variant = self.variant_name.replace("/", "_")
        out_root = self.work_dir / safe_variant

        model.run(t_end=t_end, n_steps=n_steps, out_root=out_root)
        metrics = model.get_summary_metrics()

        # Store useful info for notebook inspection/plotting.
        self.last_model = model
        self.last_metrics = metrics
        self.last_params = overrides
        try:
            self.last_simulation_data = model.get_data()
        except Exception:
            self.last_simulation_data = None

        return {
            "final_surface": float(metrics["final_surface"]),
            "final_uptake": float(metrics["final_uptake"]),
            "surface_ratio": float(metrics["surface_ratio_vs_WTss"]),
            "uptake_ratio": float(metrics["uptake_ratio_vs_WTss"]),
        }


# ---------------------------------------------------------------------
# Objective utilities
# ---------------------------------------------------------------------

def _vector_to_params_dict(params_vector):
    params_vector = np.asarray(params_vector, dtype=float)
    return dict(zip(PARAM_NAMES, params_vector))


def _loss_from_metrics(
    metrics,
    target_functional_score,
    target_abundance_score,
    weight_functional=1.0,
    weight_abundance=1.0,
):
    """
    PyPESTO objective loss.

    Functional score target corresponds to uptake_ratio.
    Abundance score target corresponds to surface_ratio.
    """
    predicted_functional = float(metrics["uptake_ratio"])
    predicted_abundance = float(metrics["surface_ratio"])

    error_functional = weight_functional * (
        predicted_functional - target_functional_score
    ) ** 2

    error_abundance = weight_abundance * (
        predicted_abundance - target_abundance_score
    ) ** 2

    return float(error_functional + error_abundance)


def create_objective_function(
    variant_name,
    target_functional_score,
    target_abundance_score,
    cluster=0,
    work_dir=None,
    weight_functional=1.0,
    weight_abundance=1.0,
):
    """
    Create a PyPESTO-compatible objective function in linear parameter space.

    Parameter vector order:
        [LDLR_init_base, k_endo_base, k_off_endo_base,
         k_recycle_ldlr_base, k_lyso_ldl_base, k_degrade_ldl_base,
         k_off_surf]
    """
    model = LDLRModelForFitting(
        variant_name=variant_name,
        cluster=cluster,
        work_dir=work_dir,
    )

    def objective(params_vector):
        params_dict = _vector_to_params_dict(params_vector)

        try:
            metrics = model.simulate(params_dict)
            return _loss_from_metrics(
                metrics,
                target_functional_score=target_functional_score,
                target_abundance_score=target_abundance_score,
                weight_functional=weight_functional,
                weight_abundance=weight_abundance,
            )
        except Exception as e:
            print(f"Simulation failed: {e}")
            return 1e10

    return objective, model


def create_objective_function_log_scale(
    variant_name,
    target_functional_score,
    target_abundance_score,
    cluster=0,
    work_dir=None,
    weight_functional=1.0,
    weight_abundance=1.0,
):
    """
    Create a PyPESTO-compatible objective function in log10 parameter space.

    log_params_vector order:
        [log10(LDLR_init_base), log10(k_endo_base), log10(k_off_endo_base),
         log10(k_recycle_ldlr_base), log10(k_lyso_ldl_base),
         log10(k_degrade_ldl_base), log10(k_off_surf)]
    """
    model = LDLRModelForFitting(
        variant_name=variant_name,
        cluster=cluster,
        work_dir=work_dir,
    )

    def objective_log(log_params_vector):
        try:
            log_params_vector = np.asarray(log_params_vector, dtype=float)
            params_vector = 10 ** log_params_vector
            params_dict = _vector_to_params_dict(params_vector)
            metrics = model.simulate(params_dict)

            return _loss_from_metrics(
                metrics,
                target_functional_score=target_functional_score,
                target_abundance_score=target_abundance_score,
                weight_functional=weight_functional,
                weight_abundance=weight_abundance,
            )

        except Exception as e:
            print(f"Simulation failed: {e}")
            return 1e10

    return objective_log, model
