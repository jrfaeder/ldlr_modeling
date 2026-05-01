"""
LDLR Model for PyPESTO Parameter Fitting
Based on ldlr_model_a7.py

This module provides a simplified interface for parameter fitting using pypesto.
The model can fit parameters to match target functional scores.
"""

import os
from contextlib import contextmanager
from pathlib import Path
import numpy as np
import pandas as pd
import bionetgen


# =========================
# Constants
# =========================
EPS = 0.01

# WT base parameters (these can be fitted)
LDLR_INIT_BASE = 1000
K_ENDO_BASE = 2.0
K_OFF_ENDO_BASE = 50.0
K_RECYCLE_LDLR_BASE = 3.0
K_LYSO_LDL_BASE = 3.0
K_DEGRADE_LDL_BASE = 5.0

# WT steady-state references
WT_SURFACE_SS = 587.858
WT_UPTAKE_SS = 621.993

# Cluster 2 specific parameter
CLUSTER2_A_RECYCLE_BETA = 0.1


@contextmanager
def pushd(new_dir: Path):
    """Context manager to temporarily change directory."""
    prev = Path.cwd()
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(prev)


class LDLRModelForFitting:
    """
    LDLR model wrapper for pypesto parameter fitting.

    This class allows fitting of model parameters to match target abundance and functional scores.
    Unlike the original model, here we fit the BASE parameters directly, not the scores.
    The scores are outputs from the model simulation.
    """

    def __init__(
        self,
        variant_name="WT",
        cluster=0,
        work_dir=None,
    ):
        """
        Initialize the LDLR model.

        Parameters
        ----------
        variant_name : str
            Name of the variant
        cluster : int
            Cluster assignment (0, 1, 2, or 3)
        work_dir : str or Path, optional
            Working directory for simulation files
        """
        self.variant_name = variant_name
        self.cluster = int(cluster)

        if work_dir is None:
            work_dir = Path("results/data/pypesto_fitting")
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

        self.result = None
        self.last_simulation_data = None

        # Cache for efficient simulation
        self.simulator = None
        self._model_loaded = False
        self.result_array = None

    def get_model_parameters(self, params_dict):
        """
        Get model parameters directly from the parameter dictionary.

        For fitting, we don't use functional/abundance scores as inputs.
        Instead, we use the base parameters directly and the scores become outputs.

        Parameters
        ----------
        params_dict : dict
            Dictionary containing base parameters:
            - LDLR_init_base
            - k_endo_base
            - k_off_endo_base
            - k_recycle_ldlr_base
            - k_lyso_ldl_base
            - k_degrade_ldl_base
            - k_on_base (optional, default 1.0)

        Returns
        -------
        dict
            Model parameters for the simulation
        """
        # Extract base parameters from input
        ldlr_init_base = params_dict.get('LDLR_init_base', LDLR_INIT_BASE)
        k_endo_base = params_dict.get('k_endo_base', K_ENDO_BASE)
        k_off_endo_base = params_dict.get('k_off_endo_base', K_OFF_ENDO_BASE)
        k_recycle_ldlr_base = params_dict.get('k_recycle_ldlr_base', K_RECYCLE_LDLR_BASE)
        k_lyso_ldl_base = params_dict.get('k_lyso_ldl_base', K_LYSO_LDL_BASE)
        k_degrade_ldl_base = params_dict.get('k_degrade_ldl_base', K_DEGRADE_LDL_BASE)
        k_on_base = params_dict.get('k_on_base', 1.0)

        # For k_off_surf, we use the inverse relationship with k_on
        # This maintains binding equilibrium behavior
        k_off_surf = params_dict.get('k_off_surf', 1.0)

        return {
            "LDLR_init": int(round(ldlr_init_base)),
            "k_on_base": k_on_base,
            "k_off_surf": k_off_surf,
            "k_endo": k_endo_base,
            "k_off_endo": k_off_endo_base,
            "k_recycle_ldlr": k_recycle_ldlr_base,
            "k_lyso_ldl": k_lyso_ldl_base,
            "k_degrade_ldl": k_degrade_ldl_base,
        }

    def get_model_string(self, model_params):
        """
        Generate the BNGL model string.

        Parameters
        ----------
        model_params : dict
            Model parameters from get_model_parameters()

        Returns
        -------
        str
            BNGL model string
        """
        p = model_params

        return f"""
# LDLR Model - {self.variant_name}
# Cluster: {self.cluster}
# Parameter fitting mode: base parameters only

begin model

begin parameters
    # Binding / uptake-side
    k_on_base          {p["k_on_base"]:.6f}
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

    def _initialize_simulator(self, params_dict):
        """
        Initialize the BioNetGen model and simulator (called once).

        Parameters
        ----------
        params_dict : dict
            Dictionary of base parameters
        """
        # Get model parameters for this variant
        model_params = self.get_model_parameters(params_dict)

        # Generate the BNGL model string
        model_str = self.get_model_string(model_params)

        # Write BNGL file
        bngl_filename = f"{self.variant_name}_fitting.bngl"
        bngl_path = self.work_dir / bngl_filename
        with open(bngl_path, "w") as f:
            f.write(model_str)

        # Load model and create simulator
        # Change to work_dir so bionetgen can find the file
        with pushd(self.work_dir):
            self.model = bionetgen.bngmodel(bngl_filename)
            self.simulator = self.model.setup_simulator()

        self._model_loaded = True

    def simulate(self, params_dict, t_end=200, n_steps=1000):
        """
        Run simulation with given parameters.

        Parameters
        ----------
        params_dict : dict
            Dictionary of base parameters
        t_end : float
            End time for simulation
        n_steps : int
            Number of time steps

        Returns
        -------
        dict
            Simulation metrics including final_uptake, final_surface, etc.
        """
        # Initialize simulator on first call
        if not self._model_loaded:
            self._initialize_simulator(params_dict)

        # Get model parameters
        model_params = self.get_model_parameters(params_dict)

        # Set parameters in simulator
        for param_name, value in model_params.items():
            if param_name == "LDLR_init":
                # Skip initial conditions - handled separately
                continue
            try:
                self.simulator[param_name] = value
            except:
                pass  # Parameter might not exist or be settable

        # Reset simulator
        self.simulator.reset()

        # Run simulation
        observables = [o for o in self.model.observables]
        selections = ['time'] + observables
        result_array = self.simulator.simulate(0, t_end, n_steps, selections=selections)

        # Store result array for extraction
        self.result_array = result_array

        # Extract results
        return self._extract_metrics()

    def _extract_metrics(self):
        """Extract metrics from simulation results."""
        if self.result_array is None:
            raise ValueError("Must run simulation first")

        # Get observables list from model
        observables = [o for o in self.model.observables]

        # Create column names: ['time'] + observables
        columns = ['time'] + observables

        # Convert result array to DataFrame
        df = pd.DataFrame(self.result_array, columns=columns)

        # Calculate internalized LDL (add as column for plotting)
        # LDL_endo and LDL_lyso are observables in the model
        df["LDL_internalized"] = df["LDL_endo"] + df["LDL_lyso"]

        # Get final values
        final_endo = float(df["LDL_endo"].iloc[-1])
        final_lyso = float(df["LDL_lyso"].iloc[-1])
        final_uptake = final_endo + final_lyso

        # Get final surface value
        final_surface = float(df["LDLR_surface"].iloc[-1])

        # Calculate ratios vs WT
        surface_ratio = final_surface / WT_SURFACE_SS
        uptake_ratio = final_uptake / WT_UPTAKE_SS

        self.last_simulation_data = df

        return {
            "final_surface": final_surface,
            "final_uptake": final_uptake,
            "surface_ratio": surface_ratio,
            "uptake_ratio": uptake_ratio,
        }


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
    Create a pypesto-compatible objective function.

    This function fits model BASE parameters to match BOTH target functional and abundance scores.
    The scores are OUTPUTS from the model (surface_ratio and uptake_ratio), not inputs.

    Parameters
    ----------
    variant_name : str
        Name of the variant
    target_functional_score : float
        Target functional score to fit (uptake_ratio relative to WT)
    target_abundance_score : float
        Target abundance score to fit (surface_ratio relative to WT)
    cluster : int
        Cluster assignment (0, 1, 2, or 3)
    work_dir : str or Path, optional
        Working directory for simulations
    weight_functional : float, optional
        Weight for functional score in objective (default 1.0)
    weight_abundance : float, optional
        Weight for abundance score in objective (default 1.0)

    Returns
    -------
    tuple
        (objective_function, model) where objective_function takes parameter vector
        and returns weighted sum of squared errors
    """
    model = LDLRModelForFitting(
        variant_name=variant_name,
        cluster=cluster,
        work_dir=work_dir,
    )

    def objective(params_vector):
        """
        Objective function for pypesto.

        Parameters
        ----------
        params_vector : array-like
            Parameter vector with order:
            [LDLR_init_base, k_endo_base, k_off_endo_base,
             k_recycle_ldlr_base, k_lyso_ldl_base, k_degrade_ldl_base, k_off_surf]

        Returns
        -------
        float
            Weighted sum of squared errors for both functional and abundance scores
        """
        # Create parameter dictionary
        param_names = [
            'LDLR_init_base',
            'k_endo_base',
            'k_off_endo_base',
            'k_recycle_ldlr_base',
            'k_lyso_ldl_base',
            'k_degrade_ldl_base',
            'k_off_surf',
        ]

        params_dict = dict(zip(param_names, params_vector))

        try:
            # Run simulation
            metrics = model.simulate(params_dict)

            # Calculate errors for both scores
            # Functional score = uptake_ratio (final uptake / WT uptake)
            predicted_functional = metrics['uptake_ratio']
            error_functional = weight_functional * (predicted_functional - target_functional_score) ** 2

            # Abundance score = surface_ratio (final surface / WT surface)
            predicted_abundance = metrics['surface_ratio']
            error_abundance = weight_abundance * (predicted_abundance - target_abundance_score) ** 2

            # Total error (weighted sum of squared errors)
            total_error = error_functional + error_abundance

            return total_error

        except Exception as e:
            # Return large penalty if simulation fails
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
    Create a pypesto-compatible objective function that operates in LOG SPACE.

    This improves optimization for parameters spanning multiple orders of magnitude.
    The optimizer works with log10(parameters), which are then transformed back
    to linear scale before simulation.

    This function fits model BASE parameters to match BOTH target functional and abundance scores.
    The scores are OUTPUTS from the model (surface_ratio and uptake_ratio), not inputs.

    Parameters
    ----------
    variant_name : str
        Name of the variant
    target_functional_score : float
        Target functional score to fit (uptake_ratio relative to WT)
    target_abundance_score : float
        Target abundance score to fit (surface_ratio relative to WT)
    cluster : int
        Cluster assignment (0, 1, 2, or 3)
    work_dir : str or Path, optional
        Working directory for simulations
    weight_functional : float, optional
        Weight for functional score in objective (default 1.0)
    weight_abundance : float, optional
        Weight for abundance score in objective (default 1.0)

    Returns
    -------
    tuple
        (objective_function, model) where objective_function takes LOG-SCALE parameter vector
        and returns weighted sum of squared errors
    """
    model = LDLRModelForFitting(
        variant_name=variant_name,
        cluster=cluster,
        work_dir=work_dir,
    )

    def objective_log(log_params_vector):
        """
        Objective function for pypesto operating in log space.

        Parameters
        ----------
        log_params_vector : array-like
            LOG10-scaled parameter vector with order:
            [log10(LDLR_init_base), log10(k_endo_base), log10(k_off_endo_base),
             log10(k_recycle_ldlr_base), log10(k_lyso_ldl_base),
             log10(k_degrade_ldl_base), log10(k_off_surf)]

        Returns
        -------
        float
            Weighted sum of squared errors for both functional and abundance scores
        """
        # Transform from log space to linear space
        params_vector = 10 ** log_params_vector

        # Create parameter dictionary
        param_names = [
            'LDLR_init_base',
            'k_endo_base',
            'k_off_endo_base',
            'k_recycle_ldlr_base',
            'k_lyso_ldl_base',
            'k_degrade_ldl_base',
            'k_off_surf',
        ]

        params_dict = dict(zip(param_names, params_vector))

        try:
            # Run simulation
            metrics = model.simulate(params_dict)

            # Calculate errors for both scores
            # Functional score = uptake_ratio (final uptake / WT uptake)
            predicted_functional = metrics['uptake_ratio']
            error_functional = weight_functional * (predicted_functional - target_functional_score) ** 2

            # Abundance score = surface_ratio (final surface / WT surface)
            predicted_abundance = metrics['surface_ratio']
            error_abundance = weight_abundance * (predicted_abundance - target_abundance_score) ** 2

            # Total error (weighted sum of squared errors)
            total_error = error_functional + error_abundance

            return total_error

        except Exception as e:
            # Return large penalty if simulation fails
            print(f"Simulation failed: {e}")
            return 1e10

    return objective_log, model
