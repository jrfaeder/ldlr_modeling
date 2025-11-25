"""
LDLR model using PyBioNetGen
Simplified version with 4 LA modules
"""

import bionetgen
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

class LDLRModel:
    """LDLR-LDL uptake model"""
    
    def __init__(self, variant_name='WT', functional_score=1.0):
        """
        Initialize model
        
        Parameters:
        -----------
        variant_name : str
            Variant identifier (e.g., 'WT', 'C52Y')
        functional_score : float
            Functional score from experiment (0-1.5)
        """
        self.variant_name = variant_name
        self.functional_score = functional_score
        self.result = None
        
    def get_model_string(self):
        """Generate BioNetGen model string"""
        
        # Scale parameters by functional score
        binding_strength = max(self.functional_score, 0.01)
        
        model = f"""
# LDLR Model - {self.variant_name}
# Functional score: {self.functional_score:.3f}

begin model

begin parameters
    # Binding (scaled by functional score)
    k_on_base      1.0
    k_off_surf     {1.0 / binding_strength:.4f}
    k_off_endo     50.0
    
    # Module strengths
    strength_LA3   1.0
    strength_LA4   1.0
    strength_LA5   1.0
    strength_LA7   0.8
    
    # Trafficking
    k_endo         2.0
    k_recycle      3.0
    k_degrade      5.0
    
    # Initial conditions
    LDLR_init      1000
    LDL_conc       100
end parameters

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
    Molecules  LDL_internalized  LDL(loc~endo) LDL(loc~lyso)
    Molecules  Complex_surface   LDLR(la3!+,loc~surface)
end observables

begin reaction rules
    # LDLR-LDL Binding
    LDLR(la3,loc~surface) + LDL(ldlr,loc~extra) <-> \\
        LDLR(la3!1,loc~surface).LDL(ldlr!1,loc~extra) \\
        k_on_base*strength_LA3, k_off_surf

    LDLR(la4,loc~surface) + LDL(ldlr,loc~extra) <-> \\
        LDLR(la4!1,loc~surface).LDL(ldlr!1,loc~extra) \\
        k_on_base*strength_LA4, k_off_surf

    LDLR(la5,loc~surface) + LDL(ldlr,loc~extra) <-> \\
        LDLR(la5!1,loc~surface).LDL(ldlr!1,loc~extra) \\
        k_on_base*strength_LA5, k_off_surf

    LDLR(la7,loc~surface) + LDL(ldlr,loc~extra) <-> \\
        LDLR(la7!1,loc~surface).LDL(ldlr!1,loc~extra) \\
        k_on_base*strength_LA7, k_off_surf

    # Endocytosis
    LDLR(la3!1,loc~surface).LDL(ldlr!1,loc~extra) -> \\
        LDLR(la3!1,loc~endosome).LDL(ldlr!1,loc~endo)  k_endo

    LDLR(la4!1,loc~surface).LDL(ldlr!1,loc~extra) -> \\
        LDLR(la4!1,loc~endosome).LDL(ldlr!1,loc~endo)  k_endo

    LDLR(la5!1,loc~surface).LDL(ldlr!1,loc~extra) -> \\
        LDLR(la5!1,loc~endosome).LDL(ldlr!1,loc~endo)  k_endo

    LDLR(la7!1,loc~surface).LDL(ldlr!1,loc~extra) -> \\
        LDLR(la7!1,loc~endosome).LDL(ldlr!1,loc~endo)  k_endo

    # LDL Release (pH-dependent)
    LDLR(la3!1,loc~endosome).LDL(ldlr!1,loc~endo) -> \\
        LDLR(la3,loc~endosome) + LDL(ldlr,loc~endo)  k_off_endo

    LDLR(la4!1,loc~endosome).LDL(ldlr!1,loc~endo) -> \\
        LDLR(la4,loc~endosome) + LDL(ldlr,loc~endo)  k_off_endo

    LDLR(la5!1,loc~endosome).LDL(ldlr!1,loc~endo) -> \\
        LDLR(la5,loc~endosome) + LDL(ldlr,loc~endo)  k_off_endo

    LDLR(la7!1,loc~endosome).LDL(ldlr!1,loc~endo) -> \\
        LDLR(la7,loc~endosome) + LDL(ldlr,loc~endo)  k_off_endo

    # Receptor Recycling
    LDLR(la3,la4,la5,la7,loc~endosome) -> \\
        LDLR(la3,la4,la5,la7,loc~surface)  k_recycle

    # LDL Degradation
    LDL(ldlr,loc~endo) -> LDL(ldlr,loc~lyso)  k_recycle
    LDL(loc~lyso) -> 0  k_degrade

end reaction rules

end model
"""
        return model
    
    def run(self, t_end=200, n_steps=200, out_dir='results/data'):
        """Run the simulation"""

        model_str = self.get_model_string()

        # Create output directory if it doesn't exist
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        # Write model to temporary file
        temp_model_file = out_path / f"{self.variant_name}_temp.bngl"
        with open(temp_model_file, 'w') as f:
            f.write(model_str)
            # Add simulation action
            f.write(f"\ngenerate_network({{overwrite=>1}})\n")
            f.write(f"simulate({{method=>\"ode\", t_end=>{t_end}, n_steps=>{n_steps}}})\n")

        print(f"Running {self.variant_name}...")
        self.result = bionetgen.run(
            str(temp_model_file),
            out=str(out_path),
            suppress=True
        )

        return self.result
    
    def get_data(self):
        """Load simulation data"""
        if self.result is None:
            raise ValueError("Must run simulation first")

        # Get the gdat data - result.gdats is a dict with model name as key
        model_name = list(self.result.gdats.keys())[0]
        gdat_array = self.result.gdats[model_name]

        # Convert to pandas DataFrame
        return pd.DataFrame(gdat_array)
    
    def get_final_uptake(self):
        """Get final LDL internalized count"""
        data = self.get_data()
        return data['LDL_internalized'].iloc[-1]

# Test if run directly
if __name__ == "__main__":
    print("Testing LDLR model...")
    
    model = LDLRModel('WT', 1.0)
    model.run()
    
    data = model.get_data()
    print(f"\nFinal LDL internalized: {data['LDL_internalized'].iloc[-1]:.1f}")
    print("✓ Model test passed!")
