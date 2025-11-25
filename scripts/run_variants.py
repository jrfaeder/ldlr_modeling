"""
Batch runner for LDLR variants
Runs 5 variants and saves results
"""

import sys
sys.path.append('models')

import pandas as pd
from pathlib import Path
from ldlr_model import LDLRModel
from scipy.stats import pearsonr

# Create directories
Path('results/data').mkdir(parents=True, exist_ok=True)

# Define variants
VARIANTS = {
    'WT': {'score': 1.00, 'domain': 'WT', 'clinvar': 'WT'},
    'C52Y': {'score': 0.05, 'domain': 'LA1', 'clinvar': 'pathogenic'},
    'D147N': {'score': 0.15, 'domain': 'LA3', 'clinvar': 'pathogenic'},
    'P526L': {'score': 0.95, 'domain': 'EGF-A', 'clinvar': 'benign'},
    'T705I': {'score': 0.92, 'domain': 'beta_propeller', 'clinvar': 'benign'}
}

def main():
    print("="*60)
    print("LDLR VARIANT SIMULATION")
    print("="*60)
    print(f"\nRunning {len(VARIANTS)} variants...\n")
    
    results = []
    
    for i, (name, info) in enumerate(VARIANTS.items(), 1):
        print(f"[{i}/{len(VARIANTS)}] {name}...", end=' ')
        
        try:
            model = LDLRModel(name, info['score'])
            model.run()
            uptake = model.get_final_uptake()
            
            print(f"✓ (uptake: {uptake:.1f})")
            
            results.append({
                'variant': name,
                'exp_score': info['score'],
                'model_uptake': uptake,
                'domain': info['domain'],
                'clinvar': info['clinvar']
            })
        except Exception as e:
            print(f"❌ FAILED: {e}")
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Normalize to WT
    wt_uptake = df.loc[df['variant']=='WT', 'model_uptake'].values[0]
    df['model_score'] = df['model_uptake'] / wt_uptake
    
    # Save
    df.to_csv('results/simulation_results.csv', index=False)
    print(f"\n💾 Results saved")
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Variants: {len(df)}")
    print(f"WT uptake: {wt_uptake:.1f}")
    
    r, p = pearsonr(df['model_score'], df['exp_score'])
    print(f"\nCorrelation: r = {r:.3f}, p = {p:.3e}")
    
    if r > 0.6:
        print("✓ STRONG correlation")
    elif r > 0.4:
        print("~ MODERATE correlation")
    else:
        print("✗ WEAK correlation")
    
    print("\n✅ Complete! Next: python scripts/analyze_results.py")

if __name__ == "__main__":
    main()
