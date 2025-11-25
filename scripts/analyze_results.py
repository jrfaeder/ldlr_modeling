"""
Analyze LDLR variant results
Creates validation and separation plots
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from pathlib import Path

# Setup
FIGURES_DIR = Path('results/figures')
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

sns.set_style("whitegrid")

def load_results():
    """Load simulation results"""
    df = pd.read_csv('results/simulation_results.csv')
    return df

def plot_validation(df):
    """Create validation scatter plot"""
    fig, ax = plt.subplots(figsize=(8, 7))
    
    colors = {'pathogenic': 'red', 'benign': 'blue', 'WT': 'green'}
    
    for clinvar in df['clinvar'].unique():
        mask = df['clinvar'] == clinvar
        ax.scatter(
            df.loc[mask, 'exp_score'],
            df.loc[mask, 'model_score'],
            c=colors.get(clinvar, 'gray'),
            label=clinvar,
            s=150,
            alpha=0.7,
            edgecolors='black',
            linewidths=1
        )
    
    ax.plot([0, 1.5], [0, 1.5], 'k--', alpha=0.5, linewidth=2)
    
    r, p = pearsonr(df['model_score'], df['exp_score'])
    
    ax.text(
        0.05, 0.95,
        f"r = {r:.3f}\np = {p:.2e}",
        transform=ax.transAxes,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
        fontsize=12,
        fontweight='bold'
    )
    
    ax.set_xlabel('Experimental Score', fontsize=13, fontweight='bold')
    ax.set_ylabel('Model Score', fontsize=13, fontweight='bold')
    ax.set_title('Model Validation', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right')
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'validation.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved validation.png")
    
    return r, p

def plot_separation(df):
    """Plot pathogenic vs benign"""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    classified = df[df['clinvar'].isin(['pathogenic', 'benign'])]
    
    data_to_plot = [
        classified[classified['clinvar']=='pathogenic']['model_score'],
        classified[classified['clinvar']=='benign']['model_score']
    ]
    
    bp = ax.boxplot(
        data_to_plot,
        labels=['Pathogenic', 'Benign'],
        patch_artist=True,
        showmeans=True
    )
    
    bp['boxes'][0].set_facecolor('lightcoral')
    bp['boxes'][1].set_facecolor('lightblue')
    
    for i, data in enumerate(data_to_plot, 1):
        y = data
        x = np.random.normal(i, 0.04, size=len(y))
        ax.scatter(x, y, alpha=0.6, s=100, edgecolors='black')
    
    ax.set_ylabel('Model Score', fontsize=13, fontweight='bold')
    ax.set_title('Pathogenic vs. Benign Separation', fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'separation.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved separation.png")

def print_summary(df, r, p):
    """Print results summary"""
    print("\n" + "="*60)
    print("RESULTS SUMMARY")
    print("="*60)
    
    print(f"\nCorrelation: r = {r:.4f}, p = {p:.2e}")
    
    path_mean = df[df['clinvar']=='pathogenic']['model_score'].mean()
    benign_mean = df[df['clinvar']=='benign']['model_score'].mean()
    
    print(f"\nPathogenic mean: {path_mean:.3f}")
    print(f"Benign mean: {benign_mean:.3f}")
    print(f"Separation: {benign_mean - path_mean:.3f}")
    
    mae = np.abs(df['model_score'] - df['exp_score']).mean()
    print(f"\nMean Absolute Error: {mae:.4f}")
    
    print("\n" + "="*60)

def main():
    print("="*60)
    print("LDLR VARIANT ANALYSIS")
    print("="*60)
    
    df = load_results()
    print(f"\nLoaded {len(df)} variants\n")
    
    print("Generating plots...")
    r, p = plot_validation(df)
    plot_separation(df)
    
    print_summary(df, r, p)
    
    print(f"\n✅ Complete! Figures in {FIGURES_DIR}")

if __name__ == "__main__":
    main()
