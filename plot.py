#!/usr/bin/env python3
"""
NexRig Harmonic Analysis Plotter
Generates publication-quality plots from harmonic_analysis.csv
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def load_results(filename='harmonic_analysis.csv'):
    """Load analysis results from CSV file"""
    if not Path(filename).exists():
        print(f"Error: {filename} not found. Run ./harmonic_analyzer first.")
        return None
    
    df = pd.read_csv(filename)
    return df

def plot_harmonic_spectrum(df, output_dir='plots'):
    """Plot harmonic content for each band (nominal corner only)"""
    Path(output_dir).mkdir(exist_ok=True)
    
    # Filter for nominal corner only
    df_nom = df[df['Corner'] == 'Nominal'].copy()
    
    bands = df_nom['Band'].unique()
    
    fig, axes = plt.subplots(2, 4, figsize=(16, 10))
    axes = axes.flatten()
    
    for idx, band in enumerate(bands):
        ax = axes[idx]
        band_data = df_nom[df_nom['Band'] == band]
        
        # Get bottom and top of band
        freq_low = band_data.iloc[0]
        freq_high = band_data.iloc[1] if len(band_data) > 1 else freq_low
        
        # Plot harmonics for bottom of band
        harmonics = range(2, 11)
        h_cols = [f'H{n}_dBc' for n in harmonics]
        
        values_low = [freq_low[col] for col in h_cols]
        values_high = [freq_high[col] for col in h_cols]
        
        ax.plot(harmonics, values_low, 'b-o', label=f'{freq_low["Frequency_MHz"]:.2f} MHz', linewidth=2)
        ax.plot(harmonics, values_high, 'r-s', label=f'{freq_high["Frequency_MHz"]:.2f} MHz', linewidth=2)
        
        # Add FCC requirement line
        ax.axhline(-43, color='g', linestyle='--', linewidth=1.5, label='FCC Limit')
        
        ax.set_xlabel('Harmonic Number', fontsize=10)
        ax.set_ylabel('Level (dBc)', fontsize=10)
        ax.set_title(f'{band} Band', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        ax.set_ylim(-100, -40)
        ax.set_xticks(harmonics)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/harmonic_spectrum_all_bands.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/harmonic_spectrum_all_bands.png")
    plt.close()

def plot_insertion_loss(df, output_dir='plots'):
    """Plot insertion loss vs frequency for all corners"""
    Path(output_dir).mkdir(exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    corners = df['Corner'].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(corners)))
    
    for corner, color in zip(corners, colors):
        corner_data = df[df['Corner'] == corner].sort_values('Frequency_MHz')
        
        # Plot with markers
        marker = 'o' if corner == 'Nominal' else 's'
        linewidth = 2.5 if corner == 'Nominal' else 1.5
        
        ax.plot(corner_data['Frequency_MHz'], 
                corner_data['Insertion_Loss_dB'],
                marker=marker, 
                linewidth=linewidth,
                label=corner,
                color=color,
                markersize=6)
    
    ax.set_xlabel('Frequency (MHz)', fontsize=12)
    ax.set_ylabel('Insertion Loss (dB)', fontsize=12)
    ax.set_title('LPF Insertion Loss vs Frequency (All Tolerance Corners)', 
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best')
    
    # Format y-axis to show losses as positive numbers
    ax.invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/insertion_loss_all_corners.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/insertion_loss_all_corners.png")
    plt.close()

def plot_component_stress(df, output_dir='plots'):
    """Plot component current stress for nominal corner"""
    Path(output_dir).mkdir(exist_ok=True)
    
    df_nom = df[df['Corner'] == 'Nominal'].copy()
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Current stress
    freq = df_nom['Frequency_MHz']
    
    ax1.plot(freq, df_nom['C1_Irms_A'], 'b-o', label='C1', linewidth=2, markersize=6)
    ax1.plot(freq, df_nom['C2_Irms_A'], 'r-s', label='C2', linewidth=2, markersize=6)
    ax1.plot(freq, df_nom['C3_Irms_A'], 'g-^', label='C3', linewidth=2, markersize=6)
    ax1.plot(freq, df_nom['L1_Irms_A'], 'm-d', label='L1', linewidth=2, markersize=6)
    ax1.plot(freq, df_nom['L2_Irms_A'], 'c-v', label='L2', linewidth=2, markersize=6)
    
    # Add rating lines
    ax1.axhline(1.0, color='orange', linestyle='--', linewidth=1.5, 
                label='Single Cap Rating (1A)')
    ax1.axhline(0.5, color='red', linestyle='--', linewidth=1.5, 
                label='Paralleled Cap Rating (0.5A each)')
    
    ax1.set_xlabel('Frequency (MHz)', fontsize=12)
    ax1.set_ylabel('RMS Current (A)', fontsize=12)
    ax1.set_title('Component Current Stress (Nominal Corner)', 
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10, loc='best')
    ax1.set_ylim(0, 1.2)
    
    # Voltage stress
    ax2.plot(freq, df_nom['C1_Vpk_V'], 'b-o', label='C1', linewidth=2, markersize=6)
    ax2.plot(freq, df_nom['C2_Vpk_V'], 'r-s', label='C2', linewidth=2, markersize=6)
    ax2.plot(freq, df_nom['C3_Vpk_V'], 'g-^', label='C3', linewidth=2, markersize=6)
    
    # Add rating line
    ax2.axhline(1000, color='red', linestyle='--', linewidth=1.5, 
                label='Voltage Rating (1000V)')
    
    ax2.set_xlabel('Frequency (MHz)', fontsize=12)
    ax2.set_ylabel('Peak Voltage (V)', fontsize=12)
    ax2.set_title('Capacitor Voltage Stress (Nominal Corner)', 
                  fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10, loc='best')
    ax2.set_ylim(0, 250)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/component_stress_nominal.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/component_stress_nominal.png")
    plt.close()

def plot_filter_dissipation(df, output_dir='plots'):
    """Plot filter power dissipation"""
    Path(output_dir).mkdir(exist_ok=True)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    corners = df['Corner'].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(corners)))
    
    for corner, color in zip(corners, colors):
        corner_data = df[df['Corner'] == corner].sort_values('Frequency_MHz')
        
        marker = 'o' if corner == 'Nominal' else 's'
        linewidth = 2.5 if corner == 'Nominal' else 1.5
        
        ax.plot(corner_data['Frequency_MHz'], 
                corner_data['Filter_Diss_mW'],
                marker=marker,
                linewidth=linewidth,
                label=corner,
                color=color,
                markersize=6)
    
    ax.set_xlabel('Frequency (MHz)', fontsize=12)
    ax.set_ylabel('Power Dissipation (mW)', fontsize=12)
    ax.set_title('LPF Power Dissipation vs Frequency (All Tolerance Corners)', 
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/filter_dissipation_all_corners.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/filter_dissipation_all_corners.png")
    plt.close()

def plot_worst_case_comparison(df, output_dir='plots'):
    """Compare best vs worst case harmonic performance"""
    Path(output_dir).mkdir(exist_ok=True)
    
    # Find worst H3 for each frequency
    h3_data = df.groupby('Frequency_MHz').agg({
        'H3_dBc': ['min', 'max'],
        'Band': 'first',
        'Corner': lambda x: x.iloc[0]  # Just need one
    }).reset_index()
    
    h3_data.columns = ['Frequency_MHz', 'H3_min', 'H3_max', 'Band', 'Corner']
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot range
    freq = h3_data['Frequency_MHz']
    ax.fill_between(freq, h3_data['H3_min'], h3_data['H3_max'], 
                     alpha=0.3, color='blue', label='Tolerance Range')
    
    # Plot nominal (from middle of range as approximation)
    ax.plot(freq, (h3_data['H3_min'] + h3_data['H3_max']) / 2, 
            'b-o', linewidth=2.5, markersize=6, label='Nominal (approx)')
    
    # FCC limit
    ax.axhline(-43, color='red', linestyle='--', linewidth=2, 
               label='FCC Requirement (-43 dBc)')
    
    # Design target
    ax.axhline(-60, color='green', linestyle='--', linewidth=2, 
               label='Design Target (-60 dBc)')
    
    ax.set_xlabel('Frequency (MHz)', fontsize=12)
    ax.set_ylabel('3rd Harmonic Level (dBc)', fontsize=12)
    ax.set_title('3rd Harmonic Suppression: Best vs Worst Case', 
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11, loc='best')
    ax.set_ylim(-85, -35)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/h3_worst_case_comparison.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/h3_worst_case_comparison.png")
    plt.close()

def generate_summary_table(df, output_dir='plots'):
    """Generate summary statistics table"""
    summary = []
    
    bands = df['Band'].unique()
    
    for band in bands:
        band_data = df[df['Band'] == band]
        nom_data = band_data[band_data['Corner'] == 'Nominal'].iloc[0]
        
        # Get worst case H3 across all corners for this band
        worst_h3 = band_data['H3_dBc'].max()  # max is least negative = worst
        
        summary.append({
            'Band': band,
            'Freq (MHz)': f"{nom_data['Frequency_MHz']:.2f}",
            'Fund Power (W)': f"{nom_data['Fund_Power_W']:.2f}",
            'Insertion Loss (dB)': f"{nom_data['Insertion_Loss_dB']:.3f}",
            'H2 (dBc)': f"{nom_data['H2_dBc']:.1f}",
            'H3 (dBc) Nom': f"{nom_data['H3_dBc']:.1f}",
            'H3 (dBc) Worst': f"{worst_h3:.1f}",
            'H5 (dBc)': f"{nom_data['H5_dBc']:.1f}",
            'Dissipation (mW)': f"{nom_data['Filter_Diss_mW']:.1f}",
            'C2 Current (A)': f"{nom_data['C2_Irms_A']:.3f}"
        })
    
    summary_df = pd.DataFrame(summary)
    
    # Save to CSV
    summary_df.to_csv(f'{output_dir}/summary_table.csv', index=False)
    print(f"\nSaved: {output_dir}/summary_table.csv")
    
    # Print to console
    print("\n" + "="*80)
    print("PERFORMANCE SUMMARY TABLE (Bottom of Each Band)")
    print("="*80)
    print(summary_df.to_string(index=False))
    print("="*80 + "\n")

def main():
    print("\n" + "="*80)
    print("NexRig Harmonic Analysis Plotter")
    print("="*80 + "\n")
    
    # Load results
    df = load_results()
    if df is None:
        return
    
    print(f"Loaded {len(df)} analysis results\n")
    
    # Generate all plots
    print("Generating plots...")
    plot_harmonic_spectrum(df)
    plot_insertion_loss(df)
    plot_component_stress(df)
    plot_filter_dissipation(df)
    plot_worst_case_comparison(df)
    generate_summary_table(df)
    
    print("\n" + "="*80)
    print("All plots generated successfully!")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
