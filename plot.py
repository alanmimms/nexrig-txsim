#!/usr/bin/env python3
"""
NexRig Harmonic Analysis Plotter
Generates publication-quality plots from harmonics.csv
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def load_results(filename='harmonics.csv'):
    """Load analysis results from CSV file"""
    if not Path(filename).exists():
        print(f"Error: {filename} not found. Run ./harmonics first.")
        return None

    df = pd.read_csv(filename)
    return df

def plot_harmonic_spectrum(df, output_dir='plots'):
    """Plot harmonic content for each band (nominal corner only)"""
    Path(output_dir).mkdir(exist_ok=True)

    # Filter for nominal corner only
    df_nom = df[df['Corner'] == 'Nominal'].copy()

    bands = df_nom['Band'].unique()

    # Use 2 rows x 5 columns to accommodate 10 bands
    fig, axes = plt.subplots(2, 5, figsize=(20, 10))
    axes = axes.flatten()

    for idx, band in enumerate(bands):
        ax = axes[idx]
        band_data = df_nom[df_nom['Band'] == band]

        # Get bottom and top of band
        freq_low = band_data.iloc[0]
        freq_high = band_data.iloc[1] if len(band_data) > 1 else freq_low

        # Plot harmonics for bottom of band
        harmonics = range(2, 11)
        h_cols = [f'H{n}dBc' for n in harmonics]

        values_low = [freq_low[col] for col in h_cols]
        values_high = [freq_high[col] for col in h_cols]

        ax.plot(harmonics, values_low, 'b-o', label=f'{freq_low["fMHz"]:.2f} MHz', linewidth=2)
        ax.plot(harmonics, values_high, 'r-s', label=f'{freq_high["fMHz"]:.2f} MHz', linewidth=2)

        # Add FCC requirement line
        ax.axhline(-43, color='g', linestyle='--', linewidth=1.5, label='FCC Limit')

        ax.set_xlabel('Harmonic Number', fontsize=10)
        ax.set_ylabel('Level (dBc)', fontsize=10)
        ax.set_title(f'{band} Band', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        ax.set_ylim(-160, -40)
        ax.set_xticks(harmonics)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/harmonic-spectrum-all-bands.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/harmonic-spectrum-all-bands.png")
    plt.close()

def plot_insertion_loss(df, output_dir='plots'):
    """Plot insertion loss vs frequency for all corners"""
    Path(output_dir).mkdir(exist_ok=True)

    fig, ax = plt.subplots(figsize=(14, 8))

    corners = df['Corner'].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(corners)))

    for corner, color in zip(corners, colors):
        corner_data = df[df['Corner'] == corner].sort_values('fMHz')

        # Plot with markers
        marker = 'o' if corner == 'Nominal' else 's'
        linewidth = 2.5 if corner == 'Nominal' else 1.5

        ax.plot(corner_data['fMHz'],
                corner_data['InsLossdB'],
                marker=marker,
                linewidth=linewidth,
                label=corner,
                color=color,
                markersize=6)

    # Add vertical lines and labels for band boundaries
    # Extract band info from nominal corner
    df_nom = df[df['Corner'] == 'Nominal'].sort_values('fMHz')
    bands = df_nom['Band'].unique()

    band_ranges = {}
    for band in bands:
        band_data = df_nom[df_nom['Band'] == band]
        fmin = band_data['fMHz'].min()
        fmax = band_data['fMHz'].max()
        band_ranges[band] = (fmin, fmax)

    # Draw vertical lines at band boundaries and label regions
    ylim = ax.get_ylim()
    for i, (band, (fmin, fmax)) in enumerate(band_ranges.items()):
        # Vertical lines at boundaries
        ax.axvline(fmin, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        ax.axvline(fmax, color='gray', linestyle='--', alpha=0.5, linewidth=1)

        # Band label in the middle
        fmid = (fmin + fmax) / 2
        ax.text(fmid, ylim[0] + (ylim[1] - ylim[0]) * 0.95, band,
                ha='center', va='top', fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))

    ax.set_xlabel('Frequency (MHz)', fontsize=12)
    ax.set_ylabel('Insertion Loss (dB)', fontsize=12)
    ax.set_title('LPF Insertion Loss vs Frequency (All Tolerance Corners)',
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best')

    # Format y-axis to show losses as positive numbers
    ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(f'{output_dir}/insertion-loss-all-corners.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/insertion-loss-all-corners.png")
    plt.close()

def plot_component_stress(df, output_dir='plots'):
    """Plot component current stress for nominal corner"""
    Path(output_dir).mkdir(exist_ok=True)

    df_nom = df[df['Corner'] == 'Nominal'].copy().sort_values('fMHz')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

    # Current stress
    freq = df_nom['fMHz']

    ax1.plot(freq, df_nom['C1IrmsA'], 'b-o', label='C1', linewidth=2, markersize=6)
    ax1.plot(freq, df_nom['C2IrmsA'], 'r-s', label='C2', linewidth=2, markersize=6)
    ax1.plot(freq, df_nom['C3IrmsA'], 'g-^', label='C3', linewidth=2, markersize=6)
    ax1.plot(freq, df_nom['C4IrmsA'], 'y-x', label='C4', linewidth=2, markersize=6)
    ax1.plot(freq, df_nom['L1IrmsA'], 'm-d', label='L1', linewidth=2, markersize=6)
    ax1.plot(freq, df_nom['L2IrmsA'], 'c-v', label='L2', linewidth=2, markersize=6)
    ax1.plot(freq, df_nom['L3IrmsA'], 'k-p', label='L3', linewidth=2, markersize=6)

    # Add rating lines
    ax1.axhline(2.0, color='red', linestyle='--', linewidth=1.5,
                label='Inductor Saturation (2A)')

    # Add band boundaries and labels
    bands = df_nom['Band'].unique()
    band_ranges = {}
    for band in bands:
        band_data = df_nom[df_nom['Band'] == band]
        fmin = band_data['fMHz'].min()
        fmax = band_data['fMHz'].max()
        band_ranges[band] = (fmin, fmax)

    ylim1 = ax1.get_ylim()
    for band, (fmin, fmax) in band_ranges.items():
        ax1.axvline(fmin, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        ax1.axvline(fmax, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        fmid = (fmin + fmax) / 2
        ax1.text(fmid, ylim1[1] - (ylim1[1] - ylim1[0]) * 0.05, band,
                ha='center', va='top', fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))

    ax1.set_xlabel('Frequency (MHz)', fontsize=12)
    ax1.set_ylabel('RMS Current (A)', fontsize=12)
    ax1.set_title('Component Current Stress (Nominal Corner)',
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10, loc='best')
    ax1.set_ylim(0, 2.5)

    # Voltage stress
    ax2.plot(freq, df_nom['C1VpkV'], 'b-o', label='C1', linewidth=2, markersize=6)
    ax2.plot(freq, df_nom['C2VpkV'], 'r-s', label='C2', linewidth=2, markersize=6)
    ax2.plot(freq, df_nom['C3VpkV'], 'g-^', label='C3', linewidth=2, markersize=6)
    ax2.plot(freq, df_nom['C4VpkV'], 'y-x', label='C4', linewidth=2, markersize=6)

    # Add rating line
    ax2.axhline(1000, color='red', linestyle='--', linewidth=1.5,
                label='Voltage Rating (1000V)')

    # Add band boundaries and labels
    ylim2 = ax2.get_ylim()
    for band, (fmin, fmax) in band_ranges.items():
        ax2.axvline(fmin, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        ax2.axvline(fmax, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        fmid = (fmin + fmax) / 2
        ax2.text(fmid, ylim2[1] - (ylim2[1] - ylim2[0]) * 0.05, band,
                ha='center', va='top', fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))

    ax2.set_xlabel('Frequency (MHz)', fontsize=12)
    ax2.set_ylabel('Peak Voltage (V)', fontsize=12)
    ax2.set_title('Capacitor Voltage Stress (Nominal Corner)',
                  fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10, loc='best')
    ax2.set_ylim(0, 400)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/component-stress-nominal.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/component-stress-nominal.png")
    plt.close()

def plot_filter_dissipation(df, output_dir='plots'):
    """Plot filter power dissipation"""
    Path(output_dir).mkdir(exist_ok=True)

    fig, ax = plt.subplots(figsize=(14, 8))

    corners = df['Corner'].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(corners)))

    for corner, color in zip(corners, colors):
        corner_data = df[df['Corner'] == corner].sort_values('fMHz')

        marker = 'o' if corner == 'Nominal' else 's'
        linewidth = 2.5 if corner == 'Nominal' else 1.5

        ax.plot(corner_data['fMHz'],
                corner_data['FilterDissmW'],
                marker=marker,
                linewidth=linewidth,
                label=corner,
                color=color,
                markersize=6)

    # Add vertical lines and labels for band boundaries
    # Extract band info from nominal corner
    df_nom = df[df['Corner'] == 'Nominal'].sort_values('fMHz')
    bands = df_nom['Band'].unique()

    band_ranges = {}
    for band in bands:
        band_data = df_nom[df_nom['Band'] == band]
        fmin = band_data['fMHz'].min()
        fmax = band_data['fMHz'].max()
        band_ranges[band] = (fmin, fmax)

    # Draw vertical lines at band boundaries and label regions
    ylim = ax.get_ylim()
    for i, (band, (fmin, fmax)) in enumerate(band_ranges.items()):
        # Vertical lines at boundaries
        ax.axvline(fmin, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        ax.axvline(fmax, color='gray', linestyle='--', alpha=0.5, linewidth=1)

        # Band label in the middle
        fmid = (fmin + fmax) / 2
        ax.text(fmid, ylim[1] - (ylim[1] - ylim[0]) * 0.05, band,
                ha='center', va='top', fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))

    ax.set_xlabel('Frequency (MHz)', fontsize=12)
    ax.set_ylabel('Power Dissipation (mW)', fontsize=12)
    ax.set_title('LPF Power Dissipation vs Frequency (All Tolerance Corners)',
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best')

    plt.tight_layout()
    plt.savefig(f'{output_dir}/filter-dissipation-all-corners.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/filter-dissipation-all-corners.png")
    plt.close()

def plot_worst_case_comparison(df, output_dir='plots'):
    """Compare best vs worst case harmonic performance"""
    Path(output_dir).mkdir(exist_ok=True)

    # Find worst H3 for each frequency
    h3_data = df.groupby('fMHz').agg({
        'H3dBc': ['min', 'max'],
        'Band': 'first',
        'Corner': lambda x: x.iloc[0]  # Just need one
    }).reset_index()

    h3_data.columns = ['fMHz', 'H3min', 'H3max', 'Band', 'Corner']

    fig, ax = plt.subplots(figsize=(14, 8))

    # Plot range
    freq = h3_data['fMHz']
    ax.fill_between(freq, h3_data['H3min'], h3_data['H3max'],
                     alpha=0.3, color='blue', label='Tolerance Range')

    # Plot nominal (from middle of range as approximation)
    ax.plot(freq, (h3_data['H3min'] + h3_data['H3max']) / 2,
            'b-o', linewidth=2.5, markersize=6, label='Nominal (approx)')

    # FCC limit
    ax.axhline(-43, color='red', linestyle='--', linewidth=2,
               label='FCC Requirement (-43 dBc)')

    # Design target
    ax.axhline(-60, color='green', linestyle='--', linewidth=2,
               label='Design Target (-60 dBc)')

    # Add band boundaries and labels
    df_nom = df[df['Corner'] == 'Nominal'].sort_values('fMHz')
    bands = df_nom['Band'].unique()
    band_ranges = {}
    for band in bands:
        band_data = df_nom[df_nom['Band'] == band]
        fmin = band_data['fMHz'].min()
        fmax = band_data['fMHz'].max()
        band_ranges[band] = (fmin, fmax)

    ylim = ax.get_ylim()
    for band, (fmin, fmax) in band_ranges.items():
        ax.axvline(fmin, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        ax.axvline(fmax, color='gray', linestyle='--', alpha=0.5, linewidth=1)
        fmid = (fmin + fmax) / 2
        ax.text(fmid, ylim[0] + (ylim[1] - ylim[0]) * 0.95, band,
                ha='center', va='top', fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))

    ax.set_xlabel('Frequency (MHz)', fontsize=12)
    ax.set_ylabel('3rd Harmonic Level (dBc)', fontsize=12)
    ax.set_title('3rd Harmonic Suppression: Best vs Worst Case',
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11, loc='best')
    ax.set_ylim(-90, -35)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/h3-worst-case-comparison.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/h3-worst-case-comparison.png")
    plt.close()

def generate_summary_table(df, output_dir='plots'):
    """Generate summary statistics table"""
    summary = []

    bands = df['Band'].unique()

    for band in bands:
        band_data = df[df['Band'] == band]
        nom_data = band_data[band_data['Corner'] == 'Nominal'].iloc[0]

        # Get worst case H3 across all corners for this band
        worst_h3 = band_data['H3dBc'].max()  # max is least negative = worst

        summary.append({
            'Band': band,
            'Freq (MHz)': f"{nom_data['fMHz']:.2f}",
            'Fund Power (W)': f"{nom_data['FundW']:.2f}",
            'Insertion Loss (dB)': f"{nom_data['InsLossdB']:.3f}",
            'H2 (dBc)': f"{nom_data['H2dBc']:.1f}",
            'H3 (dBc) Nom': f"{nom_data['H3dBc']:.1f}",
            'H3 (dBc) Worst': f"{worst_h3:.1f}",
            'H5 (dBc)': f"{nom_data['H5dBc']:.1f}",
            'H7 (dBc)': f"{nom_data['H7dBc']:.1f}",
            'Dissipation (mW)': f"{nom_data['FilterDissmW']:.1f}",
            'C2 Current (A)': f"{nom_data['C2IrmsA']:.3f}"
        })

    summary_df = pd.DataFrame(summary)

    # Save to CSV
    summary_df.to_csv(f'{output_dir}/summary-table.csv', index=False)
    print(f"\nSaved: {output_dir}/summary-table.csv")

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
