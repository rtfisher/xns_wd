#!/usr/bin/env python3

"""
Script: parse_logfiles_and_plot.py

This script parses XNS output log files to extract angular velocity (OMG) and 
angular momentum values, converts them from geometric units (G=c=Msun=1) to 
CGS units, and creates a clean plot of angular velocity vs angular momentum.

The script looks for files matching "LogFile_*.dat" in the current directory.
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import namedtuple

# Physical constants for unit conversion
G_CGS = 6.67430e-8      # Gravitational constant [cm^3 g^-1 s^-2]
C_CGS = 2.99792458e10   # Speed of light [cm/s]
MSUN_CGS = 1.98847e33   # Solar mass [g]

# Conversion factors from geometric units to CGS
# In geometric units: G = c = Msun = 1
# Length unit: GM_sun/c^2 = 1.477 km
# Time unit: GM_sun/c^3 = 4.926 μs
# Angular velocity: [1/time] -> [1/s]
# Angular momentum: [M*L^2/T] -> [g*cm^2/s]

TIME_GEOM_TO_CGS = (G_CGS * MSUN_CGS) / (C_CGS**3)  # Convert geometric time to seconds
LENGTH_GEOM_TO_CGS = (G_CGS * MSUN_CGS) / (C_CGS**2)  # Convert geometric length to cm

# Angular velocity: 1/time_geom -> 1/s
OMEGA_GEOM_TO_CGS = 1.0 / TIME_GEOM_TO_CGS

# Angular momentum: M_sun * L_geom^2 / T_geom -> g * cm^2 / s
L_GEOM_TO_CGS = MSUN_CGS * (LENGTH_GEOM_TO_CGS**2) / TIME_GEOM_TO_CGS

print(f"Conversion factors:")
print(f"Time (geom -> s): {TIME_GEOM_TO_CGS:.6e}")
print(f"Length (geom -> cm): {LENGTH_GEOM_TO_CGS:.6e}")
print(f"Angular velocity (geom -> 1/s): {OMEGA_GEOM_TO_CGS:.6e}")
print(f"Angular momentum (geom -> g⋅cm²/s): {L_GEOM_TO_CGS:.6e}")
print()

# Data structure to hold parsed values
ModelData = namedtuple('ModelData', ['filename', 'omega_geom', 'L_geom', 'omega_cgs', 'L_cgs'])

def parse_logfile(filepath):
    """
    Parse a single log file to extract OMG and ANGUL. MOMENT. values.
    
    Args:
        filepath: Path to the log file
        
    Returns:
        tuple: (omega_geom, L_geom) or (None, None) if parsing fails
    """
    omega_geom = None
    L_geom = None
    
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            
        # Search for OMG value
        omega_match = re.search(r'OMG\s*=\s*([+-]?\d+\.?\d*[eE]?[+-]?\d*)', content)
        if omega_match:
            omega_geom = float(omega_match.group(1))
            
        # Search for angular momentum value
        L_match = re.search(r'ANGUL\.\s*MOMENT\.\s*\(E\)\s*=\s*([+-]?\d+\.?\d*[eE]?[+-]?\d*)', content)
        if L_match:
            L_geom = float(L_match.group(1))
            
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
        return None, None
        
    return omega_geom, L_geom

def main():
    """Main function to parse all log files and create the plot."""
    
    # Find all LogFile_*.dat files in current directory
    current_dir = Path('.')
    logfiles = list(current_dir.glob('LogFile_*.dat'))
    
    if not logfiles:
        print("No LogFile_*.dat files found in current directory!")
        return
        
    print(f"Found {len(logfiles)} log files to process:")
    for f in logfiles:
        print(f"  {f}")
    print()
    
    # Parse all files
    models = []
    failed_files = []
    
    for filepath in logfiles:
        omega_geom, L_geom = parse_logfile(filepath)
        
        if omega_geom is not None and L_geom is not None:
            # Convert to CGS units
            omega_cgs = omega_geom * OMEGA_GEOM_TO_CGS
            L_cgs = L_geom * L_GEOM_TO_CGS
            
            models.append(ModelData(
                filename=filepath.name,
                omega_geom=omega_geom,
                L_geom=L_geom,
                omega_cgs=omega_cgs,
                L_cgs=L_cgs
            ))
            
            print(f"{filepath.name}:")
            print(f"  OMG (geom) = {omega_geom:.6e}")
            print(f"  L (geom)   = {L_geom:.6e}")
            print(f"  OMG (CGS)  = {omega_cgs:.6e} s⁻¹")
            print(f"  L (CGS)    = {L_cgs:.6e} g⋅cm²⋅s⁻¹")
            print()
        else:
            failed_files.append(filepath.name)
            print(f"Failed to parse {filepath.name}")
    
    if failed_files:
        print(f"Failed to parse {len(failed_files)} files: {failed_files}")
        
    if not models:
        print("No valid data found!")
        return
    
    print(f"Successfully parsed {len(models)} models")
    
    # Sort models by angular momentum for clean plotting
    models.sort(key=lambda x: x.L_cgs)
    
    # Extract data for plotting
    omega_cgs_values = [m.omega_cgs for m in models]
    L_cgs_values = [m.L_cgs for m in models]
    
    # Create the plot
    plt.figure(figsize=(10, 7))
    
    # Main plot
    plt.plot(L_cgs_values, omega_cgs_values, 'bo-', linewidth=2, markersize=6, 
             markerfacecolor='lightblue', markeredgecolor='blue', markeredgewidth=1)
    
    # Formatting
    plt.xlabel('Angular Momentum [g⋅cm²⋅s⁻¹]', fontsize=12, fontweight='bold')
    plt.ylabel('Angular Velocity [s⁻¹]', fontsize=12, fontweight='bold')
    plt.title('Angular Velocity vs Angular Momentum\nRotating White Dwarf Models', 
              fontsize=14, fontweight='bold', pad=20)
    
    # Use scientific notation for axes
    plt.ticklabel_format(style='scientific', axis='both', scilimits=(0,0))
    
    # Grid
    plt.grid(True, alpha=0.3, linestyle='--')
    
    # Add annotation with number of models
    plt.text(0.02, 0.98, f'N = {len(models)} models', 
             transform=plt.gca().transAxes, fontsize=10,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Add conversion info as text box
    conversion_text = (f'Converted from geometric units (G=c=M☉=1)\n'
                      f'ω: ×{OMEGA_GEOM_TO_CGS:.2e} → s⁻¹\n'
                      f'L: ×{L_GEOM_TO_CGS:.2e} → g⋅cm²⋅s⁻¹')
    plt.text(0.98, 0.02, conversion_text, 
             transform=plt.gca().transAxes, fontsize=9,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.7))
    
    # Tight layout
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('angular_velocity_vs_momentum.png', dpi=300, bbox_inches='tight')
    plt.savefig('angular_velocity_vs_momentum.pdf', bbox_inches='tight')
    
    print(f"\nPlot saved as 'angular_velocity_vs_momentum.png' and '.pdf'")
    
    # Show the plot
    plt.show()
    
    # Print summary statistics
    print("\n" + "="*50)
    print("SUMMARY STATISTICS")
    print("="*50)
    print(f"Number of models: {len(models)}")
    print(f"Angular velocity range: {min(omega_cgs_values):.3e} to {max(omega_cgs_values):.3e} s⁻¹")
    print(f"Angular momentum range: {min(L_cgs_values):.3e} to {max(L_cgs_values):.3e} g⋅cm²⋅s⁻¹")
    
    # Check if omega=0 case exists
    zero_omega_models = [m for m in models if abs(m.omega_cgs) < 1e-10]
    if zero_omega_models:
        print(f"Non-rotating model (ω≈0): L = {zero_omega_models[0].L_cgs:.3e} g⋅cm²⋅s⁻¹")

if __name__ == "__main__":
    main()
