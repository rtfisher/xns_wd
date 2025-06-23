"""
This python routine reads a CGS format EoS table (wd_eos_helmholtz.dat) 
and performs a resample of the EoS with uniform spacing in Log10-Log10 scale
both as a function of density, pressure and enthalpy.

Input table format: density[g/cm³], pressure[dyne/cm²], internal_energy_density[erg/cm³]
Lines starting with # are treated as comments and ignored.

Output is in geometrized units suitable for neutron star structure calculations.

OUTPUT FORMAT (XNS 4.0 compatible):
The output file contains three sections for different parameterizations:
1. Section 1 (1000 points): Uniform log10(ρ) spacing with P(ρ), e(ρ), H(ρ) and power-law indices
2. Section 2 (1000 points): Uniform log10(P) spacing with ρ(P) and power-law index  
3. Section 3 (1000 points): Uniform log10(H) spacing with ρ(H) and power-law index
All quantities are in geometrized units (G=c=1) with dimensions expressible as powers of length.
The format is compatible with XNS 4.0 neutron star structure code.

PSEUDO-ENTHALPY CALCULATION:
The pseudo-enthalpy H(P) is computed using the exact thermodynamic result:
    H(P) = c² ln(h/c²)
where h is the physical specific enthalpy:
    h = e/ρ + P/ρ + c²
    
This replaces the complex numerical integration approach with a direct calculation
based on fundamental thermodynamic relations.

ASSUMPTIONS:
This exact result applies under the following conditions:
1. ISENTROPIC FLOW: The fluid is assumed to have constant specific entropy (ds = 0)
2. PERFECT FLUID: No viscosity, heat conduction, or bulk viscosity effects
3. LOCAL THERMODYNAMIC EQUILIBRIUM: Standard thermodynamic relations apply locally
4. BAROTROPIC EoS: Pressure depends only on density P = P(ρ)
5. BOUNDARY CONDITION: h(0) = c² (physical specific enthalpy approaches rest mass energy at zero pressure)

These assumptions are typically valid for:
- Cold stellar matter (white dwarfs, neutron stars)
- Isentropic stellar cores
- Single-phase fluids without composition gradients

The assumptions may break down for:
- Multi-temperature systems
- Multi-component fluids with varying composition
- Systems with significant entropy gradients
- Phase transitions or mixed phases
- Non-equilibrium processes
"""
import os
import numpy as np
import numpy.ma as ma
import scipy as sp
from scipy.interpolate import interp1d
from scipy.integrate import quad

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.ion()

# Input file name
eosname = 'wd_eos_helmholtz'
input_file = eosname + '.dat'

# Physical Constants
Msun = 1.99e33
G = 6.6727e-8
c = 3e10

# Geometrized unit definitions
rhounit = Msun/(G*Msun/c**2)**3  # Density unit
punit = rhounit*c**2             # Pressure/energy unit

# Read the CGS data
density_cgs = []
pressure_cgs = []
energy_internal_cgs = []

print(f"Reading EoS data from {input_file}...")

with open(input_file, 'r') as f:
    line_count = 0
    data_count = 0
    for line in f:
        line_count += 1
        line = line.strip()
        
        # Skip empty lines and comments
        if not line or line.startswith('#'):
            continue
            
        try:
            data = line.split()
            if len(data) < 3:
                print(f"Warning: Line {line_count} has insufficient data, skipping: {line}")
                continue
                
            density_cgs.append(float(data[0]))      # g/cm³
            pressure_cgs.append(float(data[1]))     # dyne/cm²
            energy_internal_cgs.append(float(data[2]))  # erg/cm³
            data_count += 1
            
        except ValueError:
            print(f"Warning: Could not parse line {line_count}, skipping: {line}")
            continue

print(f"Successfully read {data_count} data points from {line_count} lines")

# Convert to numpy arrays
density_cgs = np.array(density_cgs)
pressure_cgs = np.array(pressure_cgs)
energy_internal_cgs = np.array(energy_internal_cgs)

# Print data ranges for verification
print(f"\nData ranges:")
print(f"Density: {np.min(density_cgs):.3e} to {np.max(density_cgs):.3e} g/cm³")
print(f"Pressure: {np.min(pressure_cgs):.3e} to {np.max(pressure_cgs):.3e} dyne/cm²")
print(f"Internal Energy: {np.min(energy_internal_cgs):.3e} to {np.max(energy_internal_cgs):.3e} erg/cm³")

# Physical reasonableness checks
nuclear_density = 2.8e14  # g/cm³
if np.max(density_cgs) > nuclear_density:
    print(f"Note: Maximum density ({np.max(density_cgs):.2e} g/cm³) exceeds nuclear density")

# Convert to geometrized units
density = density_cgs / rhounit
pressure = pressure_cgs / punit
energy = energy_internal_cgs / punit

print(f"\nGeometrized unit conversion factors:")
print(f"rhounit = {rhounit:.3e} g/cm³")
print(f"punit = {punit:.3e} dyne/cm²")

# Validation check for internal energy
if np.any(energy < 0):
    negative_indices = np.where(energy < 0)[0]
    print(f"ERROR: Negative internal energy detected at {len(negative_indices)} points!")
    print(f"First occurrence at density = {density_cgs[negative_indices[0]]:.3e} g/cm³")
    print("Check your input data - internal energy should always be positive")
else:
    print("✓ All internal energy values are positive")

# Get the range of the EoS    
minlogrho = np.min(np.log10(density))    
maxlogrho = np.max(np.log10(density))

minlogprs = np.min(np.log10(pressure))    
maxlogprs = np.max(np.log10(pressure))

print(f"\nLog10 ranges in geometrized units:")
print(f"Log10(density): {minlogrho:.3f} to {maxlogrho:.3f}")
print(f"Log10(pressure): {minlogprs:.3f} to {maxlogprs:.3f}")

# Define the resampling points for density in Log10 scale
rhopoints = np.linspace(minlogrho, maxlogrho, 1000)

# Define the interpolating function for Log10(P) and Log10(e)
rho2prs = interp1d(np.log10(density), np.log10(pressure), kind='linear')
rho2ein = interp1d(np.log10(density), np.log10(energy), kind='linear')

# Interpolate over the density points
prspoints = rho2prs(rhopoints)
einpoints = rho2ein(rhopoints)

# Compute the local powerlaw index of the EoS (p ~ rho^g; e~ rho^l) in each density interval
prsindex = np.zeros_like(rhopoints)
einindex = np.zeros_like(rhopoints)

prsindex[0:999] = (prspoints[1:1000] - prspoints[0:999]) / (rhopoints[1:1000] - rhopoints[0:999])
prsindex[999] = prsindex[998]

einindex[0:999] = (einpoints[1:1000] - einpoints[0:999]) / (rhopoints[1:1000] - rhopoints[0:999])
einindex[999] = einindex[998]

# Define dp/drho and dh/drho
dprsdrho = 10**prspoints[0:1000] * prsindex[0:1000] / (10**rhopoints[0:1000])
dhdrho = dprsdrho / (10**rhopoints[0:1000] + 10**prspoints[0:1000] + 10**einpoints[0:1000])

# Compute pseudo-enthalpy using exact thermodynamic result
print("\nComputing pseudo-enthalpy using exact thermodynamic relation...")
print("Using H(P) = c² ln(h/c²) where h = e/ρ + P/ρ + c²")
print("This assumes isentropic, perfect fluid conditions")

# Compute physical specific enthalpy h = e/ρ + P/ρ + c²
# Note: e is internal energy density, so e/ρ gives specific internal energy
specific_internal_energy = energy_internal_cgs / density_cgs  # erg/g
specific_pressure_work = pressure_cgs / density_cgs           # erg/g (equivalent to cm²/s²)
specific_rest_mass_energy = c**2                              # cm²/s²

# Physical specific enthalpy (including rest mass energy)
h_physical = specific_internal_energy + specific_pressure_work + specific_rest_mass_energy

# Validate physical specific enthalpy
if np.any(h_physical <= 0):
    print("ERROR: Non-positive physical specific enthalpy detected!")
    negative_indices = np.where(h_physical <= 0)[0]
    print(f"First occurrence at density = {density_cgs[negative_indices[0]]:.3e} g/cm³")
    print("This violates the fundamental requirement h > 0")
    raise ValueError("Invalid physical specific enthalpy")

if np.any(h_physical < c**2):
    print("WARNING: Physical specific enthalpy less than c² detected!")
    print("This may indicate issues with the boundary condition h(0) = c²")
    low_h_indices = np.where(h_physical < c**2)[0]
    print(f"Minimum h/c² ratio: {np.min(h_physical/c**2):.6f}")

# Compute pseudo-enthalpy (enthalpy potential) H(P) = c² ln(h/c²)
H_pseudo = c**2 * np.log(h_physical / c**2)

print(f"Physical specific enthalpy range: {np.min(h_physical):.3e} to {np.max(h_physical):.3e} erg/g")
print(f"h/c² ratio range: {np.min(h_physical/c**2):.6f} to {np.max(h_physical/c**2):.6f}")
print(f"Pseudo-enthalpy H(P) range: {np.min(H_pseudo):.3e} to {np.max(H_pseudo):.3e} erg/g")

# Convert pseudo-enthalpy to geometrized units and take log10
# Note: H_pseudo has units of energy per mass (erg/g), same as c²
H_pseudo_geometrized = H_pseudo / c**2  # Now dimensionless

# Create interpolation for pseudo-enthalpy on the resampled density grid
H_pseudo_interp = interp1d(np.log10(density), np.log10(H_pseudo_geometrized), kind='linear')
H_pseudo_points = H_pseudo_interp(rhopoints)

# For consistency with original format, we compute log10 of the pseudo-enthalpy
entpoints = H_pseudo_points  # Already in log10 form

# Compute the powerlaw index of the EoS (H ~ rho^k)
entindex = np.zeros_like(rhopoints)
entindex[0:999] = (entpoints[1:1000] - entpoints[0:999]) / (rhopoints[1:1000] - rhopoints[0:999])
entindex[999] = entindex[998]

# Validation checks
valid_ent = np.isfinite(entpoints)
if not np.all(valid_ent):
    print(f"Warning: {np.sum(~valid_ent)} pseudo-enthalpy values are not finite")
else:
    print("✓ All pseudo-enthalpy values are finite and valid")

print("✓ Pseudo-enthalpy computation completed using exact thermodynamic relation")

# Create diagnostic plots
# Create diagnostic plots and save them
print("\nGenerating diagnostic plots...")

# Create subdirectory for diagnostic plots
import os
diag_dir = "_resampled_eos_diagnostics"
if not os.path.exists(diag_dir):
    os.makedirs(diag_dir)
    print(f"Created directory: {diag_dir}")
else:
    print(f"Using existing directory: {diag_dir}")

plt.figure('Resampled EoS', figsize=(12, 8))

# Main EoS plot
plt.subplot(2, 2, 1)
plt.plot(rhopoints, prspoints, label='Log10(P)', linewidth=2)
plt.plot(rhopoints, entpoints, label='Log10(H)', linewidth=2)
plt.plot(rhopoints, einpoints, label='Log10(e)', linewidth=2)
plt.xlabel('Log10[density] (geometrized units)')
plt.ylabel('Log10[pressure], Log10[energy], Log10[H]')
plt.legend()
plt.grid(True, alpha=0.3)
plt.title('Resampled EoS in Log-Log Scale')

# Physical specific enthalpy validation plot
plt.subplot(2, 2, 2)
plt.semilogx(density_cgs, h_physical/c**2, 'b-', linewidth=2)
plt.xlabel('Density [g/cm³]')
plt.ylabel('h/c² (dimensionless)')
plt.title('Physical Specific Enthalpy Ratio')
plt.grid(True, alpha=0.3)
plt.axhline(y=1, color='r', linestyle='--', alpha=0.7, label='h = c²')
plt.legend()

# Pseudo-enthalpy plot
plt.subplot(2, 2, 3)
plt.semilogx(density_cgs, H_pseudo/c**2, 'g-', linewidth=2)
plt.xlabel('Density [g/cm³]')
plt.ylabel('H/c² (dimensionless)')
plt.title('Pseudo-enthalpy H(P)/c²')
plt.grid(True, alpha=0.3)

# Power-law indices
plt.subplot(2, 2, 4)
plt.plot(rhopoints, prsindex, label='γ (pressure)', linewidth=2)
plt.plot(rhopoints, einindex, label='λ (energy)', linewidth=2)
plt.plot(rhopoints, entindex, label='κ (pseudo-enthalpy)', linewidth=2)
plt.xlabel('Log10[density] (geometrized units)')
plt.ylabel('Power-law index')
plt.legend()
plt.grid(True, alpha=0.3)
plt.title('Local Power-law Indices')

plt.tight_layout()

# Save the comprehensive diagnostic plot
main_plot_path = os.path.join(diag_dir, f"{eosname}_diagnostics.png")
plt.savefig(main_plot_path, dpi=300, bbox_inches='tight')
print(f"✓ Saved main diagnostic plot: {main_plot_path}")

# Save high-quality PDF version
main_plot_pdf = os.path.join(diag_dir, f"{eosname}_diagnostics.pdf")
plt.savefig(main_plot_pdf, bbox_inches='tight')
print(f"✓ Saved PDF version: {main_plot_pdf}")

plt.show()

# Create and save individual detailed plots
print("Generating additional detailed plots...")

# Individual P-rho plot
plt.figure('Pressure vs Density', figsize=(10, 6))
plt.loglog(density_cgs, pressure_cgs, 'bo-', markersize=2, linewidth=1, label='Original data')
plt.loglog(10**rhopoints * rhounit, 10**prspoints * punit, 'r-', linewidth=2, label='Resampled')
plt.xlabel('Density [g/cm³]')
plt.ylabel('Pressure [dyne/cm²]')
plt.title('Pressure vs Density (CGS units)')
plt.legend()
plt.grid(True, alpha=0.3)
p_rho_path = os.path.join(diag_dir, f"{eosname}_pressure_density.png")
plt.savefig(p_rho_path, dpi=300, bbox_inches='tight')
print(f"✓ Saved P-ρ plot: {p_rho_path}")
plt.close()

# Individual enthalpy validation plot
plt.figure('Enthalpy Validation', figsize=(10, 6))
plt.subplot(1, 2, 1)
plt.semilogx(density_cgs, h_physical, 'b-', linewidth=2)
plt.xlabel('Density [g/cm³]')
plt.ylabel('Physical Specific Enthalpy [erg/g]')
plt.title('Physical Specific Enthalpy h')
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
# Plot the pseudo-enthalpy : H(P)/c² 
# Since H_pseudo = c² ln(h/c²), we have H_pseudo/c² = ln(h/c²)
plt.loglog(density_cgs, H_pseudo/c**2, 'g-', linewidth=2)
plt.xlabel('Density [g/cm³]')
plt.ylabel('H(P)/c² (dimensionless)')
plt.title('Pseudo-enthalpy')
plt.grid(True, alpha=0.3)
#plt.axhline(y=-1, color='r', linestyle='--', alpha=0.7, label='H(P) = 0')
#plt.axhline(y=0, color='orange', linestyle='--', alpha=0.7, label='H(P) = c²')
#plt.legend()

plt.tight_layout()
enthalpy_path = os.path.join(diag_dir, f"{eosname}_enthalpy_validation.png")
plt.savefig(enthalpy_path, dpi=300, bbox_inches='tight')
print(f"✓ Saved enthalpy validation plot: {enthalpy_path}")
plt.close()

# Power-law indices detailed plot
plt.figure('Power-law Indices', figsize=(10, 6))
plt.plot(10**rhopoints * rhounit, prsindex, label='γ = d(log P)/d(log ρ)', linewidth=2)
plt.plot(10**rhopoints * rhounit, einindex, label='λ = d(log e)/d(log ρ)', linewidth=2)
plt.plot(10**rhopoints * rhounit, entindex, label='κ = d(log H)/d(log ρ)', linewidth=2)
plt.xlabel('Density [g/cm³]')
plt.ylabel('Power-law Index')
plt.title('Local Power-law Indices vs Density')
plt.legend()
plt.grid(True, alpha=0.3)
plt.semilogx()
indices_path = os.path.join(diag_dir, f"{eosname}_powerlaw_indices.png")
plt.savefig(indices_path, dpi=300, bbox_inches='tight')
print(f"✓ Saved power-law indices plot: {indices_path}")
plt.close()

print(f"✓ All diagnostic plots saved to directory: {diag_dir}")


#---------------------------------------------------------------------------------
# Define the resampling points for pressure in Log10 scale 
prspoints2 = np.linspace(minlogprs, maxlogprs, 1000)

# Define the interpolating function for Log10(rho)
prs2rho = interp1d(np.log10(pressure), np.log10(density), kind='linear')

# Interpolate over the pressure points
rhopoints2 = prs2rho(prspoints2)

# Compute the powerlaw index of the EoS (rho ~ p^m) in each pressure interval
rhoindex2 = np.zeros_like(prspoints2)
rhoindex2[0:999] = (rhopoints2[1:1000] - rhopoints2[0:999]) / (prspoints2[1:1000] - prspoints2[0:999])
rhoindex2[999] = rhoindex2[998]

#---------------------------------------------------------------------------------
# Define the resampling points for H(P) in Log10 scale
print("\nGenerating pseudo-enthalpy-based resampling...")

# Use the computed pseudo-enthalpy range
ent_min = np.min(entpoints)
ent_max = np.max(entpoints)
entpoints3 = np.linspace(ent_min, ent_max, 1000)

# Define the interpolating function for Log10(rho) from Log10(H)
ent2rho = interp1d(entpoints, rhopoints, kind='linear', fill_value='extrapolate')

# Interpolate over the pseudo-enthalpy points
rhopoints3 = ent2rho(entpoints3)

# Compute the powerlaw index of the EoS (rho ~ H^p) in each pseudo-enthalpy interval
rhoindex3 = np.zeros_like(entpoints3)
rhoindex3[0:999] = (rhopoints3[1:1000] - rhopoints3[0:999]) / (entpoints3[1:1000] - entpoints3[0:999])
rhoindex3[999] = rhoindex3[998]

enthalpy_section_valid = True
print("✓ Pseudo-enthalpy-based resampling completed successfully")

#---------------------------------------------------------------------------------
# Write the resampled tabulated EoS

npoints = 1000
output_file = eosname + "_resampled.dat"

print(f"\nWriting resampled EoS to {output_file}...")

with open(output_file, "w") as f:
    f.write(f"# Resampled EoS from {input_file}\n")
    f.write(f"# All quantities in geometrized units (G=c=1)\n")
    f.write(f"# Conversion factors: rhounit = {rhounit:.6e} g/cm³, punit = {punit:.6e} dyne/cm²\n")
    f.write(f"#\n")
    f.write(f"{npoints}\n")
    f.write("\n")
    f.write(f"{rhopoints[0]:.15f} {rhopoints[999]:.15f}\n")
    f.write(f"{prspoints[0]:.15f} {prspoints[999]:.15f} {prsindex[0]:.15f} {prsindex[999]:.15f}\n")
    f.write(f"{einpoints[0]:.15f} {einpoints[999]:.15f} {einindex[0]:.15f} {einindex[999]:.15f}\n")
    f.write(f"{entpoints[0]:.15f} {entpoints[999]:.15f} {entindex[0]:.15f} {entindex[999]:.15f}\n")
    f.write("\n")
    
    # Section 1: Uniform Log10(rho) spacing + powerlaw indexes
    f.write("# Section 1: Uniform Log10(rho) spacing\n")
    f.write("# Log10(rho)  Log10(p)  Log10(e)  Log10(H)  gamma  lambda  kappa\n")
    f.write("# where H(P) is the pseudo-enthalpy computed from H = c² ln(h/c²)\n")
    f.write("# gamma = d(log P)/d(log rho), lambda = d(log e)/d(log rho), kappa = d(log H)/d(log rho)\n")
    for i in range(1000):
        f.write(f"{rhopoints[i]:.15f}\t{prspoints[i]:.15f}\t{einpoints[i]:.15f}\t{entpoints[i]:.15f}\t{prsindex[i]:.15f}\t{einindex[i]:.15f}\t{entindex[i]:.15f}\n")
    
    # Section 2: Uniform Log10(p) spacing + powerlaw indexes
    f.write("\n")
    f.write("# Section 2: Uniform Log10(p) spacing\n")
    f.write("# Log10(p)  Log10(rho)  rho_index\n")
    for i in range(1000):
        f.write(f"{prspoints2[i]:.15f}\t{rhopoints2[i]:.15f}\t{rhoindex2[i]:.15f}\n")
    
    # Section 3: Uniform Log10(H) spacing + powerlaw indexes
    f.write("\n")
    f.write("# Section 3: Uniform Log10(H) spacing (pseudo-enthalpy based)\n")
    f.write("# Log10(H)  Log10(rho)  rho_index\n")
    f.write("# where H(P) is the pseudo-enthalpy H = c² ln(h/c²)\n")
    for i in range(1000):
        f.write(f"{entpoints3[i]:.15f}\t{rhopoints3[i]:.15f}\t{rhoindex3[i]:.15f}\n")

print(f"✓ Resampled EoS written to {output_file}")
print(f"✓ Processing complete!")

# Final summary
print(f"\nSummary:")
print(f"- Input: {len(density_cgs)} data points from {input_file}")
print(f"- Output: 3000 total points (1000 each in ρ, P, H parameterizations)")
print(f"- Density range: {10**rhopoints[0]*rhounit:.2e} to {10**rhopoints[999]*rhounit:.2e} g/cm³")
print(f"- Pressure range: {10**prspoints[0]*punit:.2e} to {10**prspoints[999]*punit:.2e} dyne/cm²")
print(f"- Pseudo-enthalpy range: {10**entpoints[0]:.2e} to {10**entpoints[999]:.2e} (dimensionless in geometrized units)")

plt.show()
