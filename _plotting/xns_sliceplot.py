# plot_hydroeq.py
#
# Author: Robert Fisher
# Date: 062425
# Description:
# This script visualizes XNS 4.0 Hydroeq.dat output using coordinates from Grid.dat.
# It:
# - Converts all dimensional quantities to CGS units
# - Uses log or linear scale depending on variable
# - Labels axes and units meaningfully
# - Saves a PNG image to disk
#
# Usage:
#
# python plot_hydroeq.py -c <column_index>
#
# where <column_index> is the column number to plot (1-9). 
# Default is 1 (Rest Mass Density), and the remainder of the fields are
# defined in the script.

import numpy as np
import matplotlib.pyplot as plt
import argparse

# ========= Physical constants for CGS conversion ========= #
G_cgs = 6.67430e-8             # cm^3 g^-1 s^-2
c_cgs = 2.99792458e10          # cm/s
Msun_cgs = 1.98847e33          # g

# Geometrized to CGS conversion units
Lunit = G_cgs * Msun_cgs / c_cgs**2      # length unit in cm
Tunit = G_cgs * Msun_cgs / c_cgs**3      # time unit in s
RhoUnit = Msun_cgs / Lunit**3            # density unit in g/cm^3
PUnit = RhoUnit * c_cgs**2               # pressure unit in dyn/cm^2

# ========= File paths ========= #
hydro_file = 'Hydroeq.dat'
grid_file = 'Grid.dat'

# ========= Argument Parsing ========= #
parser = argparse.ArgumentParser(description="Polar plot of XNS 4.0 Hydroeq.dat in CGS units using Grid.dat")
parser.add_argument('-c', '--column', type=int, default=1,
                    help="Which column to plot (1-9). Default is 1 (Rest Mass Density).")
args = parser.parse_args()
col_index = args.column - 1
if not (0 <= col_index < 9):
    raise ValueError("Column must be between 1 and 9.")

# ========= Read Hydroeq Header ========= #
with open(hydro_file, 'r') as f:
    header_line = f.readline()
NTH, NR, OMG = map(float, header_line.split())
NTH, NR = int(NTH), int(NR)
OMG = float(OMG)

# ========= Read Grid.dat ========= #
with open(grid_file, 'r') as f:
    grid_header = f.readline()
    nth_check, nr_check, _, _ = map(float, grid_header.split())
    nth_check, nr_check = int(nth_check), int(nr_check)
    if nth_check != NTH or nr_check != NR:
        raise ValueError("Grid dimensions in Grid.dat do not match Hydroeq.dat.")
    theta_vals = np.array([float(f.readline()) for _ in range(NTH)])
    radius_vals_geom = np.array([float(f.readline()) for _ in range(NR)])
    radius_vals = radius_vals_geom * Lunit  # convert to cm

# ========= Load and reshape data ========= #
data_flat = np.loadtxt(hydro_file, skiprows=1, usecols=col_index)
if data_flat.size != NTH * NR:
    raise ValueError("Mismatch between data size and grid dimensions.")
data_grid = data_flat.reshape((NTH, NR))

# ========= Field names and scale settings ========= #
field_names = {
    1: ("Rest Mass Density", RhoUnit, "g/cm^3", "log"),
    2: ("Pressure", PUnit, "dyn/cm^2", "log"),
    3: ("Conformal Factor (psi)", 1.0, "", "linear"),
    4: ("Normalized Azimuthal Speed (|v^φ| / α)", 1.0, "", "linear"),
    5: ("Lapse Function (alpha = PSL/psi)", 1.0, "", "log"),
    6: ("Shift Vector (beta^phi)", 1.0 / Tunit, "rad/s", "linear"),
    7: ("Lapse Numerator (chi)", 1.0, "", "linear"),
    8: ("Quadrupole Scalar (Radial)", 1.0, "", "log"),
    9: ("Quadrupole Scalar (Theta)", 1.0, "", "log")
}
field_label, unit_scale, unit_str, scale_type = field_names.get(
    args.column, (f"Field #{args.column}", 1.0, "", "log")
)

# ========= Convert and scale data ========= #
data_cgs = data_grid * unit_scale
if scale_type == "log":
    plot_data = np.where(data_cgs > 0, np.log10(data_cgs), np.nan)
    colorbar_label = f"log10({field_label} [{unit_str}])" if unit_str else f"log10({field_label})"
else:
    plot_data = data_cgs
    colorbar_label = f"{field_label} [{unit_str}]" if unit_str else field_label

# ========= Coordinate Grid ========= #
Theta, Radius = np.meshgrid(theta_vals, radius_vals, indexing='ij')

# ========= Plotting ========= #
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 6))
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_thetamin(0)
ax.set_thetamax(180)

c = ax.pcolormesh(Theta, Radius, plot_data, cmap='seismic', shading='auto')
cbar = fig.colorbar(c, ax=ax, pad=0.1)
cbar.set_label(colorbar_label)

plot_title = f"XNS 4.0 Plot – Omega = {OMG:.3g} rad/s\n{field_label}"
ax.set_title(plot_title, fontsize=13)
plt.tight_layout()

filename_safe = field_label.lower().replace(" ", "_").replace("(", "").replace(")", "").replace("^", "").replace("=", "").replace("/", "").replace(".", "")
outfile = f"xns_plot_{filename_safe}.png"
plt.savefig(outfile, dpi=300)
print(f"Plot saved as: {outfile}")

plt.show()

