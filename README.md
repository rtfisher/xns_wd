Last version of XNS as of 13/04/2022.
Last correction: error ffor the electric filed

## EoS Resampling Tools

Two Python scripts are provided for resampling equation of state data:

- `resample.py`: Original script for COMPOSE library format (reads `eos.nb` and `eos.thermo` files)
- `resample_helmholtz.py`: Modified script for Helmholtz EoS format (reads `wd_eos_helmholtz.dat` file)

The Helmholtz version reads EoS data in CGS units:
- Column 1: density [g/cm³]
- Column 2: pressure [dyn/cm²] 
- Column 3: specific internal energy [erg/g]

Both scripts convert the data to geometrized units and resample with uniform spacing in log₁₀ scale.
