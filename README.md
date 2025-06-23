# XNS_WD: White Dwarf Modeling with XNS

_Last correction: error for the electric field_

## Overview

**XNS_WD** is a customized version of the XNS code adapted for modeling rotating, possibly magnetized white dwarfs. It supports equations of state (EoS) in both COMPOSE and Helmholtz formats, and computes solutions to the equations of hydrostatic equilibrium in General Relativity (GR) or Scalar-Tensor Theories (STT).

## Features

- Handles rotating and magnetized compact objects
- Supports tabulated EoS with resampling utilities
- Supports white dwarf modeling via Helmholtz EoS
- Flexible build system: serial, parallel, debug versions

---

## 🔧 Build System

### Main XNS Executables

```bash
make serial     # Standard XNS executable (XNS-s)
make nwtrps     # With global Newton-Raphson solver
make debug      # Debug build with extra checks
make parspace   # MPI parallel version
make clean      # Remove object files
make cleanall   # Remove all generated files
make help       # List all build targets
```

### Helmholtz EOS Generator (`_helmholtz/` directory)

```bash
cd _helmholtz
make resample   # Full workflow: build, generate, convert
make eos        # Generate white dwarf EOS table
make all        # Compile without generating tables
make clean      # Clean build files
make clean-all  # Remove all output
```

---

## ▶️ Run the Code

```bash
./XNS-s         # Run the main XNS solver
```

---

## 🧾 Key Configuration: `SYSTEMXNS.f90`

Edit this file to configure runs.

### EOS Setup

- `FILEEOS`: EOS filename (e.g., `'wd_eos_resampled.dat'`)
- `EOSINT = .TRUE.`: Use tabulated EOS

### Physical Parameters

- `RHOINI`: Central density (geometrized units)
- `OMG`: Angular velocity
- `GR`: `.TRUE.` for GR, `.FALSE.` for STT
- `IMAG`: Magnetic field switch

### Grid Parameters

- `NR`: Radial resolution (e.g., 900)
- `NTH`: Angular resolution (1 = spherical symmetry)

---

## 📊 Equation of State (EoS) Resampling Tools

Two Python scripts are available to prepare EoS data:

### 1. `resample.py` — for COMPOSE format

- Input: `eos.nb` and `eos.thermo` (baryon density, energy per baryon)
- Output: Resampled EoS in geometrized units

### 2. `resample_helmholtz.py` — for Helmholtz (CGS) format

- Input: `wd_eos_helmholtz.dat`
  - Column 1: density [g/cm³]
  - Column 2: pressure [dyn/cm²]
  - Column 3: specific internal energy [erg/g]
- Output: Resampled EoS in geometrized units

Both scripts:
- Convert data to geometrized units
- Resample with uniform spacing in log₁₀ space

### Conversion Notes

- `density = data[:,0] / rhounit`
- `pressure = data[:,1] / punit`
- `internal_energy_density = data[:,0] * data[:,2] / punit`
- `total_energy_density = internal + rest mass energy`

---

## 📁 Output Files

- `LogFile.dat`: Global properties (mass, radius)
- `Hydroeq.dat`: Radial profiles
- `Grid.dat`: Grid configuration
- `Surf.dat`: Surface values
- `TOVINIMOD_PROFILES.dat`: TOV profiles

---

## 🔄 Workflow for Custom EoS

1. Prepare EoS file (COMPOSE or Helmholtz)
2. Resample:
   ```bash
   python resample.py              # for COMPOSE
   python resample_helmholtz.py   # for Helmholtz
   ```
3. Set `FILEEOS` in `SYSTEMXNS.f90`
4. Set other simulation parameters
5. Compile: `make serial`
6. Run: `./XNS-s`

---

## 🧮 Unit System

Geometrized units (G = c = 1):
- `rhounit = Msun / (G*Msun/c²)³`
- `punit = rhounit * c²`

---

## 🧠 Code Architecture

Main source files (build order matters):

1. `SYSTEMXNS.f90`: Parameters and constants
2. `PHYSICS.f90`: EoS and physical models
3. `ROTATION.f90`: Rotation and magnetism
4. `FUNCTIONS.f90`: Numerical utilities
5. `XNS.f90`: Core solver
6. `XNSMAIN.f90`: Program entry point
7. `TOVINIMOD.f90`: TOV initial setup
8. `HYDROEQ.f90`: Hydrostatic equilibrium solver

---

## 🧪 Examples

Find sample configurations under `documentation/examples/` including:
- General Relativity and Scalar-Tensor cases
- Rotating and non-rotating models
- Magnetic and non-magnetic setups
