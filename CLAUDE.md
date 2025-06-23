# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# XNS Code Overview

## What is XNS?
XNS is a numerical code for computing neutron star models in both General Relativity and Scalar-Tensor theories. It solves the hydrostatic equilibrium equations for rotating, magnetized neutron stars.

## Build System

### Main XNS Code:
```bash
make serial     # Standard XNS executable (XNS-s)
make nwtrps     # XNS with global Newton-Raphson (XNS-nr)
make debug      # Debug version with extensive checks
make parspace   # MPI parallel version (XNS-mpi)
make clean      # Remove object files and executables
make cleanall   # Remove all generated files including .dat files
make help       # Show all available targets
```

### Helmholtz EOS Generator (_helmholtz/ directory):
```bash
cd _helmholtz
make resample   # Complete workflow: build, generate EOS, convert to XNS format
make eos        # Generate white dwarf EOS table only
make all        # Build executable only
make clean      # Remove build files
make clean-all  # Remove all generated files
```

### Run the code:
```bash
./XNS-s         # Main executable
```

## Key Configuration Files

### SYSTEMXNS.f90
Main configuration file containing:

**EOS Configuration:**
- `FILEEOS`: EOS filename (e.g., 'APR_resampled.dat')
- `EOSINT = .TRUE.`: Use tabulated EOS

**Physical Parameters:**
- `RHOINI`: Central density in geometrized units
- `OMG`: Rotation rate
- `GR = .TRUE.`: General Relativity (.FALSE. for Scalar-Tensor)

**Grid Parameters:**
- `NR = 900`: Radial grid points
- `NTH = 100`: Angular grid points (1 = spherical symmetry)

**Physics Options:**
- `IMAG = .FALSE.`: Magnetic fields on/off

## EoS File Preparation

### Two resampling tools available:

1. **resample.py**: For COMPOSE library format
   - Reads: `eos.nb` (baryon density) and `eos.thermo` (thermodynamic data)
   - COMPOSE format uses baryon number density and energy per baryon

2. **resample_helmholtz.py**: For CGS format files
   - Reads: files with density [g/cm³], pressure [dyn/cm²], specific internal energy [erg/g]
   - Example: `wd_eos_helmholtz.dat`

### EoS Data Conversion:
- **COMPOSE**: Uses baryon density and energy per baryon, converts via nuclear physics constants
- **Helmholtz**: Direct CGS units conversion
  - Density: `data[:,0] / rhounit`
  - Pressure: `data[:,1] / punit`
  - Internal energy density: `data[:,0] * data[:,2] / punit`
  - Total energy density: internal + rest mass energy

### Key Energy Distinctions:
- **Specific internal energy**: Energy per unit mass [erg/g]
- **Internal energy density**: Energy per unit volume [erg/cm³]
- **Total energy density**: Internal + rest mass energy density

## Output Files

- `LogFile.dat`: Summary of stellar properties (mass, radius, etc.)
- `Hydroeq.dat`: Radial profiles of physical quantities
- `Grid.dat`: Grid information
- `Surf.dat`: Surface data
- `TOVINIMOD_PROFILES.dat`: TOV profiles

## Workflow for New EoS

1. Prepare EoS data in appropriate format
2. Run resampling script: `python resample_helmholtz.py` (for CGS data)
3. Edit `SYSTEMXNS.f90` to set `FILEEOS` to your resampled file
4. Set other parameters (central density, physics options)
5. Recompile: `make serial`
6. Run: `./XNS-s`

## Unit System
Code uses geometrized units where G = c = 1:
- `rhounit = Msun/(G*Msun/c²)³`
- `punit = rhounit * c²`

## Code Architecture

### Main Source Files (compilation order matters):
1. **SYSTEMXNS.f90**: Global parameters, constants, and configuration
2. **PHYSICS.f90**: Physical models and EOS handling
3. **ROTATION.f90**: Rotation and magnetic field computations
4. **FUNCTIONS.f90**: Mathematical utilities and numerical methods
5. **XNS.f90**: Core stellar structure solver
6. **XNSMAIN.f90**: Main program entry point
7. **TOVINIMOD.f90**: TOV (Tolman-Oppenheimer-Volkoff) initial conditions
8. **HYDROEQ.f90**: Hydrostatic equilibrium solver

### Development Pattern:
- Configuration changes: Edit `SYSTEMXNS.f90` 
- Physics modifications: Edit `PHYSICS.f90`
- New numerical methods: Add to `FUNCTIONS.f90`
- Main solver changes: Edit `XNS.f90` or `HYDROEQ.f90`
- Always recompile after any source changes: `make clean && make serial`

### EOS Integration:
- Use `_helmholtz/` directory for white dwarf EOS generation
- Use `resample_helmholtz.py` for CGS format EOS files
- Use `resample.py` for COMPOSE library format EOS files
- Generated EOS files must be resampled before use with XNS

## Examples
Check `documentation/examples/` for various configurations (GR/STT, magnetic/non-magnetic, rotating/non-rotating).