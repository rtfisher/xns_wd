#!/usr/bin/env python3

"""
Script: generate_wd_rotating_sequence.py

This script automates the generation of rotating white dwarf models using the XNS code,
constructing a sequence of fixed-mass models by incrementally increasing the rotation rate
(Omega) from zero up to the breakup limit. The central density is iteratively adjusted using
bisection to maintain the target gravitational mass within a tolerance.

Key features:
- Edits SYSTEMXNS.f90 to set RHOINI and Omega
- Compiles and runs XNS
- Parses the output mass and refines RHOINI
- Stores output files tagged by central density and Omega

Sample usage:
$ python3 generate_wd_rotating_sequence.py
"""

import os
import shutil
import subprocess
import re
import numpy as np
from pathlib import Path
from datetime import datetime

# ======================== USER CONFIGURATION ========================
XNS_DIR = Path("/Users/rfisher1/Dropbox/Research/Programs/_xns_model_sequence")  # Default XNS directory
N_OMEGA = 20                             # Number of Omega steps
OMEGA_MAX = 2.3e-5                       # Approximate breakup
MASS_TOL = 0.001                         # Tolerance in mass (0.1%)
MAX_ITER = 100                           # Max iterations for bisection
LOG_TO_STDOUT = True                     # Whether to print log to stdout
MODEL_DIR = Path("_model_directory")    # Output base directory
RHO_C_LIST = [2.0e9]                     # List of central densities in g/cm^3 to loop over

# ======================== GLOBAL VARIABLES ==========================
LOGFILE = None

# ======================== UTILITY FUNCTIONS =========================
def log(msg):
    timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    line = f"{timestamp} {msg}"
    if LOG_TO_STDOUT:
        print(line)
    if LOGFILE:
        LOGFILE.write(line + "\n")
        LOGFILE.flush()

def edit_systemxns(rho_ini, omega):
    """Edit SYSTEMXNS.f90 to set RHOINI and OMG values with robust pattern matching"""
    sys_path = XNS_DIR / "SYSTEMXNS.f90"
    
    # Read the file
    with open(sys_path, "r") as f:
        lines = f.readlines()
    
    new_lines = []
    rhoini_updated = False
    omg_updated = False
    
    for line in lines:
        stripped = line.strip()
        original_line = line
        
        # More robust RHOINI pattern matching
        # Look for lines containing RHOINI that aren't comments
        if not stripped.startswith("!") and "RHOINI" in line:
            # Try multiple possible patterns
            if re.search(r"REAL\s*,?\s*PARAMETER\s*::\s*RHOINI", line, re.IGNORECASE):
                newline = f"  REAL,PARAMETER :: RHOINI = {rho_ini:.6e} / RHO_CGS_TO_GEOM ! Central density in geom units\n"
                rhoini_updated = True
                log(f"DIAGNOSTIC: Updated RHOINI line: {newline.strip()}")
            elif re.search(r"RHOINI\s*=", line):
                # Alternative pattern: just RHOINI = ...
                newline = f"  RHOINI = {rho_ini:.6e} / RHO_CGS_TO_GEOM ! Central density in geom units\n"
                rhoini_updated = True
                log(f"DIAGNOSTIC: Updated RHOINI assignment: {newline.strip()}")
            else:
                newline = original_line
        # More robust OMG pattern matching  
        elif not stripped.startswith("!") and "OMG" in line and re.search(r"REAL\s*::\s*OMG", line, re.IGNORECASE):
            newline = f"  REAL           :: OMG = {omega:.7e} ! Central Rotation rate\n"
            omg_updated = True
            log(f"DIAGNOSTIC: Updated OMG line: {newline.strip()}")
        else:
            newline = original_line
            
        new_lines.append(newline)
    
    # Write the file
    with open(sys_path, "w") as f:
        f.writelines(new_lines)
    
    # Verify the updates
    if not rhoini_updated:
        log("ERROR: RHOINI line not found or updated in SYSTEMXNS.f90")
        log("DIAGNOSTIC: Searching for RHOINI patterns in file:")
        with open(sys_path, "r") as f:
            for i, line in enumerate(f, 1):
                if "RHOINI" in line:
                    log(f"DIAGNOSTIC: Line {i}: {line.strip()}")
    
    if not omg_updated:
        log("ERROR: OMG line not found or updated in SYSTEMXNS.f90")
        log("DIAGNOSTIC: Searching for OMG patterns in file:")
        with open(sys_path, "r") as f:
            for i, line in enumerate(f, 1):
                if "OMG" in line and "REAL" in line:
                    log(f"DIAGNOSTIC: Line {i}: {line.strip()}")
    
    log(f"DIAGNOSTIC: File edit complete. RHOINI updated: {rhoini_updated}, OMG updated: {omg_updated}")

def compile_xns():
    os.chdir(XNS_DIR)
    subprocess.run(["make", "clean"], stdout=subprocess.DEVNULL)
    result = subprocess.run(["make", "serial"], capture_output=True)
    if result.returncode != 0 or not (XNS_DIR / "XNS-s").exists():
        raise RuntimeError("Compilation failed:\n" + result.stderr.decode())

def run_xns():
    """Run XNS with better error checking and diagnostics"""
    log("DIAGNOSTIC: Starting XNS execution")
    
    # Ensure we're in the right directory
    original_dir = os.getcwd()
    os.chdir(XNS_DIR)
    
    try:
        # Check if executable exists
        xns_exe = XNS_DIR / "XNS-s"
        if not xns_exe.exists():
            raise RuntimeError("XNS-s executable not found")
        
        # Run XNS and capture all output
        result = subprocess.run(["./XNS-s"], 
                              capture_output=True, 
                              text=True, 
                              timeout=60)  # 60 second timeout
        
        log(f"DIAGNOSTIC: XNS return code: {result.returncode}")
        
        if result.stdout:
            log(f"DIAGNOSTIC: XNS stdout (first 500 chars): {result.stdout[:500]}")
        if result.stderr:
            log(f"DIAGNOSTIC: XNS stderr: {result.stderr}")
        
        if result.returncode != 0:
            raise RuntimeError(f"XNS execution failed with return code {result.returncode}")
            
        # Check if LogFile.dat was created/updated
        logfile = XNS_DIR / "LogFile.dat"
        if not logfile.exists():
            raise RuntimeError("LogFile.dat was not created by XNS")
            
        log("DIAGNOSTIC: XNS execution completed successfully")
        
    except subprocess.TimeoutExpired:
        raise RuntimeError("XNS execution timed out after 60 seconds")
    finally:
        os.chdir(original_dir)

def parse_mass():
    with open(XNS_DIR / "LogFile.dat", "r") as f:
        for line in f:
            if "REST    MASS" in line:
                val = line.split("=")[-1].strip()
                if "NaN" in val:
                    return float("nan")
                return float(val)
    raise ValueError("Could not find REST MASS in LogFile.dat")

def store_output_files(mass_solar, rho_c, omega):
    subdir = MODEL_DIR / f"_{mass_solar:.5f}"
    subdir.mkdir(parents=True, exist_ok=True)
    rho_tag = rho_c * 1e-9  # Convert to 10^9 g/cm^3 units
    for fname in ["Grid.dat", "Hydroeq.dat", "LogFile.dat", "Rhovec.dat", "Surf.dat", "TOVINIMOD_DEBUG.dat", "TOVINIMOD_PROFILES.dat"]:
        src = XNS_DIR / fname
        if src.exists():
            tagged = f"{src.stem}_{rho_tag:.3f}_{omega:.6e}.dat"
            shutil.copy(src, subdir / tagged)

# ======================== MODIFIED BISECTION METHOD =========================

def find_rho_for_mass(target_mass, omega, rho_initial):
    """
    Find central density that produces target_mass for given omega.
    Returns (rho_converged, mass_converged) or (None, None) if failed.
    """
    
    def evaluate_mass(rho):
        """Helper to run model and get mass, with error handling and extensive diagnostics"""
        try:
            # Diagnostic 1: Verify the edit actually happened
            log(f"=== DIAGNOSTIC: Attempting to set RHOINI = {rho:.6e} ===")
            edit_systemxns(rho, omega)
            
            # Diagnostic 2: Verify the file was actually modified
            sys_path = XNS_DIR / "SYSTEMXNS.f90"
            with open(sys_path, "r") as f:
                content = f.read()
                rho_match = re.search(r'RHOINI\s*=\s*([+-]?\d+\.?\d*[eE]?[+-]?\d*)', content)
                omega_match = re.search(r'OMG\s*=\s*([+-]?\d+\.?\d*[eE]?[+-]?\d*)', content)
                if rho_match:
                    file_rho = float(rho_match.group(1))
                    log(f"DIAGNOSTIC: SYSTEMXNS.f90 contains RHOINI = {file_rho:.6e}")
                    if abs(file_rho - rho) > rho * 1e-10:
                        log(f"ERROR: File rho ({file_rho:.6e}) != requested rho ({rho:.6e})")
                else:
                    log("ERROR: Could not find RHOINI in SYSTEMXNS.f90")
                if omega_match:
                    file_omega = float(omega_match.group(1))
                    log(f"DIAGNOSTIC: SYSTEMXNS.f90 contains OMG = {file_omega:.6e}")
                
            # Diagnostic 3: Check compilation timestamp
            xns_executable = XNS_DIR / "XNS-s"
            if xns_executable.exists():
                import os
                before_time = os.path.getmtime(xns_executable)
                log(f"DIAGNOSTIC: XNS-s timestamp before compilation: {before_time}")
            
            compile_xns()
            
            if xns_executable.exists():
                after_time = os.path.getmtime(xns_executable)
                log(f"DIAGNOSTIC: XNS-s timestamp after compilation: {after_time}")
                if before_time == after_time:
                    log("WARNING: XNS-s timestamp unchanged - compilation may have failed silently")
            
            # Diagnostic 4: Check if LogFile.dat exists and get timestamp before run
            logfile_path = XNS_DIR / "LogFile.dat"
            if logfile_path.exists():
                before_log_time = os.path.getmtime(logfile_path)
                log(f"DIAGNOSTIC: LogFile.dat timestamp before run: {before_log_time}")
            else:
                before_log_time = None
                log("DIAGNOSTIC: LogFile.dat does not exist before run")
            
            run_xns()
            
            # Diagnostic 5: Check if LogFile.dat was updated
            if logfile_path.exists():
                after_log_time = os.path.getmtime(logfile_path)
                log(f"DIAGNOSTIC: LogFile.dat timestamp after run: {after_log_time}")
                if before_log_time and before_log_time == after_log_time:
                    log("ERROR: LogFile.dat timestamp unchanged - XNS may not have run properly")
            else:
                log("ERROR: LogFile.dat does not exist after XNS run")
                return None
            
            # Diagnostic 6: Show multiple lines from LogFile.dat around the mass
            log("DIAGNOSTIC: LogFile.dat content around REST MASS:")
            with open(logfile_path, "r") as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if "REST    MASS" in line:
                        # Show context: 2 lines before, the mass line, and 2 lines after
                        start = max(0, i-2)
                        end = min(len(lines), i+3)
                        for j in range(start, end):
                            marker = ">>> " if j == i else "    "
                            log(f"DIAGNOSTIC: {marker}{lines[j].rstrip()}")
                        break
            
            mass = parse_mass()
            log(f"DIAGNOSTIC: Parsed mass = {mass:.15f} Msol")
            
            # Diagnostic 7: Check for identical masses across calls
            if hasattr(evaluate_mass, '_last_mass') and evaluate_mass._last_mass is not None:
                if abs(mass - evaluate_mass._last_mass) < 1e-14:
                    log(f"WARNING: Identical mass to previous call: {mass:.15f}")
                    log(f"Previous mass: {evaluate_mass._last_mass:.15f}")
                    log(f"Difference: {abs(mass - evaluate_mass._last_mass):.2e}")
            evaluate_mass._last_mass = mass
            
            log(f"=== END DIAGNOSTIC ===")
            return mass
            
        except Exception as e:
            log(f"Model failed at RHOINI = {rho:.6e}: {e}")
            import traceback
            log(f"Full traceback: {traceback.format_exc()}")
            return None
    
    # Step 1: Find bracketing interval
    log("Finding bracketing interval...")
    rho_high = rho_initial
    rho_low = rho_initial * 0.5
    
    mass_high = evaluate_mass(rho_high)
    if mass_high is None:
        return None, None
        
    mass_low = evaluate_mass(rho_low)
    if mass_low is None:
        return None, None
    
    # Determine mass-density relationship from initial points
    mass_increases_with_rho = mass_high > mass_low
    log(f"Mass-density relationship: mass {'increases' if mass_increases_with_rho else 'decreases'} with density")
    
    # Expand bracket if needed - keep expanding until we bracket the target
    max_bracket_attempts = 20
    for attempt in range(max_bracket_attempts):
        # Check if we have a valid bracket
        if (mass_low - target_mass) * (mass_high - target_mass) < 0:
            log(f"Bracket found: [{rho_low:.6e}, {rho_high:.6e}] -> masses [{mass_low:.6f}, {mass_high:.6f}]")
            break
        
        # Special case: if masses are essentially identical (within numerical precision), 
        # we've hit the limit where XNS gives the same result regardless of input density
        mass_diff = abs(mass_low - mass_high)
        relative_diff = mass_diff / target_mass
        
        if mass_diff < 1e-12:  # Masses are numerically identical
            log(f"Masses are numerically identical: {mass_low:.15f} ≈ {mass_high:.15f}")
            log(f"Mass difference: {mass_diff:.2e}")
            log(f"Target mass: {target_mass:.15f}")
            log(f"Difference from target: {abs(mass_low - target_mass):.15f}")
            log(f"Relative error: {abs(mass_low - target_mass)/target_mass:.8f} ({abs(mass_low - target_mass)/target_mass*100:.6f}%)")
            
            if abs(mass_low - target_mass) / target_mass < MASS_TOL:
                log("Target is within tolerance - accepting result")
                return rho_low, mass_low  # Return either point since they're identical
            else:
                log("Target is outside tolerance but no bracket possible - this may be a physical/numerical limit")
                log("Consider increasing MASS_TOL or investigating the physics at this rotation rate")
                return None, None
        
        # If masses are different but relative difference is very small, warn about potential precision issues
        if relative_diff < MASS_TOL / 10:
            log(f"Warning: Mass range [{mass_low:.15f}, {mass_high:.15f}] is very narrow")
            log(f"Relative span: {relative_diff:.2e} (much smaller than tolerance {MASS_TOL:.2e})")
            log("This may indicate numerical precision limits")
            
        log(f"Bracket attempt {attempt+1}: target={target_mass:.6f}, current range=[{mass_low:.6f}, {mass_high:.6f}]")
        
        # Determine which direction to expand based on where target lies and mass-density relationship
        if target_mass < mass_low:
            # Need lower mass than current minimum
            if mass_increases_with_rho:
                # Mass increases with rho, so decrease rho_low to get lower mass
                rho_low *= 0.7
                new_mass = evaluate_mass(rho_low)
                if new_mass is None:
                    return None, None
                mass_low = new_mass
                log(f"Expanded low (decreasing rho): rho_low={rho_low:.6e} -> mass={mass_low:.6f}")
            else:
                # Mass decreases with rho, so increase rho_high to get lower mass
                rho_high *= 1.3
                new_mass = evaluate_mass(rho_high)
                if new_mass is None:
                    return None, None
                mass_high = new_mass
                log(f"Expanded high (increasing rho): rho_high={rho_high:.6e} -> mass={mass_high:.6f}")
                
        elif target_mass > mass_high:
            # Need higher mass than current maximum
            if mass_increases_with_rho:
                # Mass increases with rho, so increase rho_high to get higher mass
                rho_high *= 1.3
                new_mass = evaluate_mass(rho_high)
                if new_mass is None:
                    return None, None
                mass_high = new_mass
                log(f"Expanded high (increasing rho): rho_high={rho_high:.6e} -> mass={mass_high:.6f}")
            else:
                # Mass decreases with rho, so decrease rho_low to get higher mass
                rho_low *= 0.7
                new_mass = evaluate_mass(rho_low)
                if new_mass is None:
                    return None, None
                mass_low = new_mass
                log(f"Expanded low (decreasing rho): rho_low={rho_low:.6e} -> mass={mass_low:.6f}")
        else:
            # Target is between mass_low and mass_high but bracket condition failed
            # This should not happen, but handle edge case
            log(f"Edge case: target {target_mass:.6f} is between [{mass_low:.6f}, {mass_high:.6f}] but no bracket")
            # Try small expansion in both directions
            old_rho_low, old_mass_low = rho_low, mass_low
            rho_low *= 0.95
            mass_low = evaluate_mass(rho_low)
            if mass_low is None:
                rho_low, mass_low = old_rho_low, old_mass_low
                rho_high *= 1.05
                mass_high = evaluate_mass(rho_high)
                if mass_high is None:
                    return None, None
        
        # Update mass-density relationship if we've flipped ordering
        if (mass_high > mass_low) != mass_increases_with_rho:
            mass_increases_with_rho = mass_high > mass_low
            log(f"Mass-density relationship changed: mass now {'increases' if mass_increases_with_rho else 'decreases'} with density")
    else:
        log("Failed to find bracketing interval after maximum attempts")
        log(f"Final range: rho=[{rho_low:.6e}, {rho_high:.6e}], mass=[{mass_low:.6f}, {mass_high:.6f}], target={target_mass:.6f}")
        return None, None
    
    # Step 2: Bisection refinement
    log("Refining with bisection...")
    for iteration in range(MAX_ITER):
        rho_mid = 0.5 * (rho_low + rho_high)
        mass_mid = evaluate_mass(rho_mid)
        
        # Handle NaN or None values from failed XNS runs
        # This can happen when the density is too extreme (too high/low) or 
        # when XNS fails to converge to a physical solution
        if mass_mid is None or np.isnan(mass_mid):
            log(f"Model failed at midpoint rho={rho_mid:.6e}, shrinking interval")
            
            # Strategy: Instead of using the exact midpoint (which failed), 
            # try quarter points within the bracket to find a working density.
            # This helps avoid "bad" regions in parameter space while maintaining
            # the bracket integrity.
            
            # Quarter point closer to rho_high (75% toward high end)
            rho_quarter = 0.25 * rho_low + 0.75 * rho_high
            # Quarter point closer to rho_low (25% toward high end)  
            rho_three_quarter = 0.75 * rho_low + 0.25 * rho_high
            
            # Try the quarter point first (often safer than midpoint)
            mass_quarter = evaluate_mass(rho_quarter)
            if mass_quarter is not None and not np.isnan(mass_quarter):
                # Found a working point - use it as our new evaluation point
                rho_mid = rho_quarter
                mass_mid = mass_quarter
                log(f"Using quarter point: rho={rho_mid:.6e}, mass={mass_mid:.10f}")
            else:
                # Quarter point also failed, try the three-quarter point
                mass_three_quarter = evaluate_mass(rho_three_quarter)
                if mass_three_quarter is not None and not np.isnan(mass_three_quarter):
                    # Found a working point at three-quarter position
                    rho_mid = rho_three_quarter
                    mass_mid = mass_three_quarter
                    log(f"Using three-quarter point: rho={rho_mid:.6e}, mass={mass_mid:.10f}")
                else:
                    # Both alternative points failed - this suggests the entire
                    # interval may be problematic. Stop the bisection.
                    log("Both quarter points failed, stopping bisection")
                    log("This suggests numerical instability in the current density range")
                    return None, None
            
        log(f"Iteration {iteration+1}: RHOINI = {rho_mid:.6e}, mass = {mass_mid:.10f}")
        
        # Check convergence
        if abs(mass_mid - target_mass) / target_mass < MASS_TOL:
            log(f"Converged! Final mass = {mass_mid:.10f} Msol")
            return rho_mid, mass_mid
        
        # Update bracket based on which side of target the current mass lies
        # This preserves the bracket property: target is always between mass_low and mass_high
        if (mass_mid - target_mass) * (mass_low - target_mass) < 0:
            # mass_mid and mass_low are on opposite sides of target
            # Replace the high end with our new point
            rho_high = rho_mid
            mass_high = mass_mid
        else:
            # mass_mid and mass_high are on opposite sides of target  
            # Replace the low end with our new point
            rho_low = rho_mid
            mass_low = mass_mid
            
        # Safety check: prevent infinite loops due to numerical precision limits
        # If the density bracket becomes extremely small, we've hit the limits
        # of floating-point precision and should stop
        if abs(rho_high - rho_low) / rho_low < 1e-12:
            log(f"Bracket too small: [{rho_low:.6e}, {rho_high:.6e}], stopping")
            log(f"Best mass achieved: {mass_mid:.10f}, target: {target_mass:.10f}")
            log("This indicates we've reached numerical precision limits")
            # Return the best result even if not within tolerance
            return rho_mid, mass_mid
    
    log("Failed to converge within max iterations")
    log(f"Final: rho={rho_mid:.6e}, mass={mass_mid:.10f}, target={target_mass:.10f}")
    return None, None

# ======================== MAIN DRIVER ===============================

def main():
    global LOGFILE
    logname = Path("model_sequence.log")
    LOGFILE = open(logname, "w")

    for rho_c in RHO_C_LIST:
        log(f"==== Starting outer loop for rho_c = {rho_c:.2e} g/cm^3 ====")
        omegas = [OMEGA_MAX * i / (N_OMEGA - 1) for i in range(N_OMEGA)]
        target_mass = None
        last_rho = rho_c  # Track last successful density for better initial guess

        for i, omega in enumerate(omegas):
            log(f"==== Starting model {i+1}/{N_OMEGA}, Omega = {omega:.6e} ====")

            if i == 0:
                # Non-rotating case - use initial rho_c
                edit_systemxns(rho_c, omega)
                compile_xns()
                run_xns()
                target_mass = parse_mass()
                log(f"Target rest mass (Omega=0): {target_mass:.10f} Msol")
                store_output_files(target_mass, rho_c, omega)
                last_rho = rho_c  # Initialize for next iteration
                continue

            # Find rho that gives target mass at this omega
            rho_converged, mass_converged = find_rho_for_mass(target_mass, omega, last_rho)
            
            if rho_converged is None:
                log("Failed to find solution; stopping sequence")
                break

            store_output_files(target_mass, rho_converged, omega)
            last_rho = rho_converged  # Use as starting guess for next omega

    LOGFILE.close()

if __name__ == "__main__":
    main()
