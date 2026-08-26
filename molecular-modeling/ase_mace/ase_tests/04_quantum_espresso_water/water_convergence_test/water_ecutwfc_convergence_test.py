#
# Convergence test for water molecule with Quantum ESPRESSO
# Tests different ecutwfc values using static (single point) calculations
#
import os
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.build import molecule
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.units import Ry, Bohr
from datetime import datetime

# Set OpenMP threads to 1
os.environ['OMP_NUM_THREADS'] = '1'

# Set the path to your pseudopotential directory
PSEUDO_DIR = '/usr/share/espresso/pseudo/'

# Create the water molecule with a fixed geometry
h2o = molecule('H2O')
h2o.set_cell([15, 15, 15])  # Larger vacuum box for convergence test
h2o.center()

# Fix the geometry for all calculations (no optimization)
# The molecule will stay at its initial geometry

# Pseudopotentials
pseudopotentials = {
    'O': 'O.pbe-kjpaw.UPF',
    'H': 'H.pbe-kjpaw.UPF'
}

# Base input data (will be modified for each calculation)
base_input_data = {
    'control': {
        'calculation': 'scf',  # Single point calculation
        'prefix': 'h2o_conv',
        'outdir': './conv_outdir',
        'verbosity': 'low',
        'tstress': True,
        'tprnfor': True
    },
    'system': {
        'ecutwfc': 46.0,  # Will be varied
        'ecutrho': 200.0,  # Will be set based on ecutwfc
        'ibrav': 0,
        'nosym': True,
        'noinv': True,
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.02,
    },
    'electrons': {
        'conv_thr': 1e-8,
        'mixing_beta': 0.7,
        'diagonalization': 'david'
    }
}

# Command to run pw.x
command = 'mpirun -np 4 pw.x'
profile = EspressoProfile(command, pseudo_dir=PSEUDO_DIR)


def run_single_point(atoms, ecutwfc, ecutrho=None):
    """
    Run a single point (static) calculation with given cutoffs
    
    Parameters:
    -----------
    atoms : ASE Atoms object
        The atoms to calculate (geometry is fixed)
    ecutwfc : float
        Wavefunction cutoff in Ry
    ecutrho : float or None
        Charge density cutoff in Ry. If None, uses 4 * ecutwfc
    
    Returns:
    --------
    energy : float
        Total energy in eV
    """
    # Create a copy of the atoms
    atoms_copy = atoms.copy()
    
    # Set the cutoffs
    if ecutrho is None:
        ecutrho = 4.0 * ecutwfc
    
    # Create input data with current cutoffs
    input_data = base_input_data.copy()
    input_data['system']['ecutwfc'] = ecutwfc
    input_data['system']['ecutrho'] = ecutrho
    
    # Create calculator
    calc = Espresso(profile=profile,
                   pseudopotentials=pseudopotentials,
                   input_data=input_data,
                   kpts=(1, 1, 1))
    
    atoms_copy.calc = calc
    
    # Run single point calculation
    try:
        energy = atoms_copy.get_potential_energy()
        return energy
    except Exception as e:
        print(f"Error at ecutwfc = {ecutwfc} Ry: {e}")
        return None


def convergence_test_ecutwfc(atoms, cutoff_range, ecutrho_factor=4.0):
    """
    Test convergence with respect to ecutwfc using static calculations
    
    Parameters:
    -----------
    atoms : ASE Atoms object
        The atoms to test (fixed geometry)
    cutoff_range : list or array
        List of ecutwfc values to test
    ecutrho_factor : float
        Multiplier for ecutrho (ecutrho = ecutrho_factor * ecutwfc)
    
    Returns:
    --------
    results : dict
        Dictionary containing energies and details
    """
    print("\n" + "="*70)
    print(f"CONVERGENCE TEST: ecutwfc variation (Static calculations)")
    print(f"ecutrho = {ecutrho_factor} * ecutwfc")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70)
    print(f"{'ecutwfc (Ry)':>12} | {'ecutrho (Ry)':>12} | {'Energy (eV)':>15} | ΔE (meV)")
    print("-"*70)
    
    energies = []
    successful_cutoffs = []
    
    for ecut in cutoff_range:
        ecutrho = ecutrho_factor * ecut
        
        # Run single point calculation
        energy = run_single_point(atoms, ecut, ecutrho)
        
        if energy is not None:
            energies.append(energy)
            successful_cutoffs.append(ecut)
            
            # Calculate energy difference from previous
            if len(energies) > 1:
                delta_e = (energy - energies[-2]) * 1000  # Convert to meV
                print(f"{ecut:12.1f} | {ecutrho:12.1f} | {energy:15.6f} | {delta_e:8.3f} meV")
            else:
                print(f"{ecut:12.1f} | {ecutrho:12.1f} | {energy:15.6f} | {'--':>8}")
        else:
            print(f"{ecut:12.1f} | {ecutrho:12.1f} | {'FAILED':>15} | {'--':>8}")
    
    print("="*70)
    
    results = {
        'ecutwfc': successful_cutoffs,
        'energies': energies
    }
    
    return results


def plot_convergence_results(ecut_results):
    """
    Plot the convergence test results
    
    Parameters:
    -----------
    ecut_results : dict
        Results from ecutwfc convergence test
    """
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot ecutwfc convergence
    ecut_values = ecut_results['ecutwfc']
    energies = ecut_results['energies']
    
    # Convert to meV relative to highest cutoff
    energy_mev = [(e - energies[-1]) * 1000 for e in energies]
    
    ax.plot(ecut_values, energy_mev, 'bo-', linewidth=2, markersize=8)
    ax.axhline(y=1.0, color='r', linestyle='--', label='1 meV threshold', linewidth=2)
    ax.axhline(y=0.0, color='k', linestyle='-', linewidth=0.5)
    ax.set_xlabel('ecutwfc (Ry)', fontsize=14)
    ax.set_ylabel('ΔE (meV)', fontsize=14)
    ax.set_title('Convergence of Total Energy with ecutwfc (Static)', fontsize=16)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=12)
    
    # Add value labels
    for i, (x, y) in enumerate(zip(ecut_values, energy_mev)):
        ax.annotate(f'{y:.1f}', (x, y), textcoords="offset points", 
                    xytext=(0, 10), ha='center', fontsize=10)
    
    plt.tight_layout()
    plt.savefig('convergence_test_ecutwfc.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("\nConvergence plot saved as 'convergence_test_ecutwfc.png'")


def print_summary(ecut_results):
    """
    Print a summary of the convergence test results with recommendations
    """
    print("\n" + "="*70)
    print("CONVERGENCE TEST SUMMARY")
    print("="*70)
    
    # ecutwfc summary
    ecut_values = ecut_results['ecutwfc']
    energies = ecut_results['energies']
    
    print("\n📊 ecutwfc Convergence (Static Calculations):")
    print("-"*50)
    print(f"  Minimum tested: {ecut_values[0]:.0f} Ry")
    print(f"  Maximum tested: {ecut_values[-1]:.0f} Ry")
    print(f"  Total energy range: {(energies[-1] - energies[0]) * 1000:.2f} meV")
    
    # Find convergence (within 1 meV)
    if len(energies) > 1:
        energy_diff = [abs(energies[i] - energies[-1]) * 1000 for i in range(len(energies))]
        
        # Find first point within 1 meV
        converged_idx = None
        for i, diff in enumerate(energy_diff):
            if diff < 1.0:
                converged_idx = i
                break
        
        if converged_idx is not None:
            print(f"\n  ✅ Converged within 1 meV at: {ecut_values[converged_idx]:.0f} Ry")
            recommended_ecut = ecut_values[converged_idx]
        else:
            print(f"\n  ⚠️  Not converged within 1 meV at maximum tested cutoff")
            print(f"  Suggest testing higher cutoffs")
            recommended_ecut = ecut_values[-1]
        
        # Show convergence details
        print("\n  Convergence details:")
        for i, (ecut, diff) in enumerate(zip(ecut_values, energy_diff)):
            status = "✓" if diff < 1.0 else " "
            print(f"    {status} ecutwfc = {ecut:3.0f} Ry: {diff:8.3f} meV from converged value")
    
    print("\n" + "-"*50)
    print("✅ RECOMMENDED CUTOFFS:")
    print("-"*50)
    print(f"  ecutwfc = {recommended_ecut:.0f} Ry")
    print(f"  ecutrho = {4.0 * recommended_ecut:.0f} Ry (using 4x rule)")
    print("="*70)
    
    return recommended_ecut


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    
    print("="*70)
    print("WATER MOLECULE CONVERGENCE TEST")
    print("Testing convergence with respect to ecutwfc")
    print("Using STATIC (single point) calculations at fixed geometry")
    print("="*70)
    
    # Print initial geometry
    print("\n📐 Initial water geometry (fixed for all calculations):")
    print(f"  O-H1 bond length: {h2o.get_distance(0, 1):.4f} Å")
    print(f"  O-H2 bond length: {h2o.get_distance(0, 2):.4f} Å")
    print(f"  H-O-H angle: {h2o.get_angle(1, 0, 2):.2f}°")
    
    # =========================================================================
    # Test: ecutwfc convergence
    # =========================================================================
    print("\n🔬 Starting ecutwfc convergence test...")
    print("This will run static (single point) calculations for each cutoff value.")
    print("No geometry optimization will be performed.")
    
    # Define the range of cutoffs to test
    # For a comprehensive test, use more points:
    ecut_range = [20, 30, 40, 50, 60, 70, 80, 90, 100]
    
    # For a quicker test, use fewer points:
    # ecut_range = [30, 50, 70, 90]
    
    # For a very quick test:
    # ecut_range = [40, 60, 80]
    
    print(f"\nTesting ecutwfc values: {ecut_range} Ry")
    print(f"Total calculations: {len(ecut_range)}")
    print("\nStarting calculations...\n")
    
    # Run the convergence test with static calculations
    ecut_results = convergence_test_ecutwfc(h2o, ecut_range, ecutrho_factor=4.0)
    
    # =========================================================================
    # Plot and summarize results
    # =========================================================================
    print("\n📈 Generating plots and summary...")
    
    # Plot results
    plot_convergence_results(ecut_results)
    
    # Print summary with recommendations
    recommended_ecut = print_summary(ecut_results)
    
    print("\n✨ Convergence test completed successfully!")
    print(f"Recommended ecutwfc = {recommended_ecut} Ry")
    print(f"Recommended ecutrho = {4.0 * recommended_ecut:.1f} Ry")
