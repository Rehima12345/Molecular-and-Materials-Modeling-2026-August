#!/usr/bin/env python3
"""
Challenge: Compute atomization energies of diatomics using MACE and EMT
Molecules: H2, N2, O2, F2, Cl2, etc.
Compares with experimental values
"""

import os
import sys
import numpy as np
from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from ase.io import read, write
import matplotlib.pyplot as plt

# ============================================================
# Import MACE (with error handling)
# ============================================================
try:
    from mace.calculators import mace_mp
    MACE_AVAILABLE = True
    print("✓ MACE imported successfully")
except ImportError:
    MACE_AVAILABLE = False
    print("⚠ MACE not available. Install with: pip install mace-torch")
    print("  Only EMT calculations will be performed")

# ============================================================
# Experimental reference data
# ============================================================
EXP_DATA = {
    'H2': {'bond_length': 0.741, 'atomization': 4.52},
    'N2': {'bond_length': 1.098, 'atomization': 9.76},
    'O2': {'bond_length': 1.208, 'atomization': 5.12},
    'F2': {'bond_length': 1.412, 'atomization': 1.59},
    'Cl2': {'bond_length': 1.988, 'atomization': 2.48},
    'Br2': {'bond_length': 2.281, 'atomization': 1.97},
    'I2': {'bond_length': 2.665, 'atomization': 1.54},
    'CO': {'bond_length': 1.128, 'atomization': 11.09},
    'NO': {'bond_length': 1.151, 'atomization': 6.50},
    'HCl': {'bond_length': 1.275, 'atomization': 4.43},
    'HF': {'bond_length': 0.917, 'atomization': 5.87},
}

# ============================================================
# Functions
# ============================================================

def get_diatomic_symbols(molecule_name):
    """
    Convert molecule name to list of atomic symbols
    """
    symbol_map = {
        'H2': ['H', 'H'],
        'N2': ['N', 'N'],
        'O2': ['O', 'O'],
        'F2': ['F', 'F'],
        'Cl2': ['Cl', 'Cl'],
        'Br2': ['Br', 'Br'],
        'I2': ['I', 'I'],
        'CO': ['C', 'O'],
        'NO': ['N', 'O'],
        'HCl': ['H', 'Cl'],
        'HF': ['H', 'F'],
    }
    return symbol_map.get(molecule_name, [molecule_name[0], molecule_name[0]])

def compute_atomization_energy(molecule_name, calculator, initial_distance=None, method_name="EMT"):
    """
    Compute atomization energy for a diatomic molecule using given calculator
    
    Args:
        molecule_name (str): Name of molecule (e.g., 'N2', 'O2')
        calculator: ASE calculator (EMT or MACE)
        initial_distance (float): Initial bond length (Å), uses experimental if None
        method_name (str): Name of method for printing
    
    Returns:
        dict: Results including bond length, energies, atomization energy
    """
    print(f"\n{'='*60}")
    print(f"{method_name} Calculation: {molecule_name}")
    print(f"{'='*60}")
    
    symbols = get_diatomic_symbols(molecule_name)
    
    # Get initial distance
    if initial_distance is None:
        initial_distance = EXP_DATA.get(molecule_name, {}).get('bond_length', 1.2)
    
    # Single atom energy
    atom = Atoms(symbols[0], calculator=calculator)
    e_atom = atom.get_potential_energy()
    print(f"  Single {symbols[0]} atom energy: {e_atom:.6f} eV")
    
    # Diatomic molecule optimization
    molecule = Atoms(symbols, [(0., 0., 0.), (0., 0., initial_distance)], calculator=calculator)
    
    opt = BFGS(molecule, trajectory=f"{molecule_name}_{method_name}.traj")
    print(f"  Optimizing {molecule_name} (initial d = {initial_distance:.3f} Å)")
    opt.run(fmax=0.01)
    
    # Results
    bond_length = molecule.get_distance(0, 1)
    e_molecule = molecule.get_potential_energy()
    e_atomization = 2 * e_atom - e_molecule
    
    print(f"\n  Results:")
    print(f"    Optimized bond length: {bond_length:.4f} Å")
    print(f"    Molecule energy:       {e_molecule:.6f} eV")
    print(f"    Atomization energy:    {e_atomization:.4f} eV")
    print(f"    Optimization steps:    {opt.nsteps}")
    
    return {
        'bond_length': bond_length,
        'e_molecule': e_molecule,
        'e_atom': e_atom,
        'e_atomization': e_atomization,
        'n_steps': opt.nsteps
    }

def compare_with_experiment(results, molecule_name):
    """
    Compare computed results with experimental values
    """
    exp = EXP_DATA.get(molecule_name)
    if exp is None:
        return
    
    print(f"\n  Comparison with Experiment:")
    
    # Bond length
    bond_diff = results['bond_length'] - exp['bond_length']
    bond_pct = abs(bond_diff) / exp['bond_length'] * 100
    print(f"    Bond length:  {results['bond_length']:.4f} Å  (exp: {exp['bond_length']:.4f} Å)")
    print(f"      Δ = {bond_diff:+.4f} Å  ({bond_pct:.2f}%)")
    
    # Atomization energy
    energy_diff = results['e_atomization'] - exp['atomization']
    energy_pct = abs(energy_diff) / exp['atomization'] * 100
    print(f"    Atomization:  {results['e_atomization']:.4f} eV  (exp: {exp['atomization']:.2f} eV)")
    print(f"      Δ = {energy_diff:+.4f} eV  ({energy_pct:.2f}%)")

def plot_comparison(emt_results, mace_results):
    """
    Create comparison plots
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Get molecules that have both EMT and MACE results
    molecules = [m for m in emt_results.keys() if m in mace_results]
    
    if not molecules:
        print("No molecules with both EMT and MACE results for plotting")
        return
    
    # Bond length comparison
    x = np.arange(len(molecules))
    width = 0.25
    
    exp_bonds = [EXP_DATA[m]['bond_length'] for m in molecules]
    emt_bonds = [emt_results[m]['bond_length'] for m in molecules]
    mace_bonds = [mace_results[m]['bond_length'] for m in molecules]
    
    ax1.bar(x - width, exp_bonds, width, label='Experiment', color='black', alpha=0.7)
    ax1.bar(x, emt_bonds, width, label='EMT', color='blue', alpha=0.7)
    ax1.bar(x + width, mace_bonds, width, label='MACE', color='green', alpha=0.7)
    
    ax1.set_xlabel('Molecule')
    ax1.set_ylabel('Bond Length (Å)')
    ax1.set_title('Bond Length Comparison')
    ax1.set_xticks(x)
    ax1.set_xticklabels(molecules)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Atomization energy comparison
    exp_energy = [EXP_DATA[m]['atomization'] for m in molecules]
    emt_energy = [emt_results[m]['e_atomization'] for m in molecules]
    mace_energy = [mace_results[m]['e_atomization'] for m in molecules]
    
    ax2.bar(x - width, exp_energy, width, label='Experiment', color='black', alpha=0.7)
    ax2.bar(x, emt_energy, width, label='EMT', color='blue', alpha=0.7)
    ax2.bar(x + width, mace_energy, width, label='MACE', color='green', alpha=0.7)
    
    ax2.set_xlabel('Molecule')
    ax2.set_ylabel('Atomization Energy (eV)')
    ax2.set_title('Atomization Energy Comparison')
    ax2.set_xticks(x)
    ax2.set_xticklabels(molecules)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('diatomic_atomization_comparison.png', dpi=150)
    print("\n✓ Plot saved to: diatomic_atomization_comparison.png")
    plt.show()

# ============================================================
# Main Execution
# ============================================================

def main():
    print("="*70)
    print("CHALLENGE: Diatomic Atomization Energies")
    print("Methods: EMT and MACE (Machine Learning Potential)")
    print("="*70)
    
    # List of molecules to test
    molecules = ['H2', 'N2', 'O2', 'F2']  # Add more as needed
    
    # Results storage
    emt_results = {}
    mace_results = {}
    
    # ============================================================
    # EMT Calculations
    # ============================================================
    print("\n" + "="*70)
    print("PART 1: EMT CALCULATIONS")
    print("="*70)
    
    emt = EMT()
    
    for mol in molecules:
        try:
            result = compute_atomization_energy(mol, emt, method_name="EMT")
            compare_with_experiment(result, mol)
            emt_results[mol] = result
        except Exception as e:
            print(f"  ❌ Error with {mol}: {e}")
            emt_results[mol] = None
    
    # ============================================================
    # MACE Calculations (if available)
    # ============================================================
    if MACE_AVAILABLE:
        print("\n" + "="*70)
        print("PART 2: MACE CALCULATIONS")
        print("="*70)
        
        # Use float64 for better accuracy
        try:
            mace_calc = mace_mp(model="medium", device="cpu", default_dtype="float64")
            print("✓ MACE calculator initialized (float64 precision)")
        except Exception as e:
            print(f"⚠ Error initializing MACE: {e}")
            print("  Trying with default settings...")
            mace_calc = mace_mp(model="medium", device="cpu")
        
        for mol in molecules:
            try:
                result = compute_atomization_energy(mol, mace_calc, method_name="MACE")
                compare_with_experiment(result, mol)
                mace_results[mol] = result
            except Exception as e:
                print(f"  ❌ Error with {mol}: {e}")
                mace_results[mol] = None
    
    # ============================================================
    # Summary Table
    # ============================================================
    print("\n" + "="*70)
    print("SUMMARY TABLE")
    print("="*70)
    
    print(f"\n{'Molecule':<10} {'Property':<15} {'EMT':<12} {'MACE':<12} {'Exp.':<12}")
    print("-"*65)
    
    for mol in molecules:
        if mol not in emt_results or emt_results[mol] is None:
            continue
        
        exp = EXP_DATA.get(mol, {})
        
        # Bond length
        emt_bond = emt_results[mol]['bond_length'] if emt_results[mol] else 0
        mace_bond = mace_results[mol]['bond_length'] if mol in mace_results and mace_results[mol] else 0
        exp_bond = exp.get('bond_length', 0)
        
        print(f"{mol:<10} {'Bond Length':<15} {emt_bond:<12.4f} {mace_bond:<12.4f} {exp_bond:<12.4f}")
        
        # Atomization energy
        emt_energy = emt_results[mol]['e_atomization'] if emt_results[mol] else 0
        mace_energy = mace_results[mol]['e_atomization'] if mol in mace_results and mace_results[mol] else 0
        exp_energy = exp.get('atomization', 0)
        
        print(f"{mol:<10} {'Atomization':<15} {emt_energy:<12.4f} {mace_energy:<12.4f} {exp_energy:<12.2f}")
    
    # ============================================================
    # Error Analysis
    # ============================================================
    print("\n" + "="*70)
    print("ERROR ANALYSIS")
    print("="*70)
    
    emt_bond_errors = []
    mace_bond_errors = []
    emt_energy_errors = []
    mace_energy_errors = []
    
    for mol in molecules:
        if mol not in emt_results or emt_results[mol] is None:
            continue
        
        exp = EXP_DATA.get(mol)
        if exp is None:
            continue
        
        # Bond length errors
        emt_bond = emt_results[mol]['bond_length']
        mace_bond = mace_results[mol]['bond_length'] if mol in mace_results and mace_results[mol] else None
        
        emt_bond_errors.append(abs(emt_bond - exp['bond_length']))
        if mace_bond:
            mace_bond_errors.append(abs(mace_bond - exp['bond_length']))
        
        # Energy errors
        emt_energy = emt_results[mol]['e_atomization']
        mace_energy = mace_results[mol]['e_atomization'] if mol in mace_results and mace_results[mol] else None
        
        emt_energy_errors.append(abs(emt_energy - exp['atomization']))
        if mace_energy:
            mace_energy_errors.append(abs(mace_energy - exp['atomization']))
    
    if emt_bond_errors:
        print(f"\nBond Length MAE:")
        print(f"  EMT:  {np.mean(emt_bond_errors):.4f} Å")
        if mace_bond_errors:
            print(f"  MACE: {np.mean(mace_bond_errors):.4f} Å")
    
    if emt_energy_errors:
        print(f"\nAtomization Energy MAE:")
        print(f"  EMT:  {np.mean(emt_energy_errors):.4f} eV")
        if mace_energy_errors:
            print(f"  MACE: {np.mean(mace_energy_errors):.4f} eV")
    
    # ============================================================
    # Plot
    # ============================================================
    if MACE_AVAILABLE:
        plot_comparison(emt_results, mace_results)
    
    # ============================================================
    # Save Results
    # ============================================================
    with open('diatomic_atomization_results.txt', 'w') as f:
        f.write("="*70 + "\n")
        f.write("DIATOMIC ATOMIZATION ENERGY RESULTS\n")
        f.write("="*70 + "\n\n")
        
        for mol in molecules:
            if mol not in emt_results or emt_results[mol] is None:
                continue
            
            f.write(f"\n{mol}:\n")
            f.write(f"  EMT Bond Length:    {emt_results[mol]['bond_length']:.4f} Å\n")
            f.write(f"  EMT Atomization:    {emt_results[mol]['e_atomization']:.4f} eV\n")
            
            if mol in mace_results and mace_results[mol]:
                f.write(f"  MACE Bond Length:   {mace_results[mol]['bond_length']:.4f} Å\n")
                f.write(f"  MACE Atomization:   {mace_results[mol]['e_atomization']:.4f} eV\n")
            
            exp = EXP_DATA.get(mol, {})
            f.write(f"  Exp. Bond Length:   {exp.get('bond_length', 0):.4f} Å\n")
            f.write(f"  Exp. Atomization:   {exp.get('atomization', 0):.2f} eV\n")
    
    print("\n✓ Results saved to: diatomic_atomization_results.txt")
    
    print("\n" + "="*70)
    print("✓ Challenge Complete!")
    print("="*70)

if __name__ == "__main__":
    main()
