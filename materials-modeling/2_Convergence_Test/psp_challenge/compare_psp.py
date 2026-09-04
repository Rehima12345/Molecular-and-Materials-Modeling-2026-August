#!/usr/bin/env python
"""
CHALLENGE: Compare convergence for different Si pseudopotentials
"""
from ase import Atoms
from ase.calculators.espresso import Espresso, EspressoProfile
import os
import sys

# Find all UPF files in current directory
psp_files = [f for f in os.listdir('.') if f.endswith('.UPF')]

if not psp_files:
    print("No pseudopotential files found!")
    sys.exit(1)

print("=" * 85)
print("CHALLENGE: PSEUDOPOTENTIAL CONVERGENCE COMPARISON FOR Si")
print("=" * 85)

# Atomic structure
atoms = Atoms(
    symbols=['Si']*2,
    positions=[
        [1.3574500000, 4.0723500000, 4.0723500000],
        [0.0000000000, 0.0000000000, 0.0000000000]
    ],
    cell=[
        [0.0000000000, 2.7149000000, 2.7149000000],
        [2.7149000000, 0.0000000000, 2.7149000000],
        [2.7149000000, 2.7149000000, 0.0000000000]
    ],
    pbc=[True, True, True]
)

# QE settings
qe_bin = "/home/rehim/anaconda3/envs/molmatmodel/bin"
pw_command = f'mpirun -np 2 {qe_bin}/pw.x'
pw_profile = EspressoProfile(command=pw_command, pseudo_dir='./')

# Test parameters
ecutwfc_values = [20, 30, 40, 50, 60, 70, 80]
k_test = (8, 8, 8)

# Store results for table
results = []

print("\n" + "-" * 85)
print(f"{'Pseudopotential':<45} {'Functional':<12} {'Type':<10} {'Converged ecutwfc (Ry)':<20}")
print("-" * 85)

for psp in sorted(psp_files):
    # Extract info from filename
    psp_name = os.path.basename(psp)
    
    # Determine functional
    if 'pbe' in psp_name and 'pbesol' not in psp_name:
        functional = 'PBE'
    elif 'pbesol' in psp_name:
        functional = 'PBESOL'
    elif 'pz' in psp_name:
        functional = 'LDA'
    else:
        functional = 'Unknown'
    
    # Determine pseudopotential type
    if 'kjpaw' in psp_name:
        pp_type = 'PAW'
    elif 'rrkjus' in psp_name:
        pp_type = 'USPP'
    else:
        pp_type = 'NC'
    
    print(f"\nTesting: {psp_name}")
    print(f"  Functional: {functional}, Type: {pp_type}")
    print("  " + "-" * 50)
    
    base_input = {
        'control': {'calculation': 'scf', 'prefix': 'si', 'outdir': './tmp', 'verbosity': 'low'},
        'system': {
            'ibrav': 0, 'nat': 2, 'ntyp': 1,
            'occupations': 'smearing', 'smearing': 'gaussian', 'degauss': 0.01
        },
        'electrons': {'conv_thr': 1.0e-8, 'mixing_beta': 0.7}
    }
    
    prev_energy = None
    converged_ecut = None
    
    for ecut in ecutwfc_values:
        input_data = base_input.copy()
        input_data['system']['ecutwfc'] = int(ecut)
        input_data['system']['ecutrho'] = int(ecut * 10)
        
        calc = Espresso(
            profile=pw_profile,
            pseudopotentials={'Si': psp},
            input_data=input_data,
            kpts=k_test
        )
        
        atoms.calc = calc
        try:
            energy = atoms.get_potential_energy()
            if prev_energy is not None:
                diff = abs(energy - prev_energy) / 2 * 1000  # meV/atom
                status = "✓" if diff < 1.0 else " "
                print(f"    ecutwfc = {ecut:2d} Ry: {energy:.6f} eV (ΔE = {diff:.3f} meV/atom) {status}")
                if diff < 1.0 and converged_ecut is None:
                    converged_ecut = ecut
            else:
                print(f"    ecutwfc = {ecut:2d} Ry: {energy:.6f} eV")
            prev_energy = energy
        except Exception as e:
            print(f"    ecutwfc = {ecut:2d} Ry: FAILED")
            break
    
    if converged_ecut is None:
        converged_ecut = '>80'
    
    results.append([psp_name, functional, pp_type, str(converged_ecut)])
    print(f"  >>> Converged ecutwfc: {converged_ecut} Ry")

print("\n" + "=" * 85)
print("SUMMARY TABLE")
print("=" * 85)
print(f"{'Pseudopotential':<45} {'Functional':<12} {'Type':<10} {'Converged ecutwfc (Ry)':<20}")
print("-" * 85)
for r in results:
    print(f"{r[0]:<45} {r[1]:<12} {r[2]:<10} {r[3]:<20}")
print("=" * 85)

# Write table to file
with open('psp_comparison_table.txt', 'w') as f:
    f.write("=" * 85 + "\n")
    f.write("CHALLENGE: PSEUDOPOTENTIAL CONVERGENCE COMPARISON FOR Si\n")
    f.write("=" * 85 + "\n")
    f.write(f"{'Pseudopotential':<45} {'Functional':<12} {'Type':<10} {'Converged ecutwfc (Ry)':<20}\n")
    f.write("-" * 85 + "\n")
    for r in results:
        f.write(f"{r[0]:<45} {r[1]:<12} {r[2]:<10} {r[3]:<20}\n")
    f.write("=" * 85 + "\n")

print("\n✓ Results saved to: psp_comparison_table.txt")
