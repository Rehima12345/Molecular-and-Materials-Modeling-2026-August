# Save as: ase_benchmarks_02_fixed.py
"""
Challenge I.2: Enhanced Performance Tests - FIXED QE Version
"""

import os
import sys
import time
import subprocess
import shutil
import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.calculators.calculator import Calculator
from ase.calculators.mopac import MOPAC
import matplotlib.pyplot as plt

# =====================================================================
# 1. Environment & Path Resolution
# =====================================================================
conda_prefix = os.environ.get('CONDA_PREFIX')
if conda_prefix:
    os.environ['ASE_MOPAC_COMMAND'] = f"{os.path.join(conda_prefix, 'bin', 'mopac')} PREFIX.mop"
    xtb_bin = os.path.join(conda_prefix, 'bin', 'xtb')
    nwchem_bin = os.path.join(conda_prefix, 'bin', 'nwchem')
else:
    xtb_bin = "xtb"
    nwchem_bin = "nwchem"

# =====================================================================
# 2. Custom Calculator Wrappers
# =====================================================================
class CustomPySCFCalculator(Calculator):
    implemented_properties = ['energy']
    def __init__(self, method='RHF', basis='3-21g', **kwargs):
        Calculator.__init__(self, **kwargs)
        self.method = method
        self.basis = basis

    def calculate(self, atoms=None, properties=['energy'], system_changes=['positions', 'numbers']):
        Calculator.calculate(self, atoms, properties, system_changes)
        import pyscf
        xyz_coords = [f"{sym} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}" 
                      for sym, pos in zip(self.atoms.get_chemical_symbols(), self.atoms.get_positions())]
        atom_str = "; ".join(xyz_coords)
        mol = pyscf.gto.Mole()
        mol.atom = atom_str
        mol.basis = self.basis
        mol.verbose = 0
        mol.build()
        mf = pyscf.scf.RHF(mol) if self.method.upper() == 'RHF' else pyscf.scf.KS(mol)
        if self.method.upper() != 'RHF': mf.xc = self.method
        self.results['energy'] = mf.kernel() * 27.211386245988

class CustomXTBCalculator(Calculator):
    implemented_properties = ['energy']
    def __init__(self, method='2', **kwargs):
        Calculator.__init__(self, **kwargs)
        self.method = method

    def calculate(self, atoms=None, properties=['energy'], system_changes=['positions', 'numbers']):
        Calculator.calculate(self, atoms, properties, system_changes)
        xyz_file = "tmp_bench_xtb.xyz"
        with open(xyz_file, "w") as f:
            f.write(f"{len(self.atoms)}\n\n")
            for sym, pos in zip(self.atoms.get_chemical_symbols(), self.atoms.get_positions()):
                f.write(f"{sym} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
        try:
            cmd = [xtb_bin, xyz_file, "--gfn", self.method]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            for line in res.stdout.splitlines():
                if "total energy" in line.lower():
                    for token in line.split():
                        try:
                            self.results['energy'] = float(token.strip('|').strip()) * 27.211386245988
                            break
                        except ValueError: continue
                    break
        finally:
            for f in [xyz_file, "xtbopt.xyz", "xtbopt.log", "wbo", "charges", "xtbrestart", "gfnff_topo"]:
                if os.path.exists(f): os.remove(f)

class CustomNWChemCalculator(Calculator):
    implemented_properties = ['energy']
    def __init__(self, method='dft', basis='3-21g', **kwargs):
        Calculator.__init__(self, **kwargs)
        self.method = method
        self.basis = basis

    def calculate(self, atoms=None, properties=['energy'], system_changes=['positions', 'numbers']):
        Calculator.calculate(self, atoms, properties, system_changes)
        nw_file = "tmp_bench_nwchem.nw"
        
        with open(nw_file, "w") as f:
            f.write("geometry\n")
            for sym, pos in zip(self.atoms.get_chemical_symbols(), self.atoms.get_positions()):
                f.write(f"  {sym} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
            f.write(f"end\nbasis\n  * library {self.basis}\nend\ntask {self.method} energy\n")
            
        try:
            res = subprocess.run([nwchem_bin, nw_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            energy = None
            for line in res.stdout.splitlines():
                if "total dft energy" in line.lower() or "total energy" in line.lower():
                    for token in line.split():
                        try:
                            energy = float(token) * 27.211386245988
                            break
                        except ValueError: continue
                    if energy is not None: break
            
            if energy is not None:
                self.results['energy'] = energy
            else:
                raise RuntimeError(f"Could not parse energy header. NWChem exit status: {res.returncode}")
        finally:
            if os.path.exists(nw_file): os.remove(nw_file)

# =====================================================================
# 3. QE Calculator with proper setup
# =====================================================================
def get_qe_calculator(atoms=None):
    """Create QE calculator with proper configuration"""
    try:
        from ase.calculators.espresso import Espresso, EspressoProfile
        
        # Find pw.x
        pw_path = shutil.which('pw.x')
        if pw_path is None:
            return None
        
        # Use current directory for pseudopotentials
        pseudo_dir = os.getcwd()
        
        # Check if pseudopotentials exist
        pseudo_files = ['H.upf', 'O.upf']
        missing = [f for f in pseudo_files if not os.path.exists(os.path.join(pseudo_dir, f))]
        if missing:
            print(f"  Warning: Missing pseudopotentials: {missing}")
            # Create dummy pseudopotential if needed for testing
            for f in missing:
                with open(os.path.join(pseudo_dir, f), 'w') as pf:
                    pf.write("Dummy pseudopotential for testing\n")
        
        # Create profile with required arguments
        profile = EspressoProfile(
            command=pw_path,
            pseudo_dir=pseudo_dir
        )
        
        # Input parameters for a quick test
        input_data = {
            'control': {
                'calculation': 'scf',
                'restart_mode': 'from_scratch',
                'prefix': 'qe_test',
                'tprnfor': True,
                'tstress': True,
                'outdir': './qe_tmp/',
            },
            'system': {
                'ecutwfc': 30.0,
                'ecutrho': 120.0,
                'occupations': 'smearing',
                'smearing': 'gaussian',
                'degauss': 0.02,
            },
            'electrons': {
                'diagonalization': 'david',
                'conv_thr': 1e-6,
                'mixing_beta': 0.7,
            },
        }
        
        # Create calculator
        calc = Espresso(
            profile=profile,
            pseudopotentials={'H': 'H.upf', 'O': 'O.upf'},
            input_data=input_data,
            kpts=(1, 1, 1),
        )
        
        return calc
        
    except Exception as e:
        print(f"  Error setting up QE: {e}")
        return None

# =====================================================================
# 4. Enhanced Cluster Generation
# =====================================================================
def generate_cluster(molecule_type='H2O', n_molecules=1, spacing=3.5):
    """Generate clusters of different molecules"""
    cluster = Atoms()
    for i in range(n_molecules):
        mol = molecule(molecule_type)
        mol.translate([i * spacing, 0.0, 0.0])
        cluster += mol
    return cluster

# =====================================================================
# 5. Performance Analysis Functions
# =====================================================================
def analyze_performance(results, cluster_sizes, method_names):
    """Analyze and create performance plots"""
    
    # Extract timing data
    times = {}
    for method in method_names:
        times[method] = []
        for result in results[method]:
            if result != "FAILED":
                times[method].append(float(result.replace('s', '')))
            else:
                times[method].append(0)
    
    # Create plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Scaling with cluster size
    active_methods = [m for m in method_names if any(times[m])]
    
    for method in active_methods:
        if any(times[method]):
            ax1.plot(cluster_sizes, times[method], 'o-', label=method, linewidth=2, markersize=8)
    
    ax1.set_xlabel('Number of Molecules', fontsize=12)
    ax1.set_ylabel('Wall Time (seconds)', fontsize=12)
    ax1.set_title('Performance Scaling with System Size', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Bar chart for largest cluster
    largest_idx = len(cluster_sizes) - 1
    methods_with_time = []
    times_largest = []
    
    for method in method_names:
        if times[method][largest_idx] > 0:
            methods_with_time.append(method)
            times_largest.append(times[method][largest_idx])
    
    if methods_with_time:
        bars = ax2.bar(range(len(methods_with_time)), times_largest, 
                      color=['skyblue', 'lightcoral', 'lightgreen', 'gold', 'plum'],
                      edgecolor='black')
        ax2.set_xticks(range(len(methods_with_time)))
        ax2.set_xticklabels(methods_with_time, rotation=45, ha='right')
        ax2.set_ylabel('Wall Time (seconds)', fontsize=12)
        ax2.set_title(f'Performance for {cluster_sizes[largest_idx]} Molecules', fontsize=14)
        ax2.grid(True, alpha=0.3)
        
        # Add value labels
        for bar, time in zip(bars, times_largest):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                    f'{time:.3f}s', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig('performance_analysis.png', dpi=150)
    print("\n✓ Performance plot saved to: performance_analysis.png")
    plt.close()

# =====================================================================
# 6. Main Benchmark Setup - Enhanced
# =====================================================================
if __name__ == "__main__":
    print("="*60)
    print("CHALLENGE I.2: ENHANCED PERFORMANCE TESTS")
    print("Running NWChem, Quantum ESPRESSO, and MOPAC performance tests")
    print("="*60)
    
    # Test different cluster sizes
    cluster_sizes = [1, 2, 4, 8]
    
    # Methods to test
    methods = {
        'MOPAC (PM7)': lambda: MOPAC(method='PM7'),
        'xTB (GFN2)': lambda: CustomXTBCalculator(method='2'),
        'PySCF (RHF/3-21G)': lambda: CustomPySCFCalculator(method='RHF', basis='3-21g'),
        'NWChem (DFT/3-21G)': lambda: CustomNWChemCalculator(method='dft', basis='3-21g'),
    }
    
    # Try to add Quantum ESPRESSO
    try:
        from ase.calculators.espresso import Espresso, EspressoProfile
        if shutil.which('pw.x'):
            methods['QE (DFT)'] = lambda: get_qe_calculator()
            print("✓ Quantum ESPRESSO added to benchmark")
        else:
            print("⚠ Quantum ESPRESSO not found - skipping")
    except Exception as e:
        print(f"⚠ Quantum ESPRESSO not available: {e}")
    
    results = {m: [] for m in methods}
    energy_results = {m: [] for m in methods}
    
    print("\n" + "="*60)
    print("RUNNING BENCHMARKS")
    print("="*60)
    
    for n in cluster_sizes:
        atoms = generate_cluster('H2O', n)
        n_atoms = len(atoms)
        print(f"\n🚀 Cluster Size: {n} H2O molecules ({n_atoms} atoms)")
        print("-" * 58)
        
        for name, calc_init in methods.items():
            try:
                atoms.calc = calc_init()
                
                start_time = time.perf_counter()
                energy = atoms.get_potential_energy()
                end_time = time.perf_counter()
                
                wall_time = end_time - start_time
                results[name].append(f"{wall_time:.4f}s")
                energy_results[name].append(f"{energy:.2f}")
                print(f"  |-- {name:<20} : {wall_time:>8.4f} seconds (E = {energy:.2f} eV)")
            except Exception as e:
                results[name].append("FAILED")
                energy_results[name].append("FAILED")
                print(f"  |-- {name:<20} : ❌ FAILED (Reason: {e})")
            
            # Clean up
            for junk in ['mopac.out', 'mopac.arc', 'mopac.mop', 'tmp_bench_nwchem.nw']:
                if os.path.exists(junk): os.remove(junk)
    
    # =====================================================================
    # 7. Enhanced Summary
    # =====================================================================
    print("\n\n" + "="*60)
    print("                    BENCHMARK SUMMARY                     ")
    print("="*60)
    
    # Timing table
    print("\nTiming Results (seconds):")
    header_str = f"| {'Method':<20} "
    for n in cluster_sizes:
        header_str += f"| {n} H2O ({n*3} atoms) "
    print(header_str + "|")
    print("|" + "----------------------|" * (len(cluster_sizes) + 1))
    
    for name in methods:
        row_str = f"| {name:<20} "
        for idx, n in enumerate(cluster_sizes):
            row_str += f"| {results[name][idx]:^15} "
        print(row_str + "|")
    
    # Energy table
    print("\nEnergy Results (eV):")
    header_str = f"| {'Method':<20} "
    for n in cluster_sizes:
        header_str += f"| {n} H2O ({n*3} atoms) "
    print(header_str + "|")
    print("|" + "----------------------|" * (len(cluster_sizes) + 1))
    
    for name in methods:
        row_str = f"| {name:<20} "
        for idx, n in enumerate(cluster_sizes):
            row_str += f"| {energy_results[name][idx]:^15} "
        print(row_str + "|")
    
    # Speed comparison for largest cluster
    print("\n" + "="*60)
    print("SPEED COMPARISON (Largest Cluster)")
    print("="*60)
    
    largest_idx = len(cluster_sizes) - 1
    times_list = []
    method_list = []
    
    for name in methods:
        if results[name][largest_idx] != "FAILED":
            time_val = float(results[name][largest_idx].replace('s', ''))
            times_list.append(time_val)
            method_list.append(name)
    
    if times_list:
        fastest = min(times_list)
        for name, time_val in zip(method_list, times_list):
            speedup = time_val / fastest
            print(f"  {name:<20}: {time_val:.4f}s (Speedup: {speedup:.2f}x)")
    
    # Analyze and plot
    print("\nGenerating performance analysis plots...")
    analyze_performance(results, cluster_sizes, list(methods.keys()))
    
    # Save results to file
    with open('benchmark_results.txt', 'w') as f:
        f.write("="*60 + "\n")
        f.write("CHALLENGE I.2: PERFORMANCE BENCHMARK RESULTS\n")
        f.write("="*60 + "\n\n")
        f.write(f"Date: {time.ctime()}\n")
        f.write(f"System: {os.uname().nodename}\n\n")
        
        f.write("Timing Results (seconds):\n")
        f.write(header_str + "|\n")
        f.write("|" + "----------------------|" * (len(cluster_sizes) + 1) + "\n")
        for name in methods:
            row_str = f"| {name:<20} "
            for idx in range(len(cluster_sizes)):
                row_str += f"| {results[name][idx]:^15} "
            f.write(row_str + "|\n")
        
        f.write("\nEnergy Results (eV):\n")
        f.write(header_str + "|\n")
        f.write("|" + "----------------------|" * (len(cluster_sizes) + 1) + "\n")
        for name in methods:
            row_str = f"| {name:<20} "
            for idx in range(len(cluster_sizes)):
                row_str += f"| {energy_results[name][idx]:^15} "
            f.write(row_str + "|\n")
    
    print("\n✓ Results saved to: benchmark_results.txt")
    print("="*60)
    print("CHALLENGE I.2 COMPLETE!")
    print("="*60)
