Hg@C18
======

(mace_env) miroi@MIRO:~/work/projects/mol_mat_modeling_schools/Molecular-and-Materials-Modeling-2026-August/materials-modeling/8_Adsorption/Hg_on_C18/.python relaxation_Hg_C18.py
Successfully read 19 atoms from c18_hg.vasp

Starting atom relaxation by ASE

Final results
  Total Energy                  : -7481.346431 eV
  ASE-style max force (norm)    : 0.009006 eV/Å
  QE-style max force            : 0.009006 eV/Å
  Pressure                      : 3.795311 kbar

Final relaxed structure saved to: final_relaxed_structure.vasp

Relaxation complete

C18 slab
========
(mace_env) miroi@MIRO:~/work/projects/mol_mat_modeling_schools/Molecular-and-Materials-Modeling-2026-August/materials-modeling/8_Adsorption/Hg_on_C18/.python energy_c18.py
Successfully read 18 atoms from c18.vasp

Running SCF calculation...
  Total energy: -2947.320598 eV
  Fermi level: -1.698800 eV

Hg atom
=======
python energy_Hg.py
Successfully read 1 atoms from hg.vasp

Running SCF calculation...
  Total energy: -4533.787985 eV
  Fermi level: -2.274400 eV


Adsorption energy = E(Hg@C18)-E(C18)-H(Hg)=-7481.346431-(-2947.320598)-(-4533.787985)=-.237848 = -0.24 eV


