H2O molecule with QE
====================

cutoff convergence test

https://chat.deepseek.com/share/d5phgr1r1aoum7dlx5


(molmatmodel) milias@DESKTOP-7OTLCGO:~/work/projects/schools/Molecular-and-Materials-Modeling-2026-August/molecular-modeling/ase_mace/ase_tests/04_quantum_espresso_water/water_convergence_test/.python water_ecutwfc_convergence_test.py
======================================================================
WATER MOLECULE CONVERGENCE TEST
Testing convergence with respect to ecutwfc
Using STATIC (single point) calculations at fixed geometry
ecutrho determined automatically by QE from pseudopotential
======================================================================

📐 Initial water geometry (fixed for all calculations):
  O-H1 bond length: 0.9686 Å
  O-H2 bond length: 0.9686 Å
  H-O-H angle: 104.00°

🔬 Starting ecutwfc convergence test...
This will run static (single point) calculations for each cutoff value.
No geometry optimization will be performed.

Testing ecutwfc values: [20, 30, 40, 50, 60, 70, 80, 90, 100] Ry
Total calculations: 9

Starting calculations...


======================================================================
CONVERGENCE TEST: ecutwfc variation (Static calculations)
ecutrho determined automatically by QE from pseudopotential
Date: 2026-08-26 16:50:12
======================================================================
ecutwfc (Ry) |     Energy (eV) | ΔE (meV)
----------------------------------------------------------------------
        20.0 |     -596.747227 |       --
        30.0 |     -598.681798 | -1934.570 meV
        40.0 |     -598.767820 |  -86.023 meV
        50.0 |     -598.797486 |  -29.666 meV
        60.0 |     -598.808741 |  -11.255 meV
        70.0 |     -598.809060 |   -0.319 meV
        80.0 |     -598.810684 |   -1.624 meV
        90.0 |     -598.811290 |   -0.606 meV
       100.0 |     -598.811975 |   -0.685 meV
======================================================================

📈 Generating plots and summary...

Convergence plot saved as 'convergence_test_ecutwfc.png'

======================================================================
CONVERGENCE TEST SUMMARY
======================================================================

📊 ecutwfc Convergence (Static Calculations):
--------------------------------------------------
  Minimum tested: 20 Ry
  Maximum tested: 100 Ry
  Total energy range: -2064.75 meV

  ✅ Converged within 1 meV at: 90 Ry

  Convergence details:
      ecutwfc =  20 Ry: 2064.748 meV from converged value
      ecutwfc =  30 Ry:  130.177 meV from converged value
      ecutwfc =  40 Ry:   44.155 meV from converged value
      ecutwfc =  50 Ry:   14.489 meV from converged value
      ecutwfc =  60 Ry:    3.234 meV from converged value
      ecutwfc =  70 Ry:    2.914 meV from converged value
      ecutwfc =  80 Ry:    1.291 meV from converged value
    ✓ ecutwfc =  90 Ry:    0.685 meV from converged value
    ✓ ecutwfc = 100 Ry:    0.000 meV from converged value

--------------------------------------------------
✅ RECOMMENDED CUTOFFS:
--------------------------------------------------
  ecutwfc = 90 Ry
  ecutrho = NOT SET (QE determines automatically from pseudopotential)
======================================================================

✨ Convergence test completed successfully!
Recommended ecutwfc = 90 Ry
ecutrho will be determined automatically by QE


