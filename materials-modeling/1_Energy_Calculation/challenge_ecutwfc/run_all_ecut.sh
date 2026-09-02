#!/bin/bash
echo "=========================================="
echo "ecutwfc Convergence Test for Si"
echo "=========================================="
echo "ecutwfc (Ry) | Total Energy (Ry)"
echo "------------------------------------------"

for ecut in 20 30 40 50 60 80 100; do
    echo -n "Testing $ecut Ry ... "
    
    # Create input file with this ecutwfc
    cat > Si_${ecut}.in << EOT
&CONTROL
   calculation      = 'scf'
   prefix           = 'si'
   pseudo_dir       = '../'
   outdir           = './tmp'
   verbosity        = 'low'
/
&SYSTEM
   ibrav            = 2
   celldm(1)        = 10.26
   nat              = 2
   ntyp             = 1
   ecutwfc          = $ecut
   ecutrho          = $((ecut * 10))
   occupations      = 'smearing'
   smearing         = 'gaussian'
   degauss          = 0.01
/
&ELECTRONS
   conv_thr         = 1e-8
   mixing_beta      = 0.7
/
ATOMIC_SPECIES
Si 28.0855 Si.upf
ATOMIC_POSITIONS crystal
Si 0.00 0.00 0.00
Si 0.25 0.25 0.25
K_POINTS automatic
8 8 8 0 0 0
EOT
    
    # Run calculation
    mpirun -np 2 pw.x < Si_${ecut}.in > Si_${ecut}.out 2>/dev/null
    
    # Check convergence
    if grep -q "convergence has been achieved" Si_${ecut}.out; then
        energy=$(grep "!" Si_${ecut}.out | awk '{print $5}')
        echo "✓ $energy Ry"
    else
        echo "✗ FAILED to converge"
    fi
done
echo "=========================================="
