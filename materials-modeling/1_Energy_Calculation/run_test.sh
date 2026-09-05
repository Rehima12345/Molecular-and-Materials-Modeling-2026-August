#!/bin/bash
echo "=========================================="
echo "ecutwfc Convergence Test for Si"
echo "=========================================="
echo "ecutwfc (Ry) | Total Energy (Ry)"
echo "------------------------------------------"

for ecut in 20 30 40 50 60 80 100; do
    echo -n "Testing $ecut Ry ... "
    
    # Create input with this ecutwfc
    sed "s/ecutwfc = .*/ecutwfc = $ecut/" Si_base.in > Si_${ecut}.in
    
    # Run calculation
    mpirun -np 2 pw.x < Si_${ecut}.in > Si_${ecut}.out 2>/dev/null
    
    # Get energy
    energy=$(grep "!" Si_${ecut}.out | awk '{print $5}')
    
    if [ -n "$energy" ]; then
        echo "$energy Ry"
    else
        echo "FAILED"
    fi
done
echo "=========================================="
