#!/bin/bash
echo "=========================================="
echo "ecutwfc Convergence Test for Si"
echo "=========================================="
echo "ecutwfc (Ry) | Total Energy (Ry)"
echo "------------------------------------------"

for ecut in 20 30 40 50 60 80 100; do
    echo -n "Testing $ecut Ry ... "
    
    # Method: Copy Si_base.in and change ecutwfc using sed
    cp Si_base.in Si_${ecut}.in
    sed -i "s/ecutwfc = 30/ecutwfc = $ecut/" Si_${ecut}.in
    
    # Verify ecutwfc was changed
    echo "  ecutwfc in file: $(grep ecutwfc Si_${ecut}.in)"
    
    # Run calculation
    mpirun -np 2 pw.x < Si_${ecut}.in > Si_${ecut}.out 2>/dev/null
    
    # Get energy
    energy=$(grep "!" Si_${ecut}.out | awk '{print $5}')
    
    if [ -n "$energy" ]; then
        echo "  Total energy: $energy Ry"
    else
        echo "  FAILED to converge"
    fi
done
echo "=========================================="
