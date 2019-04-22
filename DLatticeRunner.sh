#!/bin/bash

GAMMA="2.0	3.0	4.0	5.0	7.5	10.0	15.0	20.0	25.0"

for gamma in $GAMMA; do
	awk -v gamma=$gamma '{if($1 == "Gamma_in:" {print gamma" GammaP_in"} else {print $0}}' lattice_inputs.txt > lattice_inputs.txt
	mpirun -n 8 ../DLattice.exe
	
	filename='Latticeout_'$gamma'_gamma'
	mv ./out.dat ./$filename
done
