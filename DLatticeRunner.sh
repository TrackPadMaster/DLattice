#!/bin/bash

GAMMA="2.0 5.0 10.0 15.0"

for gamma in $GAMMA; do
	
	awk -v gamma=$gamma '{if($3 == "GammaPin") {print gamma" ! GammaPin"} else {print $0}}' lattice_inputs0.txt > lattice_inputs.txt
	mpirun -n 4 ./DLattice.exe
	
	filename='Latticeout_'$gamma'_gamma.txt'
	mv ./out.dat ./$filename
done
