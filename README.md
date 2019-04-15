# Lattice-Simulations
One dimensional lattice, outputs average temperature given intensity and gamma


This project is coded in Fortran 90, though there will be remnants from 77.
MPICH is a necessity for running this.
I'll have another sub-section somewhere for downloading and installing that.

Remember that before running Fortran code, it will need to actually be compiled on your machine.
I bring this up because some of my machines compile slightly differently, despite being the same?

Once the code is compiled, you can run some of the shell scripts to make the whole process much less tedious.
Code should be compiled as "DLattice.exe" so the script can find it

FOR A SINGLE ITERATION
Figure out first what you want to run. Set whichever potentials you'd like to use in "lattice_potentials.txt" FIRST.
THEN go "lattice_inputs.txt" and adjust to the number of potentials you're using.
While you're there, set the other inputs as you'd like to run them.
Assuming you've already compiled, go into the command terminal, navigate to the folder that all the files are in
Type in:
    mpirun -n 8 ../DLattice.exe
And press enter. That "8" is telling MPI how many parallel processes to run. 
I recommend to set it to however many cores you have, since you're getting more data for the same amount of time as 1
All data will then be outputted to a file called "out.dat"
The file will contain three columns: final temperature, lattice potential, gamma prime
The Fortran code is designed to append this data file, not rewrite over it
That means you can run the code multiple times and it should keep track of your previous runs
THAT BEING SAID, it doesn't consider any changes in other variables.
If you change the time steps, or atom number, or dt value, it will continue appending that file
I recommend renaming the data file whenever you're done testing those variables

FOR USING THE SCRIPT
It's really very similar to a single iteration but it runs it all for you.
The first thing you'll need to do after downloading is change the permissions of the script so that it will actually run
In the terminal, navigate to the "DLatticeRunner.sh" file
Then type in:
    chmod +x ./DLatticeRunner.sh
This will give the script permission to run from now on, you won't need to do this as long as you keep using that particular script
Inside the script, you'll notice a variable named "GAMMA" with a number of different values
These are all of the Gamma prime values that will be fed into the lattice inputs
The script will alter the "lattice_inputs.txt" file, then execute the Fortran code for each value of GAMMA
There is NO EXTRA EFFORT REQUIRED.
Put in the values of gamma that you want to run, run the script, and it'll spit out the data file with all of the values you want
This is the method I'd recommend unless you're testing for bugs in the Fortran.
