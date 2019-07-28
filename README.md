# Lattice-Simulations
One dimensional lattice, outputs average temperature given intensity and gamma


This project is coded in Fortran 90, though there will be remnants from 77.
MPICH is a necessity for running this.
I'll have another sub-section somewhere for downloading and installing that.

Remember that before running Fortran code, it will need to actually be compiled on your machine.
I bring this up because some of my machines compile slightly differently, despite being the same?

Assuming you have MPICH on your computer, you can compile the data first by navigating to the directory that the Fortran file is in
Then, type:
    mpifort -o DLattice.exe DLattice.f
Once the code is compiled, you can run some of the shell scripts to make the whole process much less tedious.
In order, 'mpifort' calls mpi to compile, '-o' says output it, 'DLattice.exe' is the resulting executable, 'DLattice.f' is the Fotran you want to compile

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
It will save all of the outputs with a name saying which gamma was used
Put in the values of gamma that you want to run, run the script, and it'll spit out the data file with all of the values you want
This is the method I'd recommend unless you're testing for bugs in the Fortran.

UPDATE 06/04/19
So we've at this point spent quite a bit of time looking over the input variables for this whole thing
As a reminder, since we define hbar, m, and k as 1, that means that recoil frequency and energy are equal to 1/2
That means that values scaled by this (U0 and GammaPrime) are MULTIPLIED by 1/2
You might think that a small change in these values would only slightly change the behavior of the lattice
I seemed to also think that as I searched for bugs and explored
I was wrong. I was extremely wrong.
Besides that, the time spent IN lattice is also extremely important, more than I gave it credit for
Running 10,000 time steps and a time step of 0.01 would NOT show the rise of momentum from spontaneous emission
This would result in a pretty regular slope upwards, but not the rapid increase when U0 is very small
Then either changing to 100,000 steps or a larger step (0.1) will both fix this
My only fear for 0.1 isn't for DLattice, but for the daughter programs that might have fast-frequency terms
