	program DLattice

! Beginning a new version of the DLattice code
! This is using my naming scheme and conventions from Fortran 90
! MPI is included in this, but hardly necessary
! The code runs very quickly as it doesn't have that many iterations to go through

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start by figuring out the variables we're going to be using
! This is very different from some of the other sample code we're using
! I prefer this way even if it's really tedious, using implicit none makes sure there's no rogue variables
	implicit none
! Variables for the simulation loops
	integer :: i,j,U_num,u ! for the loops
! Remember that atom_num is multiplied by the number of processors in MPI
	integer :: nstate,atom_num,tsteps,idum
	real*8 :: gammaP_in,dt_in,U_in 
	real*8 :: state,gammaP,U0,Dfsn,kick,deltaP,jumprate
	real*8 :: p,z,t,dt
	real*8 :: ninth,cos2z,rand,gasdev
	real*8 :: p2sum,restmp,tmpavg
	real*8, allocatable :: U_Er(:),finaltemp(:)

! Variables for MPI
	integer :: ierr,myid,numprocs

! Variables for the random seed generator
	integer, allocatable :: seed(:)
	integer :: n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we need to establish the MPI stuff
! This isn't entirely necessary for just the DLattice, but we'll be
! using this later as the base for the DBiharmonic and DMultiFreq.

	include 'mpif.h'

! Set up the variables we'll need to for the MPI stuff

! Now start the MPI function, give the ID of each process
	call MPI_INIT (ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we'll read from the input files for the program

! All input files are in 10-series, out in 20-series
! Are your values given in a column rather than row?
! Then your read functions need to also be in a column
	open(11,file='lattice_inputs.txt')
	read(11,*) atom_num
	read(11,*) tsteps
	read(11,*) U_num
	read(11,*) GammaP_in
	read(11,*) dt_in
! This is essentially asking how many potentials should we be looking for
! That way we can use some arbitrary number of potentials
	allocate(U_Er(U_num))
	allocate(finaltemp(U_num))
	open(12,file='lattice_potentials.txt')	
	read(12,*) U_Er
	


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we need a random number generator
! If we just feed the processors the same seed, they'll produce the
! same results for the simulation. This is bad.
! So intead, we are going to create seeds for the number generator that
! is based on a combination of time and processor ID. That way, the
! simulation I run tomorrow is different from today, and each process
! will be running its own simulation, rather than copying.

! Special note here!
! When you ask for the random_seed generator like this, it calls from the OS
! That means it's a little dependent on how good that RNG actually is
! Here, I'm trusting the OS to give me good random numbers
! If this is in doubt, it needs to be replaced by something else in a similar manner

	call random_seed(size=n)
	allocate(seed(n))
	call random_seed(get=seed)

	call random_seed(put=seed)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We have all of the variables, now we run the loops

! Do some quick setting up
	! Initial state for atoms
	nstate = 1
	! And state for using in calculations
	state = dble(nstate)

! Division takes a long time, do it once here
	ninth = 1.d0/9.d0
!!! Apparently my logic is wrong here

! We are talking hbar, k, and m to be equal to 1
! We scale some variables according to those values
! E_r = hbar^2 k^2 / 2m & E_r = hbar * w_r ==> E_r = 1/2 and w_r = 1/2
! ~U = U0/E_r ==> ~U = 2*U0
! So these values we scale to are effectively halved in the code
! Doing potentials one at a time in the loop

	gammaP = gammaP_in * 0.5d0
	
! We also scale the time step with gamma
! Typical value for dt_in is 0.01
	dt = dt_in/gammaP
! Are three for loops bad form?
! Is running scripts for each of these better?
! I don't really care
	do u = 1,U_num
! Remember that ACTUAL potential is scaled by Er
! That means U0/Er is U0 / (1/2)
	   U0 = U_Er(u) * 0.5d0
! And we also have to set our restmp to 0
	   restmp = 0.d0
	   tmpavg = 0.d0
! The second loop will scan through our atoms
	   do i = 1,atom_num
	   
! Set the initial conditions for the new atom
! These don't really matter a ton,
	      p = gasdev(idum)
	      z = gasdev(idum)
	      t = 0.d0

! And now run through the time steps
	      do j = 1,tsteps

	         cos2z = dcos(2.d0*z)
! First take into account the diffusion of the atom
! This is analytical value assuming ONLY diffusion of same-well
! Well jump is unlikely in grand scheme of things, also scales inverse with dt
! Results in sudden VERY large jumps in momentum with very small dt
	         Dfsn=(0.1d0*ninth*gammaP)*(35.d0+(state*7.d0*cos2z))
! Then how much momentum does the atom feel from this
! Random direction times the momentum to simulate random direction
	         kick = ((2.d0*Dfsn*dt)**(0.5d0))*gasdev(idum)
! Total change in momentum is kick and the force from the potential
	         deltaP = kick + (state*U0*dsin(2.d0*z)*dt)

! Now we have to ask if we changed wells or not

! Probability we jumped in this moment
	         jumprate = ninth*dt*gammaP*(1.d0+(state*cos2z))

	         call random_number(rand)

! Compare and see if we jumped
	         if (rand .lt. jumprate) then
! This is the path if we jumped
	            state = -state
	         end if
! And update the rest of the variables
	         p = p + deltaP
	         z = z + p*dt
	         t = t + dt

! For the lattice, we're interested in taking the average of momentum
! squared. Easiest way is to take a running sum
! Scale that sum by the total number we'll be summing
		 tmpavg = tmpavg + (p*p)/(atom_num*tsteps)

	      end do ! end of time steps
! Still waiting for the tmpavg to finish...
	   end do ! end of atom loop
! Now our tmpavg for THIS potential is done
! While we have it, let's reduce it and save
	   call MPI_REDUCE(tmpavg,restmp,1,MPI_DOUBLE_PRECISION,
     & MPI_SUM,0,MPI_COMM_WORLD,ierr)
! Now we should only do this on one core
! Otherwise it would just do this 8 times
	   if (myid.eq.0) then
! Now we have the all summed together restmp
! We need to divide by number of processes to get true average
! Then converting to kb T we need to multiply by 2
	      restmp = 2.d0 * restmp/numprocs
	      open(unit=21,file='out.dat',position='append')
	      write(21,*) restmp, U_Er(u), gammaP_in
	   end if

	end do ! end of potentials loop


! Now since we just did all of the saving in the loop, no need for end section

	call MPI_FINALIZE(ierr)
	end program

	FUNCTION GASDEV(IDUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      REAL*8 GASDEV,X(2)
      SAVE ISET,GSET
      DATA ISET/0/
      IF(ISET.EQ.0) THEN
 1       CALL RANDOM_NUMBER(X)
         V1=2.0d0*X(1)-1.d0
         V2=2.0d0*X(2)-1.d0
         R=V1**2+V2**2
      IF(R.GE.1.0d0) GO TO 1
         FAC=DSQRT(-2.0d0*DLOG(R)/R)
         GSET=V1*FAC
         GASDEV=V2*FAC
         ISET=1
      ELSE
         GASDEV=GSET
         ISET=0
      ENDIF
      RETURN
      END
