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
	integer :: i,j ! for the loops
! Remember that atom_num is multiplied by the number of processors in MPI
	integer :: nstate,atom_num,tsteps,idum
	real*8 :: U_in,gammaP_in,dt_in 
	real*8 :: state,gammaP,U0,Dfsn,kick,deltaP,jumprate
	real*8 :: p,z,t,dt
	real*8 :: ninth,cos2z,rand,gasdev
	real*8 :: p2sum,restmp,tmpavg

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
	open(11,file='lattice_inputs.txt')
	read(11,*) atom_num,tsteps,U_in,GammaP_in,dt_in
	allocate(U0(U_num))
	open(12,file='lattice_potentials.txt')	
	read(12,*) U0
	allocate(finaltemp(U_num))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we need a random number generator
! If we just feed the processors the same seed, they'll produce the
! same results for the simulation. This is bad.
! So intead, we are going to create seeds for the number generator that
! is based on a combination of time and processor ID. That way, the
! simulation I run tomorrow is different from today, and each process
! will be running its own simulation, rather than copying.

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
! We are talking hbar, k, and m to be equal to 1
! We scale some variables according to those values
! E_r = hbar^2 k^2 / 2m & E_r = hbar * w_r ==> E_r = 1/2 and w_r = 1/2
! ~U = U0/E_r ==> ~U = 2*U0
! So these values we scale to are effectively halved in the code

	gammaP = gammaP_in / 2.d0
	U0 = U_in/2.d0
! We also scale the time step with gamma
	dt = dt_in/gammaP


! The first loop will scan through our atoms
	do i = 1,atom_num
	   
! Set the initial conditions for the new atom
	   p = gasdev(idum)
	   z = gasdev(idum)
	   t = 0.d0
! And now run through the time steps
	   do j = 1,tsteps

	      cos2z = dcos(2.d0*z)
! First take into account the diffusion of the atom
	      Dfsn=(0.1d0*ninth*gammaP)*(35.d0+(state*7*cos2z))
! Then how much momentum does the atom feel from this
! Random direction times the momentum to simulate random direction
	      kick = gasdev(idum)*(2.d0*Dfsn*dt)**(0.5)
! Total change in momentum is kick and the force from the potential
	      deltaP = kick + (state*U0*dsin(2.d0*z)*dt)
! Now we have to ask if we changed wells or not

! Probability we jumped in this moment
	      jumprate = ninth*dt*gammaP*(1+(state*cos2z))
	      call random_number(rand)
! Compare and see if we jumped
	      if (rand .lt. jumprate) then
	         nstate = -nstate
	         state = dble(nstate)
	      end if
! And update the rest of the variables
	      p = p + deltaP
	      z = z + p*dt
	      t = t + dt

! For the lattice, we're interested in taking the average of momentum
! squared. Easiest way is to take a running sum
	      p2sum = p*p
	   end do

	end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we take the pieces from the last section and analyze them

! Now we have this massive number p2sum that we can use
! We'll take that and divide it by the number of atoms and total steps
	tmpavg = p2sum/(atom_num*tsteps)

! Then, combine from all processes
	call MPI_REDUCE(tmpavg,restmp,1,MPI_DOUBLE_PRECISION,//&
       &	MPI_SUM,0,MPI_COMM_WORLD,ierr)

	if (myid.eq.0) then
! Now we have restmp with the sum of each averaged temperature
	   restmp = restmp/numprocs

! This might be its own section later, but just write the results out
	   open(unit=21,file='out.dat',position='append')
	   write(21,*) restmp, U_in

	end if

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
