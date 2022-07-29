! ************************************************    !
!		          Diffusion Monte Carlo						!
!		            Sample Program     						!
!				      												!
!   This program pefforms Monte Carlo Integration 		!
!   wave equation for several systems, including		!
!   the H.O., LJ. P., and for H, H2+, and H2				!

PROGRAM main		! main function
	
	IMPLICIT NONE
	
	! ********************************************************** !
	!	                 Variable Definitions				   		 !
	!																				 !	
	!	2 = replica matrix											 !
	!	ener1 = array used to store E0 values (for deviation)     !
	!	wz1 = 'boxes' used to sort replicas								 !
	!	itime = measures time												 !
	!	rf & ri = right & left measuring bounds for replica boxes !
	!	others = misc															 !
	!																				 !
	! ********************************************************** !
	
	
	! ***************************************************** !
	! 			    Help translations ( FORTRAN == C )			  !
	!																		  !
	!		double precision == double								  !
	!		integer(kind = 8) == long int                     !
	!		&& --> AND (C) 
	! ***************************************************** !
	
	! Paramaters
	
	REAL, PARAMETER :: crhh = 1.398
	REAL, PARAMETER :: Rhh2 = 2.00
	
	! POINTERS
	
	REAL, DIMENSION(:), ALLOCATABLE :: ener
	REAL, DIMENSION(:,:), ALLOCATABLE :: psips
	REAL, DIMENSION(:), ALLOCATABLE :: wz
	
	! VARIABLES DECLARATION
	
	DOUBLE PRECISION :: ri ! left limit for sampling data 
	DOUBLE PRECISION :: mrf ! right limit for sampling data
	DOUBLE PRECISION :: dt ! time step
	DOUBLE PRECISION :: deviation ! deviation from R 
	DOUBLE PRECISION :: rhh !  radius
	DOUBLE PRECISION :: Vref ! reference potential
	INTEGER :: status, stat
	INTEGER(KIND = 8) :: potnum ! choose the potential
	INTEGER(KIND = 8) :: Ndim ! number of dimensions we're working in (replicas)
	INTEGER(KIND = 8) :: Npsips ! Number of replicas
	INTEGER(KIND = 8) :: Nmax ! max number of replicas
	INTEGER(KIND = 8) :: seed ! ramdom number
	INTEGER(KIND = 8) :: Nchec ! Amount of time to run the simulation
	INTEGER(KIND = 8) :: Nwf ! number of boxes to sort replicas in
	INTEGER(KIND = 8) :: Ncol ! number of time steps to run the simulation
	
	!!!! FILE (line 34 on montecarlo.c)
	
	! *********** USER initialization *************** !
	
	CALL system('clear') ! clear console in terminal
	WRITE (*,*) '*****************************************************'
	WRITE (*,*) '*                                                   *'
	WRITE (*,*) '*           Diffusion Monte Carlo                   *'
	WRITE (*,*) '*                Simulator                          *'
	WRITE (*,*) '*                                                   *'
	WRITE (*,*) 'A.	One-dimensional potentials'
	WRITE (*,*) '		1. Harmonic Oscillator'
	WRITE (*,*) '		2. Morse Potential'
	WRITE (*,*) '-------------------------------'
	WRITE (*,*) 'B. 	3-dimensional potentials'
	WRITE (*,*) '		3. Hydrogen atom'
	WRITE (*,*) '		4. H2+ molecule'
	WRITE (*,*) '		5. H2 molecule'
	WRITE (*,*) '-------------------------------'
	WRITE (*,*) 'Select a potential to simulate:  '
	READ (*,*) potnum
	IF (potnum < 3) THEN ! (potnum = 2,1)
		Ndim = 1
	END IF
	IF ((potnum > 2) .AND. (potnum < 5)) THEN ! (potnum = 3,4)
		Ndim = 3
	END IF
	IF (potnum == 5) THEN ! (potnum = 5)
		Ndim = 6
	END IF
	WRITE (*,*) 'Number in parenthesis are suggestions'
	WRITE (*,*) 'Ndim = ', Ndim
	WRITE (*,*) 'Enter the Number of Replicas: (500)'
	READ (*,*) Npsips
	WRITE (*,*) 'Replicas = ', Npsips
	WRITE (*,*) 'Enter the max number of replicas: (2000)'
	READ (*,*) Nmax
	WRITE (*,*) 'Nmax = ', Nmax
	WRITE (*,*) 'Enter random number seed: '
	READ (*,*) seed
	IF (seed > 0) THEN
		seed = -1 * seed
	END IF
	WRITE (*,*) 'Seed = ', -1 * seed
	WRITE (*,*) 'Enter the amount of time to run the simulation:    (1000)'
	READ (*,*) Nchec
	IF (potnum /= 3) THEN
		WRITE (*,*) 'Enter the left limit for sampling data: (-20)'
		READ (*,*) ri
	ELSE 
		ri = 0
	END IF
	WRITE (*,*) 'left = ', ri
	WRITE (*,*) 'Enter the right limit for sampling data: (20) '
	READ (*,*) mrf
	WRITE (*,*) 'right = ', mrf
	WRITE (*,*) 'Enter the number of boxes to sort into: (200)'
	READ (*,*) Nwf
	Ncol = Nchec
	WRITE (*,*) 'boxes = ', Nwf
	WRITE (*,*) 'Enter a time step: (.1)'
	READ (*,*) dt
	
	IF ( (potnum < 6) .AND. (potnum > 3) ) THEN ! potnum = 4 e 5
		WRITE (*,*) 'Enter deviation from R ( H - H distance ) in units of bohr radius: '
		READ (*,*) deviation
		IF (potnum == 4) THEN 
			rhh = rhh2 + deviation
		END IF
		IF (potnum == 5) THEN
			rhh = crhh + deviation
			WRITE (*,*) 'Using', rhh, 'as hydrogen-hydrogen distance.'
		END IF
	END IF
	
	ALLOCATE(ener(1:Ncol), stat=status)
	IF (stat /= 0) THEN
		WRITE (*,*) "Something went wrong trying to allocate 'ener'"
		STOP 1
	END IF
	
	ALLOCATE(wz(1:Nwf), stat=status)
	IF (stat /= 0) THEN
		WRITE (*,*) "Something went wrong trying to allocate 'wz'"
		STOP 1
	END IF
	
	ALLOCATE(psips(1:Nmax,1:1+Ndim), stat=status)
	IF (stat /= 0) THEN
		WRITE (*,*) "Something went wrong trying to allocate 'psips'"
		STOP 1
	END IF
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! 
	! OPEN ALL FILES FOR OUTPUT !
	! 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
	! OPEN FILES !
	OPEN(1, FILE = 'Eo.dat')
	OPEN(2, FILE = 'waveshap.dat')
	OPEN(3, FILE = 'wavedump.dat')
	
	! FORMATING FILES !
	
	! Eo.dat !
	
	WRITE(1,*) "# Format:"
	WRITE(1,*) "# columm 1: time dimensionless units"
	WRITE(1,*) "# psips = ", Npsips, "  seed = ", seed, "  dt = ", dt, "  time = ", Nchec, "  left = ", ri
	WRITE(1,*) "# right = ", mrf, "  boxes = ", Nwf
	! CLOSE(1) ! Control
	
	!	waveshap.dat !
	
	WRITE(2,*) "# Format"
	WRITE(2,*) "# column 1: position (unitless for H.O, M.P., units of Bohr radius otherwise"
	WRITE(2,*) "# column 2: R (for hydrogen), in other cases the wave function normalized"
	WRITE(2,*) "# column 3: r * R (for hydrogen), does NOT occur otherwise"
	WRITE(2,*) "# Note:  For hydrogen we show the radial wavefunction, all other cases show a distribution"
	WRITE(2,*) "accross a single axis."
	WRITE(2,*) "# psips = ", Npsips, "  seed = ", seed, "  dt = ", dt, "  time = ", Nchec, "  left = ", ri
	WRITE(2,*) "# right = ", mrf, "  boxes = ", Nwf
	! CLOSE(2) ! Control
	
	! wavedump.dat !
	
	WRITE(3,*) "# this is a dump of the replica matrix at the end of the simulation."
	WRITE(3,*) "# x, y, z for all electrons"
	WRITE(3,*) "# psips = ", Npsips, "  seed = ", seed, "  dt = ", dt, "  time = ", Nchec, "  left = ", ri
	WRITE(3,*) "# right = ", mrf, "  boxes = ", Nwf
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	Vref = 0.
	
	WRITE (*,*) "Psips initialized"
	


	
	
	
	
	
	
END PROGRAM main
	
