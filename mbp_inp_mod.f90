!=======================================================================
! File mbp_inp_mod.f90 contains all variables that are set by the
! user at run-time.  There are scalar variables for physical parameters
! such as particle masses, and the number of particles to be used in
! a computation.  There are array variables that are used to set the
! initial conditions of particle coordinates and velocity vectors
! for each particle.  There are also parameters that control the
! numerical integration and specifications for postprocessing.
!
! Default values are provided for each variable.  The main program
! will read them from a file through namelist read statements.
!
! Author: Carl Sovinec
! Modified: Vincent Bensch
!=======================================================================

  MODULE mbp_inp_mod
  USE mbp_kind_mod
  IMPLICIT NONE

! The first set defines physical parameters that may be changed to
! consider a different set of units.  Default values are in MKS units.

  REAL(rknd), PARAMETER  :: elem_chrg=1.602176634e-19_rknd
  REAL(rknd), PARAMETER  :: electron_mass=9.1093837015e-31_rknd
  REAL(rknd), PARAMETER  :: proton_mass=1.67262192369e-27_rknd
  REAL(rknd), PARAMETER  :: neutron_mass=1.67492749804e-27_rknd

! Set the externally applied uniform fields.

  REAL(rknd) :: ex=r_0  !  3 components of electric field, Initial value: 0
  REAL(rknd) :: ey=r_0  !  Initial value: 0
  REAL(rknd) :: ez=r_0  !  Initial value: 0
  REAL(rknd) :: bx=r_0  !  3 components of magnetic field, Initial value: 0
  REAL(rknd) :: by=r_0  !  Initial value: 0
  REAL(rknd) :: bz=r_0  !  Initial value: 0

! The Penning-trap potential is defined in terms of the applied
! potential (vtrap), its radius (rtrap), and its axial dimension (ztrap).

  REAL(rknd) :: vtrap=r_0  !  Initial value: 0
  REAL(rknd) :: rtrap=r_1  !  Initial value: 1
  REAL(rknd) :: ztrap=r_1  !  Initial value: 1

! The next set defines run-time specification of physical parameters
! for a particular computation.

  INTEGER(iknd) :: npart=i_1    !   # of particles in a computation, Initial value: 1
  REAL(rknd) :: t_initial=r_0   !  Initial value: 0
  REAL(rknd) :: t_final=r_1     !  Initial value: 1

! The next set defines initial positions and velocity vectors in
! Cartesian components.  The integer parameter npart_max is used to
! declare these initial-value arrays. Arrays
! that are used to define the type of each particle (number of sub-
! atomic particles) are also declared.

  REAL(rknd), DIMENSION(:), ALLOCATABLE :: x_init
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: y_init
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: z_init
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: vx_init
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: vy_init
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: vz_init
  INTEGER(iknd), DIMENSION(:), ALLOCATABLE :: num_elecs
  INTEGER(iknd), DIMENSION(:), ALLOCATABLE :: num_prots
  INTEGER(iknd), DIMENSION(:), ALLOCATABLE :: num_neuts

! The following set of parameters influence the choice and operation
! of the numerical method used to solve the ODE system.

! Two Runge-Kutta integrators are implemented.  The heuristic "RK2" or
! midpoint method is second-order accurate.  The conventional "RK4"
! option is a fourth-order method.  In the implementation here, both
! use fixed step sizes that are set by (t_final-t_initial)/nstep.

! Multi-step methods from the ODEPACK library are also available.
! Setting the integrator input to "Adams" calls the DLSODE routine
! from ODEPACK with MF set to 10 for the non-stiff Adams method.
! Setting integrator to "BDF" calls DLSODE with MF set to 22, which
! is the implicit backward difference formula for stiff problems.
! The Jacobi matrix for the latter is computed numerically while
! the computation is running.

! The tolerance input is used for specifying the desired accuracy
! when using the DLSODE routine.

  CHARACTER(32) :: integrator="RK4"
  INTEGER(iknd) :: nstep=i_1  !  Initial value: 1
  REAL(rknd) :: tolerance=1.e-9_rknd

! The t_plot is the time interval between writing output for
! visualization.  For single-step integration methods, the number
! of intervals is determined by (t_final-t_initial)/t_plot, but 
! the data is saved when a step matches or exceeds a t_plot interval.

  REAL(rknd) :: t_plot=r_1  !  Initial value: 1

! Namelist units

  INTEGER(iknd) :: param_nml=1 !Parameter namelist input 
  INTEGER(iknd) :: part_nml_in=2 !Particle namelist input 
  INTEGER(iknd) :: part_nml_out=3 !Particle namelist output 
  
! Namelist filenames
  
  CHARACTER(32) :: param_nml_file="mbp_param_nml"  		!Parameter namelist input 
  CHARACTER(32) :: part_nml_file_in="mbp_part_nml_in"  	!Particle namelist input 
  CHARACTER(32) :: part_nml_file_out="mbp_part_nml_out" !Particle namelist output 
  
!Namelist variable assignments
  
  NAMELIST / mbp_param_nml / &			!Parameter namelist input 
      ex, ey, ez, &  					!Electic Field
      bx, by, bz, & 					!Magnetic Field
      vtrap, rtrap, ztrap, & 			!Trap parameters 
	  npart, t_initial, t_final, &		!Physical parameters
	  integrator, nstep, tolerance, &	!Integrator settings
	  t_plot							!Plot settings
	  
  NAMELIST / mbp_part_nml / &			!Particle namelist input/output
	  x_init, y_init, z_init, &			!Initial particle positions
	  vx_init, vy_init, vz_init, &		!Initial particle velocities
	  num_elecs, num_prots, num_neuts	!Particle composition

  END MODULE mbp_inp_mod
