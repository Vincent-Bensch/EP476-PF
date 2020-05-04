!=======================================================================
! File mbp_inp_mod.f90 contains all variables that are set by the
! user at run-time.  There are scalar variables for physical parameters
! such as pendulum masses, and the number of elements to be used in
! a computation.  There are array variables that are used to set the
! initial conditions of particle coordinates and velocity vectors
! for each element.  There are also parameters that control the
! numerical integration and specifications for postprocessing.
!
! Default values are provided for each variable.  The main program
! will read them from a file through namelist read statements.
!
! Author: Vincent Bensch
!=======================================================================

  MODULE mbp_inp_mod
  USE mbp_kind_mod
  IMPLICIT NONE

! The physical constant for gratitational acceleration is set to the
! World Geodetic System 1984 value. It can be changed by namelist input

  REAL(rknd) :: grav_accel = 9.80665_rknd

! Data for the fluid in which the pendulum moves. Default values for
! Dry air at STP (15 deg C and 101.3 kPa).

  REAL(rknd) :: rho_surr = 1.225_rknd   !Density in kg/m^3
  REAL(rknd) :: mu_surr  = 1.73e-5_rknd !Viscosity in Pa-s

! The next set defines run-time specification of physical parameters
! for a particular computation.

  INTEGER(iknd) :: npart     = i_1 ! # of elements, Initial value: 1
  REAL(rknd)    :: t_initial = r_0 ! Initial value: 0
  REAL(rknd)    :: t_final   = r_1 ! Initial value: 1

! The next set defines initial positions and velocity vectors in
! Cartesian components.  The integer parameter npart_max is used to
! declare these initial-value arrays. Arrays
! that are used to define the type of each particle (number of sub-
! atomic particles) are also declared.

  REAL(rknd), DIMENSION(:), ALLOCATABLE :: r
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: theta_init
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: theta_dot_init
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
