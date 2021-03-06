!=======================================================================
! File mbp_data_mod.f90 contains module variables that are used to
! specify the properties of individual elements.
!
! The contained subroutine is used to initialize the arrays.
!
! Author: Carl Sovinec
! Modified: Vincent Bensch
!=======================================================================

  MODULE mbp_data_mod
  USE mbp_inp_mod
  USE mbp_kind_mod
  USE mbp_io_mod

  IMPLICIT NONE

  REAL(rknd), DIMENSION(:), ALLOCATABLE, TARGET :: sln_vec   !Solution vector

  INTEGER(iknd) :: num_eqs !Number of equations in system, and len of sln_vec
  REAL(rknd), DIMENSION(:), ALLOCATABLE, TARGET :: mbp_ATOL !Absolute tolerance vector

  CONTAINS

!=======================================================================
! The module subroutine mbp_data_load loads namelist data, and allocates
! arrays
!=======================================================================

  SUBROUTINE mbp_data_load

  LOGICAL  :: nlparam_exsists, nlstate_in_exsists !Bools to check exsistance of namelists

!-----------------------------------------------------------------------
! Read parateter namelist (filename specified in mbp_inp_mod)
!-----------------------------------------------------------------------
  
  INQUIRE(FILE = nlparam_file, EXIST = nlparam_exsists)
  INQUIRE(FILE = nlin_file,    EXIST = nlstate_in_exsists)

  IF (nlparam_exsists) THEN
    OPEN (UNIT=nlparam_unit, FILE= nlparam_file, STATUS= "OLD", FORM= "FORMATTED")
    PRINT *, "Namelist opened: ", nlparam_file
    READ (UNIT=nlparam_unit, NML= nlparam)
    CLOSE (nlparam_unit)
    PRINT *, "Namelist read into memory: ", nlparam_file, NEW_LINE('A')
  ELSE
    PRINT *, "Namelist not found: ", nlparam_file, NEW_LINE('A')
  ENDIF

!-----------------------------------------------------------------------
! Allocate the arrays for the number of elements in the pendulum.
!-----------------------------------------------------------------------

  num_eqs = 2 * nelem

  ALLOCATE(elem_mass(nelem))
  ALLOCATE(elem_rad(nelem))
  ALLOCATE(elem_theta(nelem))
  ALLOCATE(elem_omega(nelem))

  ALLOCATE(sln_vec(1:num_eqs))
  ALLOCATE(mbp_ATOL(num_eqs))

!-----------------------------------------------------------------------
! Read array namelist (filename specified in mbp_inp_mod)
!-----------------------------------------------------------------------

  IF (nlstate_in_exsists) THEN
    OPEN (UNIT=nlin_unit, FILE= nlin_file, STATUS= "OLD", FORM= "FORMATTED")
    PRINT *, "Namelist opened: ", nlin_file
    READ (UNIT=nlin_unit, NML= nlstate)
    CLOSE (nlin_unit)
    PRINT *, "Namelist read into memory: ", nlin_file, NEW_LINE('A')
  ELSE
    PRINT *, "Namelist not found: ", nlin_file, NEW_LINE('A')
  ENDIF

  RETURN
  END SUBROUTINE mbp_data_load

!=======================================================================
! The module subroutine mbp_data_init sets initial conditions and
! tolerance vectors from loaded data
!=======================================================================

  SUBROUTINE mbp_data_init

!-----------------------------------------------------------------------
! Setup solution vector from ICs
!-----------------------------------------------------------------------

  sln_vec(1) = elem_theta(1)
  sln_vec(2) = elem_omega(1)
  sln_vec(3) = elem_theta(2)
  sln_vec(4) = elem_omega(2)

  mbp_ATOL(1) = theta_tolerance !Position tolearnce
  mbp_ATOL(3) = theta_tolerance !Position tolearnce
  mbp_ATOL(2) = omega_tolerance !Velocity tolerance
  mbp_ATOL(4) = omega_tolerance !Velocity tolerance

!Print energy
  PRINT *, "Total system energy: ", mbp_energy()

  END SUBROUTINE mbp_data_init

!=======================================================================
! The module subroutine mbp_data_term reparses the final state into 
! the original _init variables and writes them to the output namelist
!=======================================================================

  SUBROUTINE mbp_data_term

!Print energy
  PRINT *, "Total system energy: ", mbp_energy()

!-----------------------------------------------------------------------
! Create fcs from solution vector
!-----------------------------------------------------------------------

  elem_theta(1) = sln_vec(1)
  elem_theta(2) = sln_vec(3)

  elem_omega(1) = sln_vec(2)
  elem_omega(2) = sln_vec(4)

  OPEN (UNIT=nlout_unit, FILE= nlout_file, STATUS= "REPLACE", FORM= "FORMATTED")
  WRITE (UNIT=nlout_unit, NML= nlstate)
  PRINT *, "Namelist written from memory: ", nlout_file
  CLOSE (nlout_unit)
  RETURN

  END SUBROUTINE mbp_data_term

!=======================================================================
! The module subroutine mbp_welcome prints a welcome message to screen
!=======================================================================

  SUBROUTINE mbp_welcome

  PRINT *, repeat(NEW_LINE('A'), 5)
  PRINT *, repeat("=", 79)
  PRINT *, "|", repeat(" ", 9), &
    "Welcome to Vincent's multi-body pendulum integrator program", &
    repeat(" ", 9), "|"
  PRINT *, "|", repeat(" ", 19), &
    "Written for the Final Project of EP-476", &
    repeat(" ", 19), "|"
  PRINT *, "|", repeat(" ", 27), &
    "Last updated 2020-05-07", &
    repeat(" ", 27), "|"
  PRINT *, repeat("=", 79), NEW_LINE('A')

  RETURN
  END SUBROUTINE mbp_welcome

!=======================================================================
! The module function mbp_energy returns system energy
!=======================================================================

  FUNCTION mbp_energy()

  REAL(rknd) mbp_energy, pot_energy, kin_energy
  REAL(rknd), DIMENSION(2) :: theta, omega

  theta = sln_vec(1:3:2)
  omega = sln_vec(2:4:2)

  pot_energy = -(elem_mass(1) + elem_mass(2)) * grav_accel * elem_rad(1) * cos(theta(1)) &
             - elem_mass(2) * grav_accel * elem_rad(2) * cos(theta(2)) ! Potential energy

  kin_energy = 1/2 * elem_mass(1) * elem_rad(1) ** 2 * omega(1) ** 2 &
             + 1/2 * elem_mass(2) * ( elem_rad(1) ** 2 * omega(1) ** 2 &
             + elem_rad(2) ** 2 * omega(2) ** 2 &
             + 2 * elem_rad(1) * elem_rad(2) * omega(1) * omega(2) * cos(theta(1) - theta(2) ))

  mbp_energy = pot_energy + kin_energy

  END FUNCTION mbp_energy

  END MODULE mbp_data_mod