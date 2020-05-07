!=======================================================================
! File mbp_data_mod.f90 contains module variables that are used to
! specify the properties of individual particles.
!
! The contained subroutine is used to initialize the arrays.
!
! Author: Carl Sovinec
! Modified: Vincent Bensch
!=======================================================================

  MODULE mbp_data_mod
  USE mbp_inp_mod
  USE mbp_kind_mod
  IMPLICIT NONE

  REAL(rknd), DIMENSION(:), POINTER :: rptr     ! Pointer for cartesian position of current particle
  REAL(rknd), DIMENSION(:), POINTER :: vptr     ! Pointer for cartesian velocity of current particle

  REAL(rknd), DIMENSION(:), POINTER :: drptr    ! Pointer for cartesian position derivative of current particle
  REAL(rknd), DIMENSION(:), POINTER :: dvptr    ! Pointer for cartesian velocity derivative of current particle

  REAL(rknd), DIMENSION(:), POINTER :: pre_rptr ! Pointer for radius of previous element
  REAL(rknd), DIMENSION(:), POINTER :: pre_vptr ! Pointer for velocity of previous element

  REAL(rknd), DIMENSION(:), POINTER :: rel_rptr ! Pointer for radius relative to previous element
  REAL(rknd), DIMENSION(:), POINTER :: rel_vptr ! Pointer for velocity relative to previous element

  REAL(rknd), DIMENSION(:), ALLOCATABLE, TARGET :: sln_vec   !Solution vector

  REAL(rknd), DIMENSION(2) :: sys_pe !System potential energy (previous, current)
  REAL(rknd), DIMENSION(2) :: sys_ke !System potential energy (previous, current)

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

  num_eqs = nelem * 4

  ALLOCATE(elem_mass(nelem))
  ALLOCATE(elem_rad(nelem))
  ALLOCATE(elem_theta(nelem))
  ALLOCATE(elem_theta_dot(nelem))

  ALLOCATE(sln_vec(1:num_eqs))
  ALLOCATE(mbp_ATOL(num_eqs))
  allocate(p_rptr(2))
  allocate(p_vptr(2))

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
  
  INTEGER(iknd) :: ielem, ioff
  REAL(rknd), DIMENSION(:), POINTER :: tolptr !Pointer for tolerance of current particle
  REAL(rknd) :: running_theta     = 0_rknd
  REAL(rknd) :: running_theta_dot = 0_rknd

  CALL reset_ptrs

!-----------------------------------------------------------------------
! Loop over the elements and use the input specifications or default
! values to assign the data arrays.
!-----------------------------------------------------------------------

  DO ielem=1,nelem

    running_theta = runnng_theta + elem_theta(ielem)
    ioff = 4 * ielem

    rptr => sln_vec(ioff-3:ioff-2)  !Position pointer associated with current particle
    vptr => sln_vec(ioff-1:ioff)    !Velocity pointer associated with current particle
    tolptr => mbp_ATOL(ioff-3:ioff) !Tolerance pointer associated with current particle

    rptr(1) = cos(running_theta) * elem_rad + p_rprt(1)!X
    rptr(2) = sin(running_theta) * elem_rad + p_rprt(2)!Y

    vptr(1) = (-1_rknd) * elem_theta_dot(ielem) * elem_rad * sin(runing_theta) + p_vprt(1)!X
    vptr(2) = elem_theta_dot(ielem) * elem_rad * cos(runing_theta) + p_vprt(1)!Y

    tolptr(1:2) = pos_tolerance !Position tolearnce
    tolptr(3:4) = vel_tolerance !Velocity tolerance

    CALL update_pre_ptrs

  ENDDO

  END SUBROUTINE mbp_data_init

!=======================================================================
! The module subroutine mbp_data_term reparses the final state into 
! the original _init variables and writes them to the output namelist
!=======================================================================

  SUBROUTINE mbp_data_term

  INTEGER(iknd) :: ielem, ioff

  CALL reset_ptrs

!-----------------------------------------------------------------------
! Loop over the particles and reset IC variable from solution vector
!-----------------------------------------------------------------------

  DO ielem=1,nelem

    running_theta = runnng_theta + elem_theta(ielem)
    ioff = 4 * ielem

    rptr => sln_vec(ioff-3:ioff-2)  !Position pointer relative to the previous particle
    vptr => sln_vec(ioff-1:ioff)    !Position pointer relative to the previous particle

    CALL update_rel_ptrs

    elem_theta = atan2(rptr(2), rptr(1))
    elem_theta_dot = (rptr(1) * vptr(2) - vptr(1) * rptr(2)) /&
                   & (rptr(1) ** 2 + rptr(2) ** 2)

    CALL update_pre_ptrs

  ENDDO

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
  PRINT *, "|", repeat(" ", 10), &
    "Welcome to Vincent's multi-body pendulum integrator program", &
    repeat(" ", 10), "|"
  PRINT *, "|", repeat(" ", 20), &
    "Written for the Final Project of EP-476", &
    repeat(" ", 20), "|"
  PRINT *, "|", repeat(" ", 27), &
    "Last updated 2020-05-06", &
    repeat(" ", 27), "|"
  PRINT *, repeat("=", 79), NEW_LINE('A')

  RETURN
  END SUBROUTINE mbp_welcome

!=======================================================================
! The module subroutine reset_ptrs resets the solution vector pointers to 0
!=======================================================================

  SUBROUTINE reset_ptrs

  rptr = 0_rknd
  vptr = 0_rknd

  pre_rptr = 0_rknd
  pre_vptr = 0_rknd

  rel_rptr = 0_rknd
  rel_vptr = 0_rknd

  END SUBROUTINE reset_ptrs

!=======================================================================
! The module subroutine update_pre_ptr updates the relative solution vector pointers
!=======================================================================

  SUBROUTINE update_pre_ptrs

  pre_rptr = rptr
  pre_vptr = vptr

  END SUBROUTINE update_pre_ptrs

!=======================================================================
! The module subroutine update_rel_ptr updates the relative solution vector pointers
!=======================================================================

  SUBROUTINE update_rel_ptrs

  rel_rptr = rptr - pre_rptr
  rel_vptr = vptr - pre_rptr

  END SUBROUTINE update_rel_ptrs

  END MODULE mbp_data_mod