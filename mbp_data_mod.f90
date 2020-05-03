!=======================================================================
! File chp_data_mod.f90 contains module variables that are used to
! specify the properties of individual particles.
!
! The contained subroutine is used to initialize the arrays.
!
! Author: Carl Sovinec
! Modified: Vincent Bensch
!=======================================================================

  MODULE chp_data_mod
  USE chp_inp_mod
  USE chp_kind_mod
  IMPLICIT NONE

  REAL(rknd), DIMENSION(:), ALLOCATABLE :: pchrg  ! charge of each
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: pmass  ! mass of each
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: qom    ! charge to mass ratio
  
  REAL(rknd), DIMENSION(:), ALLOCATABLE, TARGET :: sln_vec	!Solution vector
  
  REAL(rknd), DIMENSION(:), POINTER :: rptr	!Pointer for position of current particle
  REAL(rknd), DIMENSION(:), POINTER :: vptr	!Pointer for velocity of current particle
 
  REAL(rknd), DIMENSION(:), POINTER :: drptr	!Pointer for position derivative of current particle
  REAL(rknd), DIMENSION(:), POINTER :: dvptr	!Pointer for velocity derivative of current particle
  
  INTEGER(iknd) :: num_eqs !Number of equations in system, and len of sln_vec
  REAL(rknd), DIMENSION(:), ALLOCATABLE, TARGET :: chp_ATOL !Absolute tolerance vector
  REAL(rknd) :: B_mag !Total magnetic field magnitude

  CONTAINS

!=======================================================================
! The module subroutine chp_data_init initializes the data arrays for
! the collection of particles.
!=======================================================================

  SUBROUTINE chp_data_init

  INTEGER(iknd) :: ipart, ioff
  LOGICAL  :: param_nml_exists, part_nml_in_exists !Bools to check exsistance of namelists
  REAL(rknd), DIMENSION(:), POINTER :: tolptr	!Pointer for tolerance of current particle
  REAL(rknd) :: v_part !Velocity of current particle
  
!-----------------------------------------------------------------------
! Read parateter namelist (filename specified in chp_inp_mod)
!-----------------------------------------------------------------------
  
  INQUIRE(FILE=param_nml_file, EXIST=param_nml_exists)
  INQUIRE(FILE=part_nml_file_in, EXIST=part_nml_in_exists)

  IF (param_nml_exists) THEN
	OPEN (UNIT=param_nml, FILE= param_nml_file, STATUS= "OLD", FORM= "FORMATTED")
    PRINT *, "Namelist opened: ", param_nml_file
	READ (UNIT=param_nml, NML= chp_param_nml)
	CLOSE (param_nml)
	PRINT *, "Namelist read into memory: ", param_nml_file, NEW_LINE('A')
  ELSE
    PRINT *, "Namelist not found: ", param_nml_file, NEW_LINE('A')
  ENDIF
  
!-----------------------------------------------------------------------
! Allocate the arrays for the number of particles in the computation.
!-----------------------------------------------------------------------

  num_eqs = npart * 6

  ALLOCATE(pchrg(npart))
  ALLOCATE(pmass(npart))
  ALLOCATE(qom(npart))
  ALLOCATE(sln_vec(1:num_eqs))
  
  ALLOCATE(x_init(npart))
  ALLOCATE(y_init(npart))
  ALLOCATE(z_init(npart))
  ALLOCATE(vx_init(npart))
  ALLOCATE(vy_init(npart))
  ALLOCATE(vz_init(npart))
  ALLOCATE(num_elecs(npart))
  ALLOCATE(num_prots(npart))
  ALLOCATE(num_neuts(npart))

  ALLOCATE(chp_ATOL(num_eqs))
!-----------------------------------------------------------------------
! Set ICs for allocatable arrays
!-----------------------------------------------------------------------

x_init=r_0  !  Initial value: 0
y_init=r_0  !  Initial value: 0
z_init=r_0  !  Initial value: 0
vx_init=r_0  !  Initial value: 0
vy_init=r_0  !  Initial value: 0
vz_init=r_0  !  Initial value: 0
num_elecs=i_1  !  Initial value: 1
num_prots=i_0  !  Initial value: 0
num_neuts=i_0  !  Initial value: 0

B_mag = sqrt(bx ** 2 + by ** 2 + bz ** 2)
   
!-----------------------------------------------------------------------
! Read array namelist (filename specified in chp_inp_mod)
!-----------------------------------------------------------------------
  
  IF (part_nml_in_exists) THEN
    OPEN (UNIT=part_nml_in, FILE= part_nml_file_in, STATUS= "OLD", FORM= "FORMATTED")
    PRINT *, "Namelist opened: ", part_nml_file_in
	READ (UNIT=part_nml_in, NML= chp_part_nml)
    CLOSE (part_nml_in)
	PRINT *, "Namelist read into memory: ", part_nml_file_in, NEW_LINE('A')
  ELSE
    PRINT *, "Namelist not found: ", part_nml_file_in, NEW_LINE('A')
  ENDIF  
  
!-----------------------------------------------------------------------
! Loop over the particles and use the input specifications or default
! values to assign the data arrays.
!-----------------------------------------------------------------------

  DO ipart=1,npart

    pchrg(ipart)=elem_chrg*(num_prots(ipart)-num_elecs(ipart))
    pmass(ipart)=proton_mass*num_prots(ipart) &
                +neutron_mass*num_neuts(ipart) &
                +electron_mass*num_elecs(ipart)
    qom(ipart)=pchrg(ipart)/pmass(ipart)
	
	ioff = 6 * ipart
	rptr => sln_vec(ioff-5:ioff-3) !Position pointer associated with current particle
	vptr => sln_vec(ioff-2:ioff)   !Velocity pointer associated with current particle
	tolptr => chp_ATOL(ioff-5:ioff)    !Tolerance pointer associated with current particle
	
	rptr(1) = x_init(ipart) !Position initial conditions
	rptr(2) = y_init(ipart)
	rptr(3) = z_init(ipart)
	
	vptr(1) = vx_init(ipart) !Velicity initial conditions
	vptr(2) = vy_init(ipart)
	vptr(3) = vz_init(ipart)
	
	v_part = sqrt(vptr(1) ** 2 + vptr(2) ** 2 + vptr(3) ** 2) !Particle velocity
	
	tolptr(1:3) = tolerance * ABS( (v_part * pmass(ipart)) /  (pchrg(ipart) / B_mag) )!Position tolearnce
	tolptr(4:6) = tolerance * v_part!Velocity tolerance

  ENDDO
  RETURN
  END SUBROUTINE chp_data_init
  
  
!=======================================================================
! The module subroutine chp_data_term reparses the final state into 
! the original _init variables and writes them to the output namelist
!=======================================================================

  SUBROUTINE chp_data_term
  
  INTEGER(iknd) :: ipart, ioff

!-----------------------------------------------------------------------
! Loop over the particles and reset IC variable from solution vector
!-----------------------------------------------------------------------

  DO ipart=1,npart
	
	ioff = 6 * ipart
	rptr => sln_vec(ioff-5:ioff-3) !Position pointer associated with current particle
	vptr => sln_vec(ioff-2:ioff)   !Velocity pointer associated with current particle
	
	x_init(ipart) = rptr(1) !Position initial conditions
	y_init(ipart) = rptr(2) 
	z_init(ipart) = rptr(3) 
	
	vx_init(ipart) = vptr(1) !Velicity initial conditions
	vy_init(ipart) = vptr(2)
	vz_init(ipart) = vptr(3)

  ENDDO  
  
  OPEN (UNIT=part_nml_out, FILE= part_nml_file_out, STATUS= "REPLACE", FORM= "FORMATTED")
  WRITE (UNIT=part_nml_out, NML= chp_part_nml)
  PRINT *, "Namelist written from memory: ", part_nml_file_out
  CLOSE (part_nml_out)
  RETURN
  
  END SUBROUTINE chp_data_term
 
!=======================================================================
! The module subroutine chp_welcome prints a welcome message to screen
!=======================================================================
 
  SUBROUTINE chp_welcome
  
  PRINT *, repeat(NEW_LINE('A'), 5)
  PRINT *, repeat("=", 79)
  PRINT *, "|", repeat(" ", 11), &
    "Welcome to Vincent's charged particle integrator program", &
    repeat(" ", 10), "|"
  PRINT *, "|", repeat(" ", 23), &
    "Written for Project 4 of EP-476", &
    repeat(" ", 23), "|"
  PRINT *, "|", repeat(" ", 27), &
    "Last updated 2020-04-10", &
    repeat(" ", 27), "|"
  PRINT *, repeat("=", 79), NEW_LINE('A')
  
  RETURN
  END SUBROUTINE chp_welcome
  
  END MODULE chp_data_mod
