!=======================================================================
! File mbp_df_eval.f90 contains the mbp_df_eval function. It is to
! be passed to the relevant integrator function by mbp_main.
!
! Author: Vincent Bensch
!========================================================================

  SUBROUTINE mbp_df_eval( nsize, time, svec, dfvec )
  
  USE mbp_inp_mod
  USE mbp_kind_mod
  USE mbp_data_mod
  IMPLICIT NONE
  INTEGER(iknd) :: ielem, ioff !current element under evalutation

  INTEGER(iknd), INTENT(IN) :: nsize ! size of the ODE system
  REAL(rknd), INTENT(IN) :: time ! current value of ind. variable
  REAL(rknd), DIMENSION(nsize), TARGET, INTENT(IN) :: svec! current ‘solution’ vector
  REAL(rknd), DIMENSION(nsize), TARGET, INTENT(OUT) :: dfvec ! will hold the derivative-function vector upon return

  REAL(rknd) :: C1, TD, N1, N2, STD, CTD

  REAL(rknd), DIMENSION(:), POINTER :: theta_in
  REAL(rknd), DIMENSION(:), POINTER :: theta_dot_out
  REAL(rknd), DIMENSION(:), POINTER :: omega_in
  REAL(rknd), DIMENSION(:), POINTER :: omega_dot_out

  theta_in => svec(1:4:2)
  omega_in => svec(2:4:2)

  theta_dot_out => dfvec(1:4:2)
  omega_dot_out => dfvec(2:4:2)

  theta_dot_out = omega_in

  TD = theta_in(1) - theta_in(2)
  STD = sin(TD)
  CTD = cos(TD)

  C1 = 2 * elem_mass(1) + elem_mass(2) - elem_mass(2) * cos(2 * TD)

  N1 = -grav_accel * (2 * elem_mass(1) + elem_mass(2) ) * sin(theta_in(1)) &
       - elem_mass(2) * grav_accel * sin(theta_in(1) - 2 * theta_in(2)) &
       - 2 * STD * elem_mass(2) * (omega_in(2) ** 2 * elem_rad(2) &
       + omega_in(1) ** 2 * elem_rad(1) * CTD)

  N2 = 2 * STD * (omega_in(1) ** 2 * elem_rad(1) * (elem_mass(1) + elem_mass(2)) &
       + grav_accel * (elem_mass(1) + elem_mass(2)) * cos(theta_in(1)) &
       + omega_in(2) ** 2 * elem_rad(2) * elem_mass(2) * CTD)

  omega_dot_out(1) = N1 / (elem_rad(1) * C1)
  omega_dot_out(2) = N2 / (elem_rad(2) * C1)

  RETURN

  END SUBROUTINE mbp_df_eval