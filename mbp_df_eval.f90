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
  REAL(rknd), DIMENSION(nsize), TARGET, INTENT(IN) :: svec ! current ‘solution’ vector
  REAL(rknd), DIMENSION(nsize), TARGET, INTENT(OUT) :: dfvec ! will hold the derivative-function vector upon return

  REAL(rknd) :: lt_prev = 0_rknd
  REAL(rknd) :: lt_next = 0_rknd

  CALL reset_ptrs



  RETURN

  END SUBROUTINE mbp_df_eval