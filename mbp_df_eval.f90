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

  REAL(rknd), DIMENSION(nelem + 1, 3) :: t_elem = 0_rknd !Tention in linkage from element n to element n-1 (x_comp, ycomp, mag)

  REAL(rknd) :: lt_prev = 0_rknd
  REAL(rknd) :: lt_next = 0_rknd

  CALL reset_ptrs

  DO ielem=nelem, 1 !Loop though elements

!-----------------------------------------------------------------------
! Index pointers are advanced.
!-----------------------------------------------------------------------

    ioff = 4 * ielem
    rptr => svec(ioff-3:ioff-2)   !Position pointer associated with current element
    vptr => svec(ioff-1:ioff)     !Velocity pointer associated with current element

    drptr => dfvec(ioff-3:ioff-2) !Derivative of position pointer associated with current element
    dvptr => dfvec(ioff-1:ioff)   !Derivative of velocity pointer associated with current element

    rptr => svec(ioff-3:ioff-2)   !Position pointer associated with current element
    vptr => svec(ioff-1:ioff)     !Velocity pointer associated with current element

!-----------------------------------------------------------------------
! The derrivative function of position is set the current 'solution' velocity
!-----------------------------------------------------------------------

    drptr(1) = vptr(1)
    drptr(2) = vptr(2)

!-----------------------------------------------------------------------
! Compute linkage tention
!-----------------------------------------------------------------------

    

!-----------------------------------------------------------------------
! Evaluate the acceleration at ielem from gravitational vector, and linkages
!-----------------------------------------------------------------------

    dvptr(1) = (t_elem(ielem, 1) - t_elem(ielem - 1, 1)) / elem_mass(ielem)
    dvptr(2) = (t_elem(ielem, 2) - t_elem(ielem - 1, 2)) / elem_mass(ielem) - grav_accel

  ENDDO

  RETURN

  END SUBROUTINE mbp_df_eval