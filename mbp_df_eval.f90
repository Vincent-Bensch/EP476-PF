!=======================================================================
! File mbp_df_eval.f90 contains the mbp_df_eval function. It is to
! be passed to the relevant integrator function by mbp_main. In the current
! configuration it cannot handle time dependancy. The time vatiable is
! passed for forwards compatability.
!
! Author: Vincent Bensch
!========================================================================



  SUBROUTINE mbp_df_eval( nsize, time, svec, dfvec )
  
  USE mbp_inp_mod
  USE mbp_kind_mod
  USE mbp_data_mod
  IMPLICIT NONE
  INTEGER(iknd) :: ipart, ioff !current particle under evalutation

  INTEGER(iknd), INTENT(IN) :: nsize ! size of the ODE system
  REAL(rknd), INTENT(IN) :: time ! current value of ind. variable
  REAL(rknd), DIMENSION(nsize), TARGET, INTENT(IN) :: svec ! current ‘solution’ vector
  REAL(rknd), DIMENSION(nsize), TARGET, INTENT(OUT) :: dfvec ! will hold the derivative-function vector upon return
  

  DO ipart=1,npart !Loop though particles
 
!-----------------------------------------------------------------------
! Index pointers are advanced.
!-----------------------------------------------------------------------
    ioff = 6 * ipart
    rptr => svec(ioff-5:ioff-3) !Position pointer associated with current particle
    vptr => svec(ioff-2:ioff)   !Velocity pointer associated with current particle
    drptr => dfvec(ioff-5:ioff-3) !Derivative of position pointer associated with current particle
    dvptr => dfvec(ioff-2:ioff)  !Derivative of velocity pointer associated with current particle

!-----------------------------------------------------------------------
! The derrivative function of position is set the current 'solution' velocity
!-----------------------------------------------------------------------

    drptr(1) = vptr(1)
    drptr(2) = vptr(2)
    drptr(3) = vptr(3)

!-----------------------------------------------------------------------
! Evaluate the acceleration at ipart from constant electic and magnetic fields
!-----------------------------------------------------------------------

    dvptr(1) = qom(ipart) * (ex + vptr(2) * bz - vptr(3) * by)
    dvptr(2) = qom(ipart) * (ey + vptr(3) * bx - vptr(1) * bz)
    dvptr(3) = qom(ipart) * (ez + vptr(1) * by - vptr(2) * bx)

  ENDDO
  
  RETURN
  
  END SUBROUTINE mbp_df_eval