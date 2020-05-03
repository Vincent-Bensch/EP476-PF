!=======================================================================
! File mbp_rk_int.f90 contains student-built second-order and fourth-order
! Runge-Kutta numerical integrators. They are general purpose subroutines 
! that has been created for a charged partical solver, but should be 
! generally applicable
!
! Author: Vincent Bensch
!========================================================================

SUBROUTINE RK2(F, NEQ, Y, T, TOUT, NUM_STEPS) !Second order Runge-Kutta
USE mbp_kind_mod
IMPLICIT NONE

EXTERNAL :: F
INTEGER(iknd), INTENT(IN) :: NUM_STEPS, NEQ
REAL(rknd), INTENT(IN) :: T, TOUT
REAL(rknd), DIMENSION(NEQ), INTENT(INOUT) :: Y
REAL(rknd), DIMENSION(NEQ) :: k1, k2, k3
REAL(rknd) :: tstep
REAL(rknd) :: htstep
INTEGER(iknd) :: istep
REAL(rknd) :: t_curr

tstep = (TOUT-T)/ NUM_STEPS
htstep = tstep / 2._rknd
t_curr = T

DO istep=1,NUM_STEPS
  CALL F( NEQ, t_curr, Y, k1)
  k2 = Y + htstep * k1
  CALL F( NEQ, t_curr + htstep, k2, k3)
  Y = Y + tstep * k3
  t_curr = t_curr + tstep
ENDDO

RETURN
END SUBROUTINE RK2


SUBROUTINE RK4(F, NEQ, Y, T, TOUT, NUM_STEPS) !Fourth order Runge-Kutta
USE mbp_kind_mod
IMPLICIT NONE

EXTERNAL :: F
INTEGER(iknd), INTENT(IN) :: NUM_STEPS, NEQ
REAL(rknd), INTENT(IN) :: T, TOUT
REAL(rknd), DIMENSION(NEQ), INTENT(INOUT) :: Y
REAL(rknd), DIMENSION(NEQ) :: k1, k2, k3
REAL(rknd), DIMENSION(NEQ) :: d1, d2, d3, d4
REAL(rknd) :: tstep
REAL(rknd) :: htstep
INTEGER(iknd) :: istep
REAL(rknd) :: t_curr

tstep = (TOUT-T)/ NUM_STEPS
htstep = tstep / 2._rknd
t_curr = T

DO istep=1,NUM_STEPS
  CALL F( NEQ, t_curr, Y, k1)
  d1 = k1 * tstep
  
  CALL F( NEQ, t_curr + htstep, Y + d1/2., k1)
  d2 = k1 * tstep  
  
  CALL F( NEQ, t_curr + htstep, Y + d2/2., k1)
  d3 = k1 * tstep  

  CALL F( NEQ, t_curr + tstep, Y + d3, k1)
  d4 = k1 * tstep  
  
  Y = Y + (d1 + 2*d2 + 2*d3 + d4) / 6.
  
  t_curr = t_curr + tstep
ENDDO

RETURN
END SUBROUTINE RK4