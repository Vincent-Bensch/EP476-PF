!=======================================================================
! File mbp_knd_mod.f90 contains basic parameters for variable definition
! in other modules.
!
! Author: Vincent Bensch
!========================================================================

 MODULE mbp_kind_mod
      IMPLICIT NONE

      INTEGER, PARAMETER :: rknd=SELECTED_REAL_KIND(14, 100)
      INTEGER, PARAMETER :: iknd=SELECTED_INT_KIND(8)

      REAL(KIND = rknd), PARAMETER :: r_0=0._rknd
      REAL(KIND = rknd), PARAMETER :: r_1=1._rknd
      REAL(KIND = rknd), PARAMETER :: r_4=4._rknd

      INTEGER(KIND = iknd), PARAMETER :: i_0=0_iknd
      INTEGER(KIND = iknd), PARAMETER :: i_1=1_iknd
END MODULE mbp_kind_mod