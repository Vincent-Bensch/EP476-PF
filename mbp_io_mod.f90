!=======================================================================
! This file contains the mbp_io_mod module which defines 
! input/output unit numbers with mnemonics to avoid confusion.
!=======================================================================

  MODULE mbp_io_mod
    USE mbp_kind_mod
    IMPLICIT NONE


!  Namelist filenames

    CHARACTER(32), PARAMETER :: nlparam_file  = "mbp_param.dat"  !Parameter input
    CHARACTER(32), PARAMETER :: nlin_file     = "mbp_input.dat"  !IC input
    CHARACTER(32), PARAMETER :: nlout_file    = "mbp_output.dat" !FC output
    CHARACTER(32), PARAMETER :: time_file     = "mbp_time.dat"   !Time mesh for plot
    CHARACTER(32), PARAMETER :: angle_file    = "mbp_anlge.dat"  !Angle mesh for plot
  
!  Unit numbers

    INTEGER(iknd), PARAMETER :: nlparam_unit  = 10 !Parameter input
    INTEGER(iknd), PARAMETER :: nlin_unit     = 11 !IC input
    INTEGER(iknd), PARAMETER :: nlout_unit    = 12 !FC output
    INTEGER(iknd), PARAMETER :: time_unit     = 13 !Time mesh for plot
    INTEGER(iknd), PARAMETER :: angle_unit    = 14 !x coord mesh for plot

  END MODULE mbp_io_mod
