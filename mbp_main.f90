!=======================================================================
! File mbp_main.f90 contains the primary driver for the
! EP 476 Charged-Particle Simulator (Application 1)
!
! Author: Vincent Bensch
!========================================================================

PROGRAM mbp_main

USE mbp_kind_mod
USE mbp_inp_mod
USE mbp_data_mod

IMPLICIT NONE

EXTERNAL :: mbp_df_eval !Externally defined derivative function for passing to integrators

INTEGER, DIMENSION(:), ALLOCATABLE :: mbp_iwork
REAL(rknd), DIMENSION(:), ALLOCATABLE :: mbp_rwork

INTEGER(iknd) :: mbp_istate = i_1 !DLSODE state, input 1 = first call
INTEGER(iknd) :: mbp_LRW, mbp_LIW, mbp_MF, mbp_JAC !length of rwork and iwork

CALL mbp_welcome

PRINT *, "Starting data loading", NEW_LINE('A')

CALL mbp_data_load

PRINT *, "Starting data initialization", NEW_LINE('A')

CALL mbp_data_init

SELECT CASE (integrator)
    CASE ("RK2")
      PRINT *, "Starting integration loop with 2nd order Runge-Kutta integrator",&
      NEW_LINE('A')
      CALL RK2(mbp_df_eval, num_eqs, sln_vec, t_initial, t_final, nstep)

    CASE ("RK4")
      PRINT *, "Starting integration loop with 4th order Runge-Kutta integrator",&
                NEW_LINE('A')
      CALL RK4(mbp_df_eval, num_eqs, sln_vec, t_initial, t_final, nstep)

    CASE ("ADAMS")
      PRINT *, "Starting integration loop with ODEPACK, DLSODE, non-stiff Adams"

      mbp_LIW = 20
      mbp_LRW = 20 + 16 * num_eqs

      PRINT *, "LIW: ", mbp_LIW, "LRW: ", mbp_LRW, "set"

      ALLOCATE(mbp_rwork(mbp_LRW)) !Allocate rworks size per docs

      PRINT *, "rwork allocated"

      mbp_rwork(5:10) = r_0 !Zeros in optional input band of rwork

      PRINT *, "rwork set"

      ALLOCATE(mbp_iwork(mbp_LIW)) !Allocate iworksize per docs

      PRINT *, "iwork allocated"

      mbp_iwork(5:10) = i_0 !Zeros in optional input band of rwork
      mbp_iwork(6) = nstep !Set maxteps

      PRINT *, "iwork set"

      mbp_MF = 10 !Set abrams

      PRINT *, "All parameters set, calling integrator"

      CALL DLSODE_caller

      PRINT *, "Number of Jacobi evaluations: ", mbp_iwork(13)
      PRINT *, "Number of Internal steps taken: ", mbp_iwork(11)
      PRINT *, "Number of Derrivative function evaluations: ", mbp_iwork(12)
      PRINT *, "Integrator exited with code: ", mbp_istate, NEW_LINE('A')

    CASE ("BDF")
      PRINT *, "Starting integration loop with ODEPACK, DLSODE, BDF"

      mbp_LIW = 20 + num_eqs
      mbp_LRW =  22 + 9 * num_eqs + num_eqs ** 2

      PRINT *, "LIW: ", mbp_LIW, "LRW: ", mbp_LRW, "set"

      ALLOCATE(mbp_rwork(mbp_LRW)) !Allocate rworks size per docs

      PRINT *, "rwork allocated"

      mbp_rwork(5:10) = r_0 !Zeros in optional input band of rwork

      PRINT *, "rwork set"

      ALLOCATE(mbp_iwork(mbp_LIW)) !Allocate iworksize per docs

      PRINT *, "iwork allocated"

      mbp_iwork(5:10) = i_0 !Zeros in optional input band of rwork
      mbp_iwork(6) = nstep !Set maxteps

      PRINT *, "iwork set"

      mbp_MF = 22 !Set BDF

      PRINT *, "All parameters set, calling integrator", mbp_istate, NEW_LINE('A')

      CALL DLSODE_caller

      PRINT *, "Number of Jacobi evaluations: ", mbp_iwork(13)
      PRINT *, "Number of Internal steps taken: ", mbp_iwork(11)
      PRINT *, "Number of Derrivative function evaluations: ", mbp_iwork(12)
      PRINT *, "Integrator exited with code: ", mbp_istate, NEW_LINE('A')


    CASE DEFAULT
      PRINT *, "Invalid integrator selection", NEW_LINE('A')
      PRINT *, NEW_LINE('A'), repeat("-", 35), " Goodbye ", &
               repeat("-", 36), repeat(NEW_LINE('A'), 5)
      RETURN
END SELECT



!Loop with integrator here

CALL mbp_data_term

PRINT *, NEW_LINE('A'), repeat("-", 35), " Goodbye ", &
  repeat("-", 36), repeat(NEW_LINE('A'), 5)
  
RETURN

CONTAINS

SUBROUTINE DLSODE_caller
  PRINT *, "External Call"
  CALL DLSODE (mbp_df_eval,&		!F = df function passed
    num_eqs,&				!NEQ = number of equations in system
    sln_vec,&				!Y = solution vector at ICs
    t_initial,&				!T = time at computation start
    t_final,&				!TOUT = time at computation end
    2_iknd,&				!ITOL = tol setting per problem statement
    rel_tolerance,&			!RTOL = user specified relative tolerance
    mbp_ATOL,&				!ATOL = absolute tolerance vector specified in data_mod
    i_1,&					!ITASK = 1 - normal computation
    mbp_istate,&			!ISTATE = 1 - first call
    i_1,&					!IOPT = 1 - optional input given
    mbp_rwork,&				!RWORK Provided real work array
    mbp_LRW,&				!LRW = len rwork
    mbp_iwork,&				!IWORK = provded int work array plus optional max-step
    mbp_LIW,&				!LIW = len iwork
    mbp_JAC,&				!JAC = Pass dummy
    mbp_MF)					!MF = specify type
END SUBROUTINE DLSODE_caller

END PROGRAM mbp_main