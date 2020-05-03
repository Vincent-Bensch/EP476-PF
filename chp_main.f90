!=======================================================================
! File chp_main.f90 contains the primary driver for the
! EP 476 Charged-Particle Simulator (Application 1)
!
! Author: Vincent Bensch
!========================================================================

PROGRAM chp_main

USE chp_kind_mod
USE chp_inp_mod
USE chp_data_mod

IMPLICIT NONE

EXTERNAL :: chp_df_eval !Externally defined derivative function for passing to integrators

INTEGER, DIMENSION(:), ALLOCATABLE :: chp_iwork
REAL(rknd), DIMENSION(:), ALLOCATABLE :: chp_rwork

INTEGER(iknd) :: chp_istate = i_1 !DLSODE state, input 1 = first call
INTEGER(iknd) :: chp_LRW, chp_LIW, chp_MF, chp_JAC !length of rwork and iwork

CALL chp_welcome

PRINT *, "Starting data initialization", NEW_LINE('A')

CALL chp_data_init

SELECT CASE (integrator)
	CASE ("RK2")
		PRINT *, "Starting integration loop with 2nd order Runge-Kutta integrator",&
    		NEW_LINE('A')
		CALL RK2(chp_df_eval, num_eqs, sln_vec, t_initial, t_final, nstep)
		
	CASE ("RK4")
		PRINT *, "Starting integration loop with 4th order Runge-Kutta integrator",&
    		NEW_LINE('A')
		CALL RK4(chp_df_eval, num_eqs, sln_vec, t_initial, t_final, nstep)
		
	CASE ("ADAMS")
		PRINT *, "Starting integration loop with ODEPACK, DLSODE, non-stiff Adams"
		
		chp_LIW = 20
		chp_LRW = 20 + 16 * num_eqs
		
		PRINT *, "LIW: ", chp_LIW, "LRW: ", chp_LRW, "set"

		ALLOCATE(chp_rwork(chp_LRW)) !Allocate rworks size per docs
		
		PRINT *, "rwork allocated"
		
		chp_rwork(5:10) = r_0 !Zeros in optional input band of rwork
		
		PRINT *, "rwork set"
		
		ALLOCATE(chp_iwork(chp_LIW)) !Allocate iworksize per docs
		
		PRINT *, "iwork allocated"
		
		chp_iwork(5:10) = i_0 !Zeros in optional input band of rwork
		chp_iwork(6) = nstep !Set maxteps
		
		PRINT *, "iwork set"
		
		chp_MF = 10 !Set abrams
		
		PRINT *, "All parameters set, calling integrator"
		
		CALL DLSODE_caller
			
		PRINT *, "Number of Jacobi evaluations: ", chp_iwork(13)
		PRINT *, "Number of Internal steps taken: ", chp_iwork(11)
		PRINT *, "Number of Derrivative function evaluations: ", chp_iwork(12)
		PRINT *, "Integrator exited with code: ", chp_istate, NEW_LINE('A')
	
	CASE ("BDF")
		PRINT *, "Starting integration loop with ODEPACK, DLSODE, BDF"
		
		chp_LIW = 20 + num_eqs
		chp_LRW =  22 + 9 * num_eqs + num_eqs ** 2
		
		PRINT *, "LIW: ", chp_LIW, "LRW: ", chp_LRW, "set"

		ALLOCATE(chp_rwork(chp_LRW)) !Allocate rworks size per docs
		
		PRINT *, "rwork allocated"
		
		chp_rwork(5:10) = r_0 !Zeros in optional input band of rwork
		
		PRINT *, "rwork set"
		
		ALLOCATE(chp_iwork(chp_LIW)) !Allocate iworksize per docs
		
		PRINT *, "iwork allocated"
		
		chp_iwork(5:10) = i_0 !Zeros in optional input band of rwork
		chp_iwork(6) = nstep !Set maxteps
		
		PRINT *, "iwork set"
		
		chp_MF = 22 !Set BDF
		
		PRINT *, "All parameters set, calling integrator", chp_istate, NEW_LINE('A')
		
		CALL DLSODE_caller
			
		PRINT *, "Number of Jacobi evaluations: ", chp_iwork(13)
		PRINT *, "Number of Internal steps taken: ", chp_iwork(11)
		PRINT *, "Number of Derrivative function evaluations: ", chp_iwork(12)
		PRINT *, "Integrator exited with code: ", chp_istate, NEW_LINE('A')
		
		
	CASE DEFAULT
		PRINT *, "Invalid integrator selection", NEW_LINE('A')
		PRINT *, NEW_LINE('A'), repeat("-", 35), " Goodbye ", &
			repeat("-", 36), repeat(NEW_LINE('A'), 5)
		RETURN
END SELECT



!Loop with integrator here

CALL chp_data_term

PRINT *, NEW_LINE('A'), repeat("-", 35), " Goodbye ", &
  repeat("-", 36), repeat(NEW_LINE('A'), 5)
  
RETURN

CONTAINS

SUBROUTINE DLSODE_caller
  PRINT *, "External Call"
  CALL DLSODE (chp_df_eval,&		!F = df function passed
	num_eqs,&				!NEQ = number of equations in system
	sln_vec,&				!Y = solution vector at ICs
	t_initial,&				!T = time at computation start
	t_final,&				!TOUT = time at computation end
	2_iknd,&				!ITOL = tol setting per problem statement
	tolerance,&				!RTOL = user specified relative tolerance
	chp_ATOL,&				!ATOL = absolute tolerance vector specified in data_mod
	i_1,&					!ITASK = 1 - normal computation
	chp_istate,&			!ISTATE = 1 - first call
	i_1,&					!IOPT = 1 - optional input given
	chp_rwork,&				!RWORK Provided real work array
	chp_LRW,&				!LRW = len rwork
	chp_iwork,&				!IWORK = provded int work array plus optional max-step
	chp_LIW,&				!LIW = len iwork
	chp_JAC,&				!JAC = Pass dummy
	chp_MF)					!MF = specify type
END SUBROUTINE DLSODE_caller

END PROGRAM chp_main