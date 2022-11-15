MODULE timing
	USE OMP_LIB
	REAL(8) :: tglobal_start, tglobal_end, tlocal_start, tlocal_end, ttemp_mark1, ttemp_mark2
	REAL(4) :: telapsed, timenow
	CONTAINS
		FUNCTION timefetch()
			IMPLICIT NONE
			REAL(8) :: timefetch
			timefetch = omp_get_wtime()
		END FUNCTION timefetch
END MODULE timing


MODULE logger_mod
	USE timing, only: tglobal_start, timefetch
	TYPE :: logger_type
		INTEGER(1) :: lu, level, info, moreinfo, trace, debug, warn, error, fatal
		CHARACTER(len=255) :: fname
		CONTAINS
			PROCEDURE :: init
			PROCEDURE :: log_real
			PROCEDURE :: log_int
			PROCEDURE :: log_str
			GENERIC :: log => log_real, log_int, log_str
	END TYPE

	CONTAINS
	INCLUDE 'logger.f90'
END MODULE logger_mod


MODULE helpers
	IMPLICIT NONE

	REAL*8 :: vanG_pars(4)

	CONTAINS
		INCLUDE 'kUS.f90'
		INCLUDE 'vanGI.f90'
END MODULE helpers


MODULE core
	USE helpers, only: vanG_pars
	USE timing
	USE logger_mod, only: logger_type

	IMPLICIT NONE

	INTEGER*4  :: elems, nts, dt
	LOGICAL, ALLOCATABLE  :: chd(:)
	REAL*8, ALLOCATABLE :: gok(:), bot(:), n(:), k(:), gw_sm_interconnectivity(:), macropore_inf_degree(:), p(:,:), et(:,:)
	REAL*8, ALLOCATABLE :: gws(:,:), sws(:,:), sm(:,:), epv(:,:), gw_dis(:,:), sw_dis(:,:), sm_dis(:,:), &
		Qin(:,:), Qout(:,:), QdIFf(:,:)
	
	CHARACTER(255) :: input_path, output_path, strbuffer
	INTEGER, PARAMETER  :: lu=42, tu=99

	TYPE(logger_type) :: logger

	CONTAINS
		INCLUDE 'build.f90'
		INCLUDE 'init.f90'
		INCLUDE 'run.f90'
END MODULE core