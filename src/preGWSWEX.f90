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



MODULE Mlogger

	USE timing, only: tglobal_start, timefetch

	TYPE Clogger
		INTEGER(1) :: unit, level, info, moreinfo, trace, debug, warn, error, fatal
		CHARACTER(len=255) :: fname
		CONTAINS
			PROCEDURE :: init
			PROCEDURE :: log_real
			PROCEDURE :: log_int
			PROCEDURE :: log_str
			GENERIC :: log => log_real, log_int, log_str
	END TYPE Clogger

	CONTAINS
		INCLUDE 'logger.f90'
END MODULE Mlogger



MODULE sm_helpers
	USE YAMLInterface
	USE YAMLRead

	IMPLICIT NONE
	REAL(8) :: alpha, n, m, theta_r, theta_s
	
	TYPE CvanG
		REAL(8) :: alpha, n, m, theta_r, theta_s
		CONTAINS
			PROCEDURE, PASS :: init
			PROCEDURE, PASS :: integrate
	END TYPE CvanG

	CONTAINS
		SUBROUTINE init(self, ymodel_vanG_args)
			IMPLICIT NONE
			CLASS(CvanG) :: self
			TYPE(YAMLMap), INTENT(INOUT) :: ymodel_vanG_args
			INTEGER :: ires

			self% alpha = ymodel_vanG_args% value_double('alpha', ires)
			write(*,*) 'alpha = ', self% alpha
			self% n = ymodel_vanG_args% value_double('n', ires)
			write(*,*) 'n = ', self% n
			self% m = (1-(1/self% n))
			self% theta_r = ymodel_vanG_args% value_double('theta_r', ires)
			write(*,*) 'theta_r = ', self% theta_r
			self% theta_s = ymodel_vanG_args% value_double('theta_s', ires)
			write(*,*) 'theta_s = ', self% theta_s
			
			alpha = self% alpha
			n = self% n
			m = self% m
			theta_r = self% theta_r
			theta_s = self% theta_s
			CALL ymodel_vanG_args% destroy()
		END SUBROUTINE init

		INCLUDE 'kUS.f90'
		INCLUDE 'vanGI.f90'

END MODULE sm_helpers



MODULE core

	USE sm_helpers, only: CvanG
	USE timing
	USE Mlogger, only: Clogger

	IMPLICIT NONE

	INTEGER(4)  :: elems, nts, dt
	LOGICAL, ALLOCATABLE  :: chd(:)
	REAL(8), ALLOCATABLE :: gok(:), bot(:), n(:), k(:), gw_sm_interconnectivity(:), macropore_inf_degree(:), p(:,:), et(:,:)
	REAL(8), ALLOCATABLE :: gws(:,:), sws(:,:), sm(:,:), epv(:,:), gw_dis(:,:), sw_dis(:,:), sm_dis(:,:), &
		Qin(:,:), Qout(:,:), QdIFf(:,:)
	
	CHARACTER(255) :: input_path, output_path, strbuffer
	INTEGER, PARAMETER  :: lu=42, tu=99

	TYPE(Clogger) :: logger
	TYPE(CvanG) :: vanG

	CONTAINS
		INCLUDE 'build.f90'
		INCLUDE 'init.f90'
		INCLUDE 'run.f90'
END MODULE core