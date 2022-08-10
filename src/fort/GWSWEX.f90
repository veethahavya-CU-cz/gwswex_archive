MODULE helpers
	USE OMP_LIB
	IMPLICIT NONE

	REAL*8 :: vanG_pars(4)
	INTEGER*2 :: intgrt_n

	contains
		include "kUS.f90"
		include "vanGI.f90"
END MODULE helpers


MODULE core
	USE OMP_LIB
	USE helpers, only: vanG_pars

	IMPLICIT NONE

	INTEGER*4  :: elems, nts, dt
	logical, allocatable  :: chd(:)
	REAL*8, allocatable :: gok(:), bot(:), n(:), k(:), gw_sw_interconnectivity(:), macropore_inf_degree(:), p(:,:), et(:,:)
	REAL*8, allocatable :: gws(:,:), sws(:,:), sm(:,:), epv(:,:), gw_dis(:,:), sw_dis(:,:), sm_dis(:,:), &
		Qin(:,:), Qout(:,:), QdIFf(:,:)
	
	CHARACTER(255) :: input_path, output_path, log_file
	INTEGER, parameter  :: lu=42, tu=99

	contains
		include "build.f90"
		include "init.f90"
		include "run.f90"
END MODULE core


program GWSWEX
	USE core, only: build, init, run

	CALL build()
	CALL init()
	CALL run()
END program GWSWEX