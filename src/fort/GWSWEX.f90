module helpers
	USE OMP_LIB
	implicit none

	real*8 :: vanG_pars(4)
	integer*2 :: intgrt_n

	contains
		include "kSM.f90"
		include "kGW.f90"
		include "vanGI.f90"
end module helpers


module core
	USE OMP_LIB
	USE helpers, only: vanG_pars

	implicit none

	integer*4  :: elems, nts, dt
	logical, allocatable  :: chd(:)
	real*8, allocatable :: gok(:), bot(:), n(:), k(:), p(:,:), et(:,:)
	real*8, allocatable :: gws(:,:), sws(:,:), sm(:,:), epv(:,:), gw_dis(:,:), sw_dis(:,:), sm_dis(:,:)&
		, Qin(:,:), Qout(:,:), Qdiff(:,:)
	
	character(255) :: input_path, output_path, log_file
	integer, parameter  :: lu=42, tu=99

	contains
		include "build.f90"
		include "init.f90"
		include "run.f90"
end module core


program GWSWEX
	USE core, only: build, init, run

	CALL build()
	CALL init()
	CALL run()
end program GWSWEX