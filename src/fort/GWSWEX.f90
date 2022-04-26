module helpers

    USE OMP_LIB
    implicit none

    real*8 :: parms_vanGI(4)

    contains
        include "kSM.f90"
        include "kGW.f90"
        include "integrators.f90"
        include "vanGI.f90"

end module helpers



module gwswex

    USE OMP_LIB
    USE stdlib_logger

    implicit none

    integer  :: elems, nts, dt
    logical, allocatable  :: chd(:)
    real*8, allocatable :: gok(:), bot(:), n(:), k(:), p(:,:), et(:,:)
    real*8 :: vanG_pars(4)

    contains
        include "build.f90"
        include "init.f90"
        include "run.f90"
        include "run_f.f90"

end module gwswex