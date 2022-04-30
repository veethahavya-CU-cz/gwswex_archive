function vanGI_simps(d, simpsons_n)
    real*8 :: vanGI_simps
    real*8, intent(inout) :: d
    integer, intent(in) :: simpsons_n
    d = -d/100
    vanGI_simps = simpsons(theta, d, dble(0), simpsons_n)*100
    contains
        function theta(h_c)
            real*8, intent(in) :: h_c
            real*8 :: theta, theta_r, theta_s, alpha, n, m
            theta_r = parms_vanGI(1)
            theta_s = parms_vanGI(2)
            alpha = parms_vanGI(3)
            n = parms_vanGI(4)
            m = (1-(1/n))
            theta = theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m)
        end function theta
end function vanGI_simps


function vanGI_rect(d, rect_n)
    real*8 :: vanGI_rect
    real*8, intent(inout) :: d
    integer, intent(in) :: rect_n
    d = -d/100
    vanGI_rect = simpsons(theta, d, dble(0), rect_n)*100
    contains
        function theta(h_c)
            real*8, intent(in) :: h_c
            real*8 :: theta, theta_r, theta_s, alpha, n, m
            theta_r = parms_vanGI(1)
            theta_s = parms_vanGI(2)
            alpha = parms_vanGI(3)
            n = parms_vanGI(4)
            m = (1-(1/n))
            theta = theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m)
        end function theta
end function vanGI_rect


function vanGI_trap(d, trap_n)
    real*8 :: vanGI_trap
    real*8, intent(inout) :: d
    integer, intent(in) :: trap_n
    d = -d/100
    vanGI_trap = simpsons(theta, d, dble(0), trap_n)*100
    contains
        function theta(h_c)
            real*8, intent(in) :: h_c
            real*8 :: theta, theta_r, theta_s, alpha, n, m
            theta_r = parms_vanGI(1)
            theta_s = parms_vanGI(2)
            alpha = parms_vanGI(3)
            n = parms_vanGI(4)
            m = (1-(1/n))
            theta = theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m)
        end function theta
end function vanGI_trap

function theta(h_c, ptr) bind(c)
    USE, intrinsic :: iso_c_binding
    real(c_double) :: h_c
    type(c_ptr), value :: ptr
    real(c_double) :: theta, theta_r, theta_s, alpha, n, m
    real(c_double), pointer :: null_ptr
    
    call c_f_pointer(ptr, null_ptr)

    theta_r = parms_vanGI(1)
    theta_s = parms_vanGI(2)
    alpha = parms_vanGI(3)
    n = parms_vanGI(4)
    m = (1-(1/n))

    theta = theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m)

end function theta

function vanGI_fgsl(d)
    USE fgsl
    USE, intrinsic :: iso_c_binding
    implicit none
    real*8 :: vanGI_fgsl
    real*8, intent(inout) :: d
    integer(fgsl_size_t), parameter :: nmax=1000
    real(fgsl_double), target :: ptr
    real(fgsl_double) :: result, error
    integer(fgsl_int) :: status
    type(c_ptr) :: null_ptr
    type(fgsl_function) :: f_obj
    type(fgsl_integration_workspace) :: wk

    ptr = 1.0D0
    null_ptr = c_loc(ptr)
    f_obj = fgsl_function_init(theta, ptr)
    wk = fgsl_integration_workspace_alloc(nmax)

    d = -d/100
    status = fgsl_integration_qags(f_obj, d, 0.0_fgsl_double, 0.0_fgsl_double, 1.0e-7_fgsl_double, nmax, wk, result, error)
    vanGI_fgsl = result

    call fgsl_function_free(f_obj)
    call fgsl_integration_workspace_free(wk)

end function vanGI_fgsl