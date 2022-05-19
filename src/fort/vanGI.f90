function theta_c(h_c, ptr_c) bind(c)
	USE, intrinsic :: iso_c_binding
	USE fgsl
	
	real(c_double), value :: h_c
	type(c_ptr), value :: ptr_c
	real(c_double), pointer :: ptr_f
	real(c_double) :: theta_c
	real*8 :: theta_r, theta_s, alpha, n, m

	call c_f_pointer(ptr_c, ptr_f)

	theta_r = vanG_pars(1)
	theta_s = vanG_pars(2)
	alpha = vanG_pars(3)
	n = vanG_pars(4)
	m = (1-(1/n))
	theta_c = theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m)
end function theta_c

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
	f_obj = fgsl_function_init(theta_c, null_ptr)
	wk = fgsl_integration_workspace_alloc(nmax)

	status = fgsl_integration_qags(f_obj, -d, 0.0_fgsl_double, 0.0_fgsl_double, 1.0e-7_fgsl_double, nmax, wk, result, error)
	vanGI_fgsl = result

	call fgsl_function_free(f_obj)
	call fgsl_integration_workspace_free(wk)

end function vanGI_fgsl