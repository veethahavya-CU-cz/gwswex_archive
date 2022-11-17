FUNCTION theta_c(h_c, ptr_c) bind(c)
	USE, INTRINSIC :: iso_c_binding
	USE fgsl
	
	REAL(c_DOuble), value :: h_c
	TYPE(c_ptr), value :: ptr_c
	REAL(c_DOuble), pointer :: ptr_f
	REAL(c_DOuble) :: theta_c
	REAL(8) :: theta_r, theta_s, alpha, n, m

	CALL c_f_pointer(ptr_c, ptr_f)

	theta_r = vanG_pars(1)
	theta_s = vanG_pars(2)
	alpha = vanG_pars(3)
	n = vanG_pars(4)
	m = (1-(1/n))
	theta_c = theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m)
END FUNCTION theta_c

FUNCTION vanGI_fgsl(d)
	USE fgsl
	USE, INTRINSIC :: iso_c_binding
	IMPLICIT NONE
	REAL(8) :: vanGI_fgsl
	REAL(8), INTENT(INOUT) :: d
	INTEGER(fgsl_size_t), PARAMETER :: nmax=1000
	REAL(fgsl_DOuble), target :: ptr
	REAL(fgsl_DOuble) :: result, error
	INTEGER(fgsl_int) :: status
	TYPE(c_ptr) :: null_ptr
	TYPE(fgsl_function) :: f_obj
	TYPE(fgsl_integration_workspace) :: wk

	ptr = 1.0D0
	null_ptr = c_loc(ptr)
	f_obj = fgsl_function_init(theta_c, null_ptr)
	wk = fgsl_integration_workspace_alloc(nmax)

	status = fgsl_integration_qags(f_obj, -d, 0.0_fgsl_DOuble, 0.0_fgsl_DOuble, 1.0e-7_fgsl_DOuble, nmax, wk, result, error)
	vanGI_fgsl = result

	CALL fgsl_function_free(f_obj)
	CALL fgsl_integration_workspace_free(wk)

END FUNCTION vanGI_fgsl