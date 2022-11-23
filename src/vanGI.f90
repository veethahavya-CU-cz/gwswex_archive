FUNCTION theta_c(h_c, ptr_c) bind(c)
	USE, INTRINSIC :: iso_c_binding
	USE fgsl
	
	REAL(c_DOuble), value :: h_c
	TYPE(c_ptr), value :: ptr_c
	REAL(c_DOuble), pointer :: ptr_f
	REAL(c_DOuble) :: theta_c

	CALL c_f_pointer(ptr_c, ptr_f)

	theta_c = theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m)
END FUNCTION theta_c


FUNCTION integrate(self, d)
	USE fgsl
	USE, INTRINSIC :: iso_c_binding
	IMPLICIT NONE
	CLASS(CvanG) :: self
	REAL(8) :: integrate
	REAL(8), INTENT(INOUT) :: d
	INTEGER(fgsl_size_t), PARAMETER :: nmax=1000
	REAL(fgsl_DOuble), target :: ptr
	REAL(fgsl_DOuble) :: result, error
	INTEGER(fgsl_int) :: status
	TYPE(c_ptr) :: cptr
	TYPE(fgsl_function) :: f_obj
	TYPE(fgsl_integration_workspace) :: wk

	ptr = 1.0D0
	cptr = c_loc(ptr)
	f_obj = fgsl_function_init(theta_c, cptr)
	wk = fgsl_integration_workspace_alloc(nmax)

	status = fgsl_integration_qags(f_obj, -d, 0.0_fgsl_DOuble, 0.0_fgsl_DOuble, 1.0e-7_fgsl_DOuble, nmax, wk, result, error)
	integrate = result

	CALL fgsl_function_free(f_obj)
	CALL fgsl_integration_workspace_free(wk)

END FUNCTION integrate