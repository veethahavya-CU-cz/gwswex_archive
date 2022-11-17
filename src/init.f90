SUBROUTINE init(chd_l, p_l, et_l)
	LOGICAL, INTENT(IN)  :: chd_l(:)
	REAL(8), INTENT(IN) :: p_l(:,:), et_l(:,:)

	tlocal_start = timefetch()
	CALL logger%log(logger%info, "initialising")

	chd = chd_l
	p = p_l
    et = et_l

	tlocal_end = timefetch()
	
	write(strbuffer,*) (tlocal_end-tlocal_start)
	strbuffer = "initialised in "//TRIM(strbuffer)//" s"//achar(10)
	CALL logger%log(logger%info, strbuffer)
	flush(logger%unit)
END SUBROUTINE