SUBROUTINE build(attribs, gok_l, bot_l, n_l, k_l, macropore_inf_degree_l, vanG_pars_l)

	INTEGER(4), INTENT(IN) :: attribs(4)
	REAL(8), INTENT(IN) :: gok_l(:), bot_l(:), n_l(:), k_l(:), macropore_inf_degree_l(:)
	REAL(8) :: vanG_pars_l(4)
	CHARACTER(255) :: wd

	tglobal_start = timefetch()
	tlocal_start = timefetch()

	elems = attribs(1)
	nts = attribs(2)
	dt = attribs(3)
	logger%level = INT(attribs(4), kind=1)

	CALL getcwd(wd)
	input_path = TRIM(wd)//'/input/'
	output_path = TRIM(wd)//'/output/'

	logger%fname = 'GWSWEX.log'
	logger%unit = lu
	CALL logger%init()
	write(*,*) 'OGu', logger%unit

	CALL logger%log(logger%info, "Program called")
	CALL logger%log(logger%info, "building model")

	allocate(gok(elems), bot(elems), n(elems), k(elems), chd(elems), gw_sm_interconnectivity(elems), macropore_inf_degree(elems), &
		p(elems,nts), et(elems,nts))
	allocate(gws(elems,nts), sws(elems,nts), sm(elems,nts), epv(elems,nts), gw_dis(elems,nts), sw_dis(elems,nts), sm_dis(elems,nts), &
		Qin(elems,nts), Qout(elems,nts), QdIFf(elems,nts))

	gok = gok_l
	bot = bot_l
	n = n_l
	k = k_l
	macropore_inf_degree = macropore_inf_degree_l
	vanG_pars = vanG_pars_l
	
	tlocal_end = timefetch()

	write(strbuffer,*) (tlocal_end-tlocal_start)
	strbuffer = "model built in "//TRIM(strbuffer)//" s"//achar(10)
	CALL logger%log(logger%info, TRIM(strbuffer))
	flush(logger%unit)
END SUBROUTINE