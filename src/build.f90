SUBROUTINE build(fyaml_path, gok_l, bot_l, n_l, k_l, macropore_inf_degree_l, vanG_pars_l)
	USE YAMLInterface
	USE YAMLRead
	!USE iso_c_binding, only: c_char

	!CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: fyaml_path
	CHARACTER(*), INTENT(IN) :: fyaml_path
	TYPE(YAMLHandler) :: fyaml  ! be sure to close
	TYPE(YAMLMap) :: ymodel_args, yutil_args, ymodel_domain_args, yutil_logger_args
	INTEGER :: ires

	REAL(8), INTENT(IN) :: gok_l(:), bot_l(:), n_l(:), k_l(:), macropore_inf_degree_l(:)
	REAL(8) :: vanG_pars_l(4)
	CHARACTER(255) :: wd

	tglobal_start = timefetch()
	tlocal_start = timefetch()

	fyaml = yaml_open_file(fyaml_path)
	ymodel_args = yaml_start_from_map(fyaml, 'model')
	! yutil_args = yaml_start_from_map(fyaml, 'utils')

	ymodel_domain_args = ymodel_args%value_map('domain')
	! yutil_logger_args = yutil_args%value_map('logger')

	elems = ymodel_domain_args%value_int("nelements", ires)
	nts = ymodel_domain_args%value_int("nts", ires)
	dt = ymodel_domain_args%value_int("dt", ires)
	! logger%level = INT(yutil_logger_args%value_int("dt", ires), kind=1)
	logger%level = INT(1, kind=1)

	call ymodel_args%destroy()
	! call yutil_args%destroy()
	call ymodel_domain_args%destroy()
	! call yutil_logger_args%destroy()
	call yaml_close_file(fyaml)

	CALL getcwd(wd)
	input_path = TRIM(wd)//'/input/'
	output_path = TRIM(wd)//'/output/'

	logger%fname = 'GWSWEX.log'
	logger%unit = lu
	CALL logger%init()

	CALL logger%log(logger%info, "Program called")
	CALL logger%log(logger%info, "building model")

	ALLOCATE(gok(elems), bot(elems), n(elems), k(elems), chd(elems), gw_sm_interconnectivity(elems), macropore_inf_degree(elems), &
		p(elems,nts), et(elems,nts))
	ALLOCATE(gws(elems,nts), sws(elems,nts), sm(elems,nts), epv(elems,nts), gw_dis(elems,nts), sw_dis(elems,nts), sm_dis(elems,nts), &
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