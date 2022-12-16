SUBROUTINE build(fyaml_path, gok_l, bot_l, n_l, k_l, macropore_inf_degree_l)
	USE YAMLInterface
	USE YAMLRead

	CHARACTER(*), INTENT(IN) :: fyaml_path
	TYPE(YAMLHandler) :: fyaml  ! be sure to close
	TYPE(YAMLMap) :: ymodel_args, yutil_args, ymodel_domain_args, ymodel_vanG_args, yutil_logger_args
	INTEGER :: ires

	REAL(8), INTENT(IN) :: gok_l(:), bot_l(:), n_l(:), k_l(:), macropore_inf_degree_l(:)
	CHARACTER(255) :: wd

	tglobal_start = timefetch()
	tlocal_start = timefetch()

	fyaml = yaml_open_file(fyaml_path)

	ymodel_args = yaml_start_from_map(fyaml, 'model')
	ymodel_domain_args = ymodel_args%value_map('domain')
	ymodel_vanG_args = ymodel_args%value_map('vanG')
	elems = ymodel_domain_args%value_int("nelements", ires)
	nts = ymodel_domain_args%value_int("nts", ires)
	dt = ymodel_domain_args%value_int("dt", ires)
	CALL vanG% init(ymodel_vanG_args)
	call ymodel_args%destroy()
	call ymodel_domain_args%destroy()

	yutil_args = yaml_start_from_map(fyaml, 'utils')
	yutil_logger_args = yutil_args%value_map('logger')

	logger% level = INT(yutil_logger_args%value_int("level", ires), kind=1)
	call yutil_args%destroy()
	call yutil_logger_args%destroy()

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
	
	tlocal_end = timefetch()

	write(strbuffer,*) (tlocal_end-tlocal_start)
	strbuffer = "model built in "//TRIM(strbuffer)//" s"//achar(10)
	CALL logger%log(logger%info, TRIM(strbuffer))
	flush(logger%unit)
END SUBROUTINE