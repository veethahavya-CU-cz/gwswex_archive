SUBROUTINE build()

	INTEGER*4 :: attribs(4)
	CHARACTER(255) :: wd, build_file, gok_file, bot_file, n_file, k_file, macropore_inf_degree_file, vanG_pars_file

	tglobal_start = timefetch()
	tlocal_start = timefetch()

	CALL getcwd(wd)
	input_path = TRIM(wd)//'/input/'
	output_path = TRIM(wd)//'/output/'

	build_file = TRIM(input_path)//"build.dat"
	gok_file = TRIM(input_path)//"gok.ip"
	bot_file = TRIM(input_path)//"bot.ip"
	n_file = TRIM(input_path)//"n.ip"
	k_file = TRIM(input_path)//"k.ip"
	macropore_inf_degree_file = TRIM(input_path)//"macropore_inf_degree.ip"
	vanG_pars_file = TRIM(input_path)//"vanG_pars.ip"

	OPEN(tu, file=build_file, action='READ')
	READ(tu, *) attribs
	elems = attribs(1)
	nts = attribs(2)
	dt = attribs(3)
	logger%level = INT(attribs(4), kind=1)
	CLOSE(tu, status='keep')

	logger%fname = TRIM(output_path)//'GWSWEX.log'
	CALL logger%init()

	CALL logger%log(logger%info, "Program called")
	CALL logger%log(logger%info, "building model")

	allocate(gok(elems), bot(elems), n(elems), k(elems), chd(elems), gw_sm_interconnectivity(elems), macropore_inf_degree(elems), &
		p(elems,nts), et(elems,nts))
	allocate(gws(elems,nts), sws(elems,nts), sm(elems,nts), epv(elems,nts), gw_dis(elems,nts), sw_dis(elems,nts), sm_dis(elems,nts), &
		Qin(elems,nts), Qout(elems,nts), QdIFf(elems,nts))

	OPEN(tu, file=gok_file, form='unformatted', action='READ')
	READ(tu) gok
	CLOSE(tu, status='keep')
	OPEN(tu, file=bot_file, form='unformatted', action='READ')
	READ(tu) bot
	CLOSE(tu, status='keep')
	OPEN(tu, file=n_file, form='unformatted', action='READ')
	READ(tu) n
	CLOSE(tu, status='keep')
	OPEN(tu, file=k_file, form='unformatted', action='READ')
	READ(tu) k
	CLOSE(tu, status='keep')
	OPEN(tu, file=macropore_inf_degree_file, form='unformatted', action='READ')
	READ(tu) macropore_inf_degree
	CLOSE(tu, status='keep')
	OPEN(tu, file=vanG_pars_file, form='unformatted', action='READ')
	READ(tu) vanG_pars
	CLOSE(tu, status='keep')

	tlocal_end = timefetch()

	write(strbuffer,*) (tlocal_end-tlocal_start)
	strbuffer = "model built in "//TRIM(strbuffer)//" s"//achar(10)
	CALL logger%log(logger%info, TRIM(strbuffer))
END SUBROUTINE