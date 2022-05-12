subroutine build()
    USE helpers, only : parms_vanGI, intgrt_n

    integer*4 :: attribs
    character(*) :: wd, build_file, gok_file, bot_file, n_file, k_file, vanG_pars_file

    CALL getcwd(wd)
    input_path = trim(wd)//'/input/'
    output_path = trim(wd)//'/output/'

    build_file = trim(input_path)//"build.dat"
    gok_file = trim(input_path)//"gok.ip"
    bot_file = trim(input_path)//"bot.ip"
    n_file = trim(input_path)//"n.ip"
    k_file = trim(input_path)//"k.ip"
    vanG_pars_file = trim(input_path)//"vanG_pars.ip"

    open(unit=lu, file=log_file, status="replace")
    write(lu,*) "build initialized"

    open(tu, file=build_file)
	read(tu, *) attribs
    close(tu, status='keep')
	elems = attribs(1)
	ts = attribs(2)
	ts_size = attribs(3)

    allocate(gok(elems), bot(elems), n(elems), k(elems), chd(elems), p(elems,nts), et(elems,nts))
    allocate(gws(elems,nts), sws(elems,nts), sm(elems,nts), epv(elems,nts), gw_dis(elems,nts), sw_dis(elems,nts), sm_dis(elems,nts)&
        , Qin(elems,nts), Qout(elems,nts), Qdiff(elems,nts))

    open(tu, file=gok_file, form='unformatted')
	read(tu) gok
    close(tu, status='keep')
    open(tu, file=bot_file, form='unformatted')
	read(tu) bot
    close(tu, status='keep')
    open(tu, file=n_file, form='unformatted')
	read(tu) n
    close(tu, status='keep')
    open(tu, file=k_file, form='unformatted')
	read(tu) k
    close(tu, status='keep')
    open(tu, file=vanG_pars_file, form='unformatted')
	read(tu) vanG_pars
    close(tu, status='keep')
    intgrt_n = 750

    write(lu,*) "built"
end subroutine