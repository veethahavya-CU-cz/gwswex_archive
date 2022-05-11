subroutine build(el, ts, ts_size, gok_l, bot_l, n_l, k_l, vanG_pars_l, intgrt)
    USE helpers, only : parms_vanGI, intgrt_n
    integer, intent(in) :: el, ts, ts_size
    integer, optional :: intgrt
    real*8, intent(in) :: gok_l(:), bot_l(:), n_l(:), k_l(:), vanG_pars_l(4)
    elems = el
    nts = ts
    dt = ts_size
    allocate(gok(elems), bot(elems), n(elems), k(elems), chd(elems), p(elems,nts), et(elems,nts))
    gok = gok_l
    bot = bot_l
    n = n_l
    k = k_l
    vanG_pars = vanG_pars_l
    parms_vanGI = vanG_pars_l
    if (present(intgrt)) then
        intgrt_n = intgrt
    else
        intgrt_n = 750
    end if
    open(unit=42, file="fort.log", status="replace")
    write(42,*) "built"
end subroutine