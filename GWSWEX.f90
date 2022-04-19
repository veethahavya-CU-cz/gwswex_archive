module helpers
    implicit none
contains
    function kSM(s, ks)
        real*8, intent(in) :: s, ks
        real*8 :: kSM, n
        n = 2.5
        kSM = ks*s*(1-(1-(s)**1/(1-(1/n)))**(1-(1/n)))**2
    end function kSM

    function kGW(s, ks)
        real*8, intent(in) :: s, ks
        real*8 :: kGW, n
        n = 2.5
        kGW = ks*s*(1-(1-(s)**1/(1-(1/n)))**(1-(1/n)))**2
    end function kGW
end module helpers

module gwswex
    USE OMP_LIB
    implicit none
    integer  :: elems, nts, dt
    logical, allocatable  :: chd(:)
    real*8, allocatable :: gok(:), bot(:), n(:), k(:), p(:,:), et(:,:)
contains

    subroutine build(el, ts, ts_size, gok_l, bot_l, n_l, k_l)
        integer, intent(in) :: el, ts, ts_size
        real*8, intent(in) :: n_l, bot_l(:), gok_l(:), k_l(:)
        elems = el
        nts = ts
        dt = ts_size
        allocate(gok(elems), bot(elems), n(elems), chd(elems), p(elems,nts), et(elems,nts))
        gok = gok_l
        bot = bot_l
        n = n_l
        k = k_l
        open(unit=42, file="fort.log", status="replace")
        write(42,*) "built"
    end subroutine

    subroutine init(chd_l, p_l, et_l)
        real*8, intent(in) :: p_l(:,:), et_l(:,:)
        logical, intent(in) :: chd_l(:)
        p = p_l
        et = et_l
        chd = chd_l
        open(unit=42, file="fort.log", status="old")
        write(42,*) "initialised"
    end subroutine

    subroutine run(vanGI, gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff)
        USE helpers
        external :: vanGI
        real*8 :: vanGI
        !f2py real*8, intent(in):: d
        !f2py real*8, intent(out) :: eq
        !f2py eq = vanGI(d)
        real*8, intent(inout) :: gws(:,:), sws(:,:), sm(:,:), epv(:,:), gw_dis(:,:), sw_dis(:,:), sm_dis(:,:), &
            Qin(:,:), Qout(:,:), Qdiff(:,:)
        real*8 :: L, sw_et_deficit, excess_gw_vol, sm_eq, k_inf, inf, excess_p, inf_deficit, sw_inf, &
            k_inf_gw, gw_inf, et_deficit, sw_et
        integer :: e, t
        open(unit=42, file="fort.log", status="old")
        write(42,*) "run entered"

        do t = 2, nts
            write(42,*) "outer loop entered. ts ", t-1
            do e = 1, elems
                write(42,*) "inner loop entered. elem", e
                write(42,*) "gok", gok(e)
                write(42,*) "bot", bot(e)
                if(.NOT. chd(e)) then
                    L = gok(e) - gws(e,t-1) !prev. GW depth
                    if(L<0 .OR. L==0) then !NO UZ case
                        write(42,*) "noUZ entered"
                        !excess GW correction
                        write(42,*) "gws is ", gws(e,t-1)
                        write(42,*) "sws is ", sws(e,t-1)
                        write(42,*) "sm is ", sm(e,t-1)
                        excess_gw_vol = -L*n(e) + sm(e,t-1)
                        gws(e,t) = gok(e)
                        sm(e,t) = 0
                        epv(e,t) = 0
                        sws(e,t) = sws(e,t-1) + excess_gw_vol + p(e,t)*dt
                        write(42,*) "excess_gw_vol ", excess_gw_vol
                        write(42,*) "gws after +p ", gws(e,t)
                        write(42,*) "sws after +p", sws(e,t)
                        write(42,*) "sm after +p", sm(e,t)
                        !ET extraction
                        if (sws(e,t)>et(e,t)*dt) then
                            sws(e,t) = sws(e,t) - et(e,t)*dt
                        else
                            sw_et_deficit = et(e,t)*dt - sws(e,t)
                            sws(e,t) = 0
                            gws(e,t) = gws(e,t) - (sw_et_deficit/n(e))
                            epv(e,t) = (gok(e) - gws(e,t))*n(e)
                        end if
                        write(42,*) "gws after -et ", gws(e,t)
                        write(42,*) "sws after -et", sws(e,t)
                        write(42,*) "sm after -et", sm(e,t)
                        !calc storage discharges
                        gw_dis(e,t) = (gws(e,t) - gws(e,t-1))*n(e)
                        sm_dis(e,t) = (sm(e,t)) - sm(e,t-1)
                        sw_dis(e,t) = sws(e,t) - sws(e,t-1)
                        Qin(e,t) = p(e,t)*dt - et(e,t)*dt
                        Qout(e,t) = gw_dis(e,t) + sw_dis(e,t) + sm_dis(e,t)
                        sw_et_deficit = 0
                    else
                        write(42,*) "UZ entered"
                        !P dist and SW push
                        write(42,*) "L is", L
                        write(42,*) "P is", p(e,t)*dt
                        write(42,*) "sm is", sm(e,t-1)
                        write(42,*) "epv is", epv(e,t-1)
                        write(42,*) "sm/epv", sm(e,t-1)/epv(e,t-1)
                        k_inf = kSM(min(sm(e,t-1)/epv(e,t-1), 1.0)*n(e), k(e)) !calc K from wetness at the begining of this dt i.e. end of last dt
                        write(42,*) "got k", k_inf
                        inf = max(min(k_inf*dt, p(e,t)*dt, epv(e,t-1)-sm(e,t-1)), 0.0)
                        write(42,*) "inf is", inf
                        excess_p = p(e,t)*dt - inf
                        write(42,*) "excess_p", excess_p
                        inf_deficit = k_inf*dt - inf
                        write(42,*) "inf_deficit", inf_deficit
                        write(42,*) "sws is", sws(e,t-1)
                        sw_inf = max(min(inf_deficit, sws(e,t-1), epv(e,t-1)-sm(e,t-1)), 0.0)
                        sws(e,t) = sws(e,t-1) - sw_inf + excess_p
                        write(42,*) "sw_inf", sw_inf
                        write(42,*) "sws calcd", sws(e,t)
                        sm(e,t) = sm(e,t-1) + inf + sw_inf
                        write(42,*) "sm calcd", sm(e,t)
                        sm_eq = vanGI(-L) !!!consider doing GW push and gw-sm bal after ET extraction
                        write(42,*) "gws is ", gws(e,t-1)
                        write(42,*) "vanGI called. sm_eq is ", sm_eq
                        k_inf_gw = kGW(min(sm(e,t-1)/epv(e,t-1), 1.0)*n(e), k(e)) !calc K from current wetness (after P and SW inf)
                        gw_inf = min(sm(e,t)-sm_eq, k_inf_gw*dt) !if sm<sm_eq, gw_inf is -ve ...
                        write(42,*) "k_inf_gw is", k_inf_gw
                        write(42,*) "gw_inf is", gw_inf
                        sm(e,t) = sm(e,t) - gw_inf !... deficit sm gets added to sm from gw
                        write(42,*) "sm recalcd ", sm(e,t)
                        gws(e,t) = gws(e,t-1) + gw_inf/n(e) !... and subtracted from gw
                        write(42,*) "gws calcd", gws(e,t)
                        if(gws(e,t)>gok(e)) then
                            excess_gw_vol = (gws(e,t)-gok(e))*n(e) + sm(e,t)
                            gws(e,t) = gok(e)
                            sm(e,t) = 0
                            sws(e,t) = sws(e,t) + excess_gw_vol
                        end if
                        write(42,*) "gws recalcd", gws(e,t)

                        !ET removal and SM-GW rebalance
                        write(42,*) "ET is", et(e,t)*dt
                        sw_et = min(sws(e,t), et(e,t)*dt)
                        sws(e,t) = sws(e,t) - sw_et
                        et_deficit = et(e,t)*dt - sw_et
                        write(42,*) "sw et removed", sw_et
                        sm(e,t) = sm(e,t) - et_deficit
                        write(42,*) "sm et removed", et_deficit
                        write(42,*) "sm is", sm(e,t)
                        sm_eq = vanGI(-(gok(e) - gws(e,t))) !!!gw-sm balancing: consider adding a convergence criteria here
                        write(42,*) "new sm_eq", sm_eq
                        k_inf_gw = kGW(min(sm(e,t-1)/epv(e,t-1), 1.0)*n(e), k(e))*dt - max(gw_inf, 0.00) !subtract k_inf_gw already utilized and allow freely capilary rise beyond k_inf_gw
                        write(42,*) "k_inf_gw remaining", k_inf_gw
                        gw_inf = min(sm(e,t)-sm_eq, k_inf_gw*dt)
                        write(42,*) "addnl gw_inf", gw_inf
                        sm(e,t) = sm(e,t) - gw_inf
                        gws(e,t) = gws(e,t) + gw_inf/n(e)
                        write(42,*) "sm-gw balanced", sm(e,t), gws(e,t)

                        epv(e,t) = (gok(e) - gws(e,t))*n(e)
                        gw_dis(e,t) = (gws(e,t) - gws(e,t-1))*n(e)
                        sw_dis(e,t) = sws(e,t) - (sws(e,t-1))
                        sm_dis(e,t) = sm(e,t) - sm(e,t-1)
                        Qin(e,t) = p(e,t)*dt - et(e,t)*dt
                        Qout(e,t) = gw_dis(e,t) + sw_dis(e,t) + sm_dis(e,t)
                    end if
                else
                    excess_gw_vol = sm(e,t-1)
                    gws(e,t) = gws(e,t-1)
                    sm(e,t) = 0
                    epv(e,t) = 0
                    sws(e,t) = sws(e,t-1) + p(e,t)*dt - et(e,t)*dt + excess_gw_vol		
                    gw_dis(e,t) = 0
                    sw_dis(e,t) = sws(e,t) - sws(e,t-1)
                    sm_dis(e,t) = 0
                    Qin(e,t) = p(e,t)*dt - et(e,t)*dt
                    Qout(e,t) = gw_dis(e,t) + sw_dis(e,t) + sm_dis(e,t)
                end if
            end do
        end do
        Qdiff = Qin - Qout
    end subroutine
end module gwswex