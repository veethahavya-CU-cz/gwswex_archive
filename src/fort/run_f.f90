subroutine run_f(gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff)
!f2py threadsafe
    USE helpers
    real*8, intent(inout) :: gws(:,:), sws(:,:), sm(:,:), epv(:,:), gw_dis(:,:), sw_dis(:,:), sm_dis(:,:), &
        Qin(:,:), Qout(:,:), Qdiff(:,:)
    real*8 :: L, sw_et_deficit, excess_gw_vol, sm_eq, k_inf, inf, excess_p, inf_deficit, sw_inf, &
        k_inf_gw, inf_gw, et_deficit, sw_et, start, finish
    integer :: e, t
    open(unit=42, file="fort.log", status="old")
    write(42,*) "run entered"

    do t = 2, nts
        write(42,*) "outer loop entered. ts ", t-1
        !$OMP PARALLEL DO SHARED(gws, sws, sm, epv) &
        !$OMP PRIVATE(L, sw_et_deficit, excess_gw_vol, sm_eq, k_inf, inf, excess_p, inf_deficit, sw_inf, k_inf_gw, inf_gw) &
        !$OMP PRIVATE(et_deficit, sw_et)
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
                    k_inf = kSM(min(sm(e,t-1)/epv(e,t-1), 1.0)*n(e), k(e), vanG_pars) !calc K from wetness at the begining of this dt i.e. end of last dt
                    write(42,*) "got k", k_inf
                    inf = min(k_inf*dt, p(e,t)*dt)
                    write(42,*) "inf aka p_sm is ", inf
                    excess_p = p(e,t)*dt - inf
                    write(42,*) "excess p aka p_sw is ", excess_p
                    write(42,*) "sws is", sws(e,t-1)
                    write(42,*) "ET is", et(e,t)*dt
                    sw_et = min(sws(e,t-1)+excess_p, et(e,t)*dt)
                    inf_deficit = k_inf*dt - inf
                    write(42,*) "inf_deficit", inf_deficit
                    sw_inf = min(inf_deficit, sws(e,t-1)+excess_p-sw_et)
                    write(42,*) "sw_inf", sw_inf
                    sws(e,t) = sws(e,t-1) - sw_inf + excess_p - sw_et
                    et_deficit = et(e,t)*dt - sw_et
                    if(gws(e,t-1) <= bot(e)) then
                        et_deficit = 0
                    end if
                    write(42,*) "sw et removed", sw_et
                    write(42,*) "sws calcd", sws(e,t)
                    sm(e,t) = sm(e,t-1) + inf + sw_inf - et_deficit
                    write(42,*) "sm et removed", et_deficit
                    write(42,*) "sm calcd", sm(e,t)
                    call cpu_time(start)
                    sm_eq = vanGI_simps(L, 750)
                    call cpu_time(finish)
                    write(42,*) "vanGI_simps time ", finish-start
                    write(42,*) "gws is ", gws(e,t-1)
                    write(42,*) "vanGI_simps called. sm_eq is ", sm_eq
                    k_inf_gw = kGW(min(sm(e,t)/epv(e,t-1), 1.0)*n(e), k(e), vanG_pars) !calc K from current wetness (after P and SW inf)
                    inf_gw = min(sm(e,t)-sm_eq, k_inf_gw*dt) !if sm<sm_eq, inf_gw is -ve ...
                    if(gws(e,t-1) + inf_gw/n(e) < bot(e)) then
                        inf_gw = - min(abs((gws(e,t-1) - bot(e)))*n(e), abs(k_inf_gw*dt))
                    end if
                    write(42,*) "k_inf_gw is", k_inf_gw
                    write(42,*) "inf_gw is", inf_gw
                    sm(e,t) = sm(e,t) - inf_gw !... deficit sm gets added to sm from gw
                    write(42,*) "sm recalcd ", sm(e,t)
                    gws(e,t) = gws(e,t-1) + inf_gw/n(e) !... and subtracted from gw
                    write(42,*) "gws calcd", gws(e,t)
                    if(gws(e,t)>gok(e)) then
                        excess_gw_vol = (gws(e,t)-gok(e))*n(e) + sm(e,t)
                        gws(e,t) = gok(e)
                        sm(e,t) = 0
                        sws(e,t) = sws(e,t) + excess_gw_vol
                        write(42,*) "gws recalcd", gws(e,t)
                    end if
                    epv(e,t) = (gok(e) - gws(e,t))*n(e)
                    if(sm(e,t)>epv(e,t)) then
                        sws(e,t) = sws(e,t) + (sm(e,t)-epv(e,t))
                        sm(e,t) = epv(e,t)
                    end if
                    L = gok(e) - gws(e,t)
                    sm_eq = vanGI_simps(L, 750) !!!gw-sm balancing: consider adding a convergence criteria here
                    write(42,*) "new sm_eq", sm_eq
                    k_inf_gw = kGW(min(sm(e,t)/epv(e,t), 1.0)*n(e), k(e), vanG_pars)*dt - max(inf_gw, 0.00) !subtract k_inf_gw already utilized and allow freely capilary rise beyond k_inf_gw
                    write(42,*) "k_inf_gw remaining", k_inf_gw
                    inf_gw = min(sm(e,t)-sm_eq, max(k_inf_gw*dt,0.0))
                    if(gws(e,t) + inf_gw/n(e) < bot(e)) then
                        inf_gw = - min(abs((gws(e,t) - bot(e)))*n(e), abs(k_inf_gw*dt))
                        if(sm(e,t)<0) then
                            sm(e,t) = 0
                        end if
                    end if
                    write(42,*) "addnl inf_gw", inf_gw
                    sm(e,t) = sm(e,t) - inf_gw
                    gws(e,t) = gws(e,t) + inf_gw/n(e)
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
        !$OMP END PARALLEL DO
    end do
    Qdiff = Qin - Qout
end subroutine