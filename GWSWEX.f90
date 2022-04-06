module helpers
    implicit none
contains
    function kSM(s)
        real*8, intent(in) :: s
        real*8 :: kSM
        real*8 :: ks, n
        ks = 3e-2
        n = 2.5
        kSM = ks*s*(1-(1-(s)**1/(1-(1/n)))**(1-(1/n)))**2
    end function kSM

    function kGW(s)
        real*8, intent(in) :: s
        real*8 :: kGW
        real*8 :: ks, n
        ks = 3e-2
        n = 2.5
        kGW = ks*s*(1-(1-(s)**1/(1-(1/n)))**(1-(1/n)))**2
    end function kGW
end module helpers

module gwswex
    USE OMP_LIB
    implicit none
    integer  :: elems, nts, dt
    logical, allocatable  :: chd(:)
    real*8 :: n
    real*8, allocatable :: gok(:), k(:), p(:,:), et(:,:)
contains

    subroutine build(el, ts, ts_size, gok_l, n_l)
        integer, intent(in) :: el, ts, ts_size, gok_l(:)
        real*8, intent(in) :: n_l
        elems = el
        nts = ts
        dt = ts_size
        allocate(gok(elems), chd(elems), p(elems,nts), et(elems,nts))
        gok = gok_l
        n = n_l
        open(unit=42, file="fort.log", status="replace")
        write(42,*) "built"
    end subroutine

    subroutine init(chd_l, p_l, et_l)
        real*8, intent(in) :: p_l(:,:), et_l(:,:)
        logical, intent(in) :: chd_l(:)
        p = p_l
        et = et_l
        chd = chd_l
        open(unit=42, file="fort.log", status="append")
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
        integer :: i, j
        open(unit=42, file="fort.log", status="old")
        write(42,*) "run entered"

        do j = 2, nts
            write(42,*) "outer loop entered. ts ", j-1
            do i = 1, elems
                write(42,*) "inner loop entered. elem", i
                if(.NOT. chd(i)) then
                    L = gok(i) - gws(i,j-1) !prev. GW depth
                    if(L<0 .OR. L==0) then !NO UZ case
                        write(42,*) "noUZ entered"
                        !excess GW correction
                        excess_gw_vol = -L*n + sm(i,j-1)
                        gws(i,j) = gok(i)
                        sm(i,j) = 0
                        epv(i,j) = 0
                        sws(i,j) = sws(i,j-1) + excess_gw_vol + p(i,j)*dt
                        !ET extraction
                        if (sws(i,j)>et(i,j)*dt) then
                            sws(i,j) = sws(i,j) - et(i,j)*dt
                        else
                            sw_et_deficit = et(i,j)*dt - sws(i,j)
                            sws(i,j) = 0
                            gws(i,j) = gws(i,j) - (sw_et_deficit/n)
                            epv(i,j) = (gok(i) - gws(i,j))*n
                        end if
                        !calc storage discharges
                        gw_dis(i,j) = (gws(i,j) - gws(i,j-1))*n
                        sm_dis(i,j) = (sm(i,j)) - sm(i,j-1)
                        sw_dis(i,j) = sws(i,j) - sws(i,j-1)
                        Qin(i,j) = p(i,j)*dt - et(i,j)*dt
                        Qout(i,j) = gw_dis(i,j) + sw_dis(i,j) + sm_dis(i,j)
                        sw_et_deficit = 0
                    else
                        write(42,*) "UZ entered"
                        !P dist and SW push
                        write(42,*) "L is", L
                        write(42,*) "P is", p(i,j)*dt
                        write(42,*) "sm is", sm(i,j-1)
                        write(42,*) "epv is", epv(i,j-1)
                        write(42,*) "sm/epv", sm(i,j-1)/epv(i,j-1)
                        k_inf = kSM(min(sm(i,j-1)/epv(i,j-1), 1.0)*n) !calc K from wetness at the begining of this dt i.e. end of last dt
                        write(42,*) "got k", k_inf*dt
                        inf = min(k_inf*dt, p(i,j)*dt)
                        write(42,*) "inf is", inf
                        excess_p = p(i,j)*dt - inf
                        write(42,*) "excess_p", excess_p
                        inf_deficit = k_inf*dt - inf
                        write(42,*) "inf_deficit", inf_deficit
                        write(42,*) "sws is", sws(i,j-1)
                        sw_inf = min(inf_deficit, sws(i,j-1))
                        sws(i,j) = sws(i,j-1) - sw_inf + excess_p
                        write(42,*) "sw_inf", sw_inf
                        write(42,*) "sws calcd", sws(i,j)
                        sm(i,j) = sm(i,j-1) + inf + sw_inf
                        write(42,*) "sm calcd", sm(i,j)
                        sm_eq = vanGI(-L) !!!consider doing GW push and gw-sm bal after ET extraction
                        write(42,*) "gws is ", gws(i,j-1)
                        write(42,*) "vanGI called. sm_eq is ", sm_eq
                        k_inf_gw = kGW(min(sm(i,j-1)/epv(i,j-1), 1.0)*n) !calc K from current wetness (after P and SW inf)
                        gw_inf = min(sm(i,j)-sm_eq, k_inf_gw*dt) !if sm<sm_eq, gw_inf is -ve ...
                        write(42,*) "k_inf_gw is", k_inf_gw*dt
                        write(42,*) "gw_inf is", gw_inf
                        sm(i,j) = sm(i,j) - gw_inf !... deficit sm gets added to sm from gw
                        write(42,*) "sm recalcd ", sm(i,j)
                        gws(i,j) = gws(i,j-1) + gw_inf/n !... and subtracted from gw
                        write(42,*) "gws calcd", gws(i,j)
                        if(gws(i,j)>gok(i)) then
                            excess_gw_vol = (gws(i,j)-gok(i))*n + sm(i,j)
                            gws(i,j) = gok(i)
                            sm(i,j) = 0
                            sws(i,j) = sws(i,j) + excess_gw_vol
                        end if
                        write(42,*) "gws recalcd", gws(i,j)

                        !ET removal and SM-GW rebalance
                        write(42,*) "ET is", et(i,j)*dt
                        sw_et = min(sws(i,j), et(i,j)*dt)
                        sws(i,j) = sws(i,j) - sw_et
                        et_deficit = et(i,j)*dt - sw_et
                        write(42,*) "sw et removed", sw_et
                        sm(i,j) = sm(i,j) - et_deficit
                        write(42,*) "sm et removed", et_deficit
                        write(42,*) "sm is", sm(i,j)
                        sm_eq = vanGI(-(gok(i) - gws(i,j))) !!!gw-sm balancing: consider adding a convergence criteria here
                        write(42,*) "new sm_eq", sm_eq
                        k_inf_gw = kGW(min(sm(i,j-1)/epv(i,j-1), 1.0)*n)*dt - max(gw_inf, 0.00) !subtract k_inf_gw already utilized and allow freely capilary rise beyond k_inf_gw
                        write(42,*) "k_inf_gw remaining", k_inf_gw
                        gw_inf = min(sm(i,j)-sm_eq, k_inf_gw*dt)
                        write(42,*) "addnl gw_inf", gw_inf
                        sm(i,j) = sm(i,j) - gw_inf
                        gws(i,j) = gws(i,j) + gw_inf/n
                        write(42,*) "sm-gw balanced", sm(i,j), gws(i,j)

                        epv(i,j) = (gok(i) - gws(i,j))*n
                        gw_dis(i,j) = (gws(i,j) - gws(i,j-1))*n
                        sw_dis(i,j) = sws(i,j) - (sws(i,j-1))
                        sm_dis(i,j) = sm(i,j) - sm(i,j-1)
                        Qin(i,j) = p(i,j)*dt - et(i,j)*dt
                        Qout(i,j) = gw_dis(i,j) + sw_dis(i,j) + sm_dis(i,j)
                    end if
                else
                    excess_gw_vol = sm(i,j-1)
                    gws(i,j) = gws(i,j-1)
                    sm(i,j) = 0
                    epv(i,j) = 0
                    sws(i,j) = sws(i,j-1) + p(i,j)*dt - et(i,j)*dt + excess_gw_vol		
                    gw_dis(i,j) = 0
                    sw_dis(i,j) = sws(i,j) - sws(i,j-1)
                    sm_dis(i,j) = 0
                    Qin(i,j) = p(i,j)*dt - et(i,j)*dt
                    Qout(i,j) = gw_dis(i,j) + sw_dis(i,j) + sm_dis(i,j)
                end if
            end do
        end do
        Qdiff = Qin - Qout
    end subroutine
end module gwswex