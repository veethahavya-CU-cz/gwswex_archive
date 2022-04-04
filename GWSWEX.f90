module gwswex
    USE OMP_LIB
    implicit none
    integer  :: elems, nts, dt
    logical, allocatable  :: chd(:)
    real*8 :: n, sy, theta_s, theta_r, alpha, beta, sw_th
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
        write(*,*) "built"
    end subroutine

    subroutine init(chd_l, p_l, et_l)
        real*8, intent(in) :: p_l(:,:), et_l(:,:)
        logical, intent(in) :: chd_l(:)
        p = p_l
        et = et_l
        chd = chd_l
        write(*,*) "initialised"
    end subroutine

    subroutine run(vanGI, kSM, gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff)
    !f2py threadsafe
        external :: vanGI
        real*8 :: vanGI
        !f2py real*8, intent(in):: d
        !f2py real*8, intent(out) :: eq
        !f2py eq = vanGI(d)
        external :: kSM
        real*8 :: kSM
        !f2py real*8, intent(in) :: sm
        !f2py real*8, intent(out):: k
        !f2py k = kSM(sm)
        real*8, intent(inout) :: gws(:,:), sws(:,:), sm(:,:), epv(:,:), gw_dis(:,:), sw_dis(:,:), sm_dis(:,:), &
            Qin(:,:), Qout(:,:), Qdiff(:,:)
        real*8 :: L, sw_et_deficit, excess_gw_vol, sm_eq, k_inf, inf, excess_p, inf_deficit, sw_inf, &
            k_inf_gw, gw_inf, et_deficit, sw_et
        integer :: i, j
        write(*,*) "run entered"

        do j = 2, nts
            write(*,*) "outer loop entered. ts ", j-1
            !$OMP PARALLEL PRIVATE(i), SHARED(j, gws,sws,sm,epv, gw_dis,sw_dis,sm_dis)
            !OMP SHARED(L,sw_et_deficit,excess_gw_vol,sm_eq,k_inf,inf,excess_p,inf_deficit,sw_inf,k_inf_gw,gw_inf,et_deficit,sw_et)
            !$OMP DO 
            do i = 1, elems
                write(*,*) "inner loop entered. elem", i
                if(.NOT. chd(i)) then
                    L = gok(i) - gws(i,j-1) !prev. GW depth
                    if(L<0 .OR. L==0) then !NO UZ case
                        write(*,*) "noUZ entered"
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
                        write(*,*) "UZ entered"
                        !P dist and SW push
                        write(*,*) "L is", L
                        write(*,*) "P is", p(i,j)*dt
                        write(*,*) "sm is", sm(i,j-1)
                        k_inf = kSM(sm(i,j-1)) !calc K from wetness at the begining of this dt i.e. end of last dt
                        write(*,*) "got k", k_inf*dt
                        inf = min(k_inf*dt, p(i,j)*dt)
                        write(*,*) "inf is", inf
                        excess_p = p(i,j)*dt - inf
                        write(*,*) "excess_p", excess_p
                        inf_deficit = k_inf*dt - inf
                        write(*,*) "inf_deficit", inf_deficit
                        write(*,*) "sws is", sws(i,j-1)
                        sw_inf = min(inf_deficit, sws(i,j-1))
                        sws(i,j) = sws(i,j-1) - sw_inf + excess_p
                        write(*,*) "sw_inf", sw_inf
                        write(*,*) "sws calcd", sws(i,j)
                        sm(i,j) = sm(i,j-1) + inf + sw_inf
                        write(*,*) "sm calcd", sm(i,j)
                        sm_eq = vanGI(-L) !!!consider doing GW push and gw-sm bal after ET extraction
                        write(*,*) "gws is ", gws(i,j-1)
                        write(*,*) "vanGI called. sm_eq is ", sm_eq
                        k_inf_gw = kSM(sm(i,j)) !calc K from current wetness (after P and SW inf)
                        gw_inf = min(sm(i,j)-sm_eq, k_inf_gw*dt) !if sm<sm_eq, gw_inf is -ve ...
                        write(*,*) "k_inf_gw is", k_inf_gw*dt
                        write(*,*) "gw_inf is", gw_inf
                        sm(i,j) = sm(i,j) - gw_inf !... deficit sm gets added to sm from gw
                        write(*,*) "sm recalcd ", sm(i,j)
                        gws(i,j) = gws(i,j-1) + gw_inf/n !... and subtracted from gw
                        write(*,*) "gws calcd", gws(i,j)
                        if(gws(i,j)>gok(i)) then
                            excess_gw_vol = (gws(i,j)-gok(i))*n + sm(i,j)
                            gws(i,j) = gok(i)
                            sm(i,j) = 0
                            sws(i,j) = sws(i,j) + excess_gw_vol
                        end if
                        write(*,*) "gws recalcd", gws(i,j)

                        !ET removal and SM-GW rebalance
                        write(*,*) "ET is", et(i,j)*dt
                        sw_et = min(sws(i,j), et(i,j)*dt)
                        sws(i,j) = sws(i,j) - sw_et
                        et_deficit = et(i,j)*dt - sw_et
                        write(*,*) "sw et removed", sw_et
                        sm(i,j) = sm(i,j) - et_deficit
                        write(*,*) "sm et removed", et_deficit
                        write(*,*) "sm is", sm(i,j)
                        sm_eq = vanGI(-(gok(i) - gws(i,j))) !!!gw-sm balancing: consider adding a convergence criteria here
                        write(*,*) "new sm_eq", sm_eq
                        k_inf_gw = kSM(sm(i,j))*dt - max(gw_inf, 0.00) !subtract k_inf_gw already utilized and allow freely capilary rise beyond k_inf_gw
                        write(*,*) "k_inf_gw remaining", k_inf_gw
                        gw_inf = min(sm(i,j)-sm_eq, k_inf_gw*dt)
                        write(*,*) "addnl gw_inf", gw_inf
                        sm(i,j) = sm(i,j) - gw_inf
                        gws(i,j) = gws(i,j) + gw_inf/n
                        write(*,*) "sm-gw balanced", sm(i,j), gws(i,j)

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
            !$OMP END DO
            !$OMP END PARALLEL
        end do
        Qdiff = Qin - Qout
    end subroutine
end module gwswex