module gwswex
    USE OMP_LIB
    implicit none

    integer  :: elems, nts, dt
    logical, allocatable  :: chd(:)
    real(8) :: n, m, n_gw, beta, alpha, sw_th
    real(8), allocatable :: gok(:), k(:), p(:), et(:), gws_ini(:), sws_ini(:), epv_ini(:), sm_ini(:)
contains

    subroutine build(el, ts, ts_size, gok_l)
        integer, intent(in) :: el, ts, ts_size, gok_l(:)
        elems = el
        nts = ts
        dt = ts_size
        allocate(gok(elems), k(elems), chd(elems), p(nts), et(nts), gws_ini(elems), sws_ini(elems), epv_ini(elems), sm_ini(elems))
        gok = gok_l
    end subroutine

    subroutine params(n_l, m_l, n_gw_l, beta_l, alpha_l, sw_th_l, k_l)
        real(8), intent(in) :: n_l, m_l, n_gw_l, beta_l, alpha_l, sw_th_l, k_l(:)
        n = n_l
        m = m_l
        n_gw = n_gw_l
        beta = beta_l
        alpha = alpha_l
        sw_th = sw_th_l
        k = k_l
    end subroutine

    subroutine init(chd_l, p_l, et_l, gws_ini_l, sws_ini_l, epv_ini_l, sm_ini_l)
        real(8), intent(in) :: p_l(:), et_l(:), gws_ini_l(:), sws_ini_l(:), epv_ini_l(:), sm_ini_l(:)
        logical, intent(in) :: chd_l(:)
        chd = chd_l
        p = p_l
        et = et_l
        gws_ini = gws_ini_l
        sws_ini = sws_ini_l
        epv_ini = epv_ini_l
        sm_ini = sm_ini_l
    end subroutine

    subroutine run(gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff)
    !f2py threadsafe
        real(8), intent(inout) :: gws(:,:), sws(:,:), sm(:,:), epv(:,:), gw_dis(:,:), sw_dis(:,:), sm_dis(:,:), &
            Qin(:,:), Qout(:,:), Qdiff(:,:)
        real(8) :: max_p_gw, pPlus, p_gw, p_sw, p_sm, pwp, gw_to_sm, L, excess_gw_vol,&
        sw_et_deficit, sm_to_et, sm_to_gw, sw_to_gw, gw_to_et
        integer :: i, j
        gws(:,1) = gws_ini(:)
        sws(:,1) = sws_ini(:)
        epv(:,1) = epv_ini(:)
        sm(:,1) = sm_ini(:)
        gw_dis(:,1) = 0
        sw_dis(:,1) = 0
        sm_dis(:,1) = 0
        Qin(:,1) = 0
        Qout(:,1) = 0

        !$OMP PARALLEL DO
        do j = 2, nts
            do i = 1, elems
                if(.NOT. chd(i)) then
                    L = gok(i) - gws(i,j-1)
                    if(L<0 .OR. L==0) then
                        excess_gw_vol = -L*n_gw + sm(i,j-1)
                        gws(i,j) = gok(i)
                        sm(i,j) = 0
                        epv(i,j) = 0
                        
                        sws(i,j) = sws(i,j-1) + excess_gw_vol + p(j)*dt
                        if (sws(i,j)>et(j)*dt) then
                            sws(i,j) = sws(i,j) - et(j)*dt
                        else
                            sw_et_deficit = et(j)*dt - sws(i,j)
                            sws(i,j) = 0
                            gws(i,j) = gws(i,j) - (sw_et_deficit/n)
                            epv(i,j) = (gok(i) - gws(i,j))*n
                        end if
                        
                        gw_dis(i,j) = (gws(i,j) - gws(i,j-1))*n
                        sm_dis(i,j) = (sm(i,j)) - sm(i,j-1)
                        sw_dis(i,j) = sws(i,j) - sws(i,j-1)
                        Qin(i,j) = p(j)*dt - et(j)*dt
                        Qout(i,j) = gw_dis(i,j) + sw_dis(i,j) + sm_dis(i,j)
                        sw_et_deficit = 0
                    else
                        epv(i,j) = L*n
                        pwp = epv(i,j)*m
                        
                        !P dist
                        pPlus = p(j)*dt
                        p_sm = pPlus * (1 - (sm(i,j-1)/epv(i,j))**beta)
                        if(p_sm<0) then
                            p_sm = 0
                        end if
                        p_gw = min((min(k(i)*dt,pPlus)*(sm(i,j-1)/epv(i,j))**alpha), pPlus-p_sm, k(i)*dt)
                        p_sw = pPlus - p_sm - p_gw
                        
                        !SW Push
                        max_p_gw = k(i)*dt*((sm(i,j-1)/epv(i,j))**alpha)
                        if(p_gw<max_p_gw .AND. sws(i,j-1)>sw_th) then
                            if(sws(i,j-1)>max_p_gw-p_gw) then
                                sw_to_gw = max_p_gw - p_gw - sw_th
                            else
                                sw_to_gw = sws(i,j-1) - sw_th
                            end if
                            sws(i,j-1) = sws(i,j-1) - sw_to_gw
                            p_gw = p_gw + sw_to_gw
                        end if
                        
                        !SM calc
                        sm(i,j) = sm(i,j-1) + p_sm
                        !SW calc and ET removal
                        sws(i,j) = sws(i,j-1) + p_sw
                        if(sws(i,j)>et(j)*dt+sw_th) then
                            sws(i,j) = sws(i,j) - et(j)*dt
                        else if(sm(i,j)>et(j)*dt) then
                            sm(i,j) = sm(i,j) - et(j)*dt
                        else
                            gw_to_et = et(j)*dt
                            gws(i,j-1) = gws(i,j-1) - (gw_to_et/n)
                        end if
                        
                        !GW calc
                        gws(i,j) = gws(i,j-1) + (p_gw/n)
                        epv(i,j) = (gok(i) - gws(i,j))*n
                        pwp = epv(i,j)*m
                        
                        !SM adjustments
                        if(sm(i,j)<pwp) then
                            gw_to_sm = pwp - sm(i,j)
                            sm(i,j) = pwp
                        else
                            gw_to_sm = 0
                        end if
                        if(sm(i,j)>epv(i,j)) then
                            sm_to_gw = sm(i,j) - epv(i,j)
                            sm(i,j) = epv(i,j)
                        else
                            sm_to_gw = 0
                        end if
                        
                        !GW re-calc
                        gws(i,j) = gws(i,j) - (gw_to_sm/n) + (sm_to_gw/n)
                        epv(i,j) = (gok(i) - gws(i,j))*n
                        if (gws(i,j)>gok(i)) then
                            sws(i,j) = sws(i,j) + (gws(i,j)-gok(i))*n_gw + sm(i,j)
                            gws(i,j) = gok(i)
                            sm(i,j) = 0
                            epv(i,j) = 0
                        end if
                        
                        gw_dis(i,j) = (gws(i,j) - (gws(i,j-1)+gw_to_et/n))*n
                        sw_dis(i,j) = sws(i,j) - (sws(i,j-1)+sw_to_gw)
                        sm_dis(i,j) = sm(i,j) - sm(i,j-1)
                        Qin(i,j) = p(j)*dt - et(j)*dt
                        Qout(i,j) = gw_dis(i,j) + sw_dis(i,j) + sm_dis(i,j)
                        gw_to_et = 0
                        sw_to_gw = 0
                    end if
                else
                    excess_gw_vol = sm(i,j-1)
                    gws(i,j) = gws(i,j-1)
                    sm(i,j) = 0
                    epv(i,j) = 0
                    sws(i,j) = sws(i,j-1) + p(j)*dt - et(j)*dt + excess_gw_vol		
                    gw_dis(i,j) = 0
                    sw_dis(i,j) = sws(i,j) - sws(i,j-1)
                    sm_dis(i,j) = 0
                    Qin(i,j) = p(j)*dt - et(j)*dt
                    Qout(i,j) = gw_dis(i,j) + sw_dis(i,j) + sm_dis(i,j)
                end if
            end do
        end do
        !$OMP END PARALLEL DO
        Qdiff = Qin - Qout
    end subroutine
end module gwswex