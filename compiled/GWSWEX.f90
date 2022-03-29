program GWSWEX
	implicit none
	real(8), dimension(9) :: args
	integer :: elems, ts, ts_size
	real(8), allocatable :: gok(:), k(:), p(:), et(:), gws(:,:), gws_ini(:), sws(:,:), sws_ini(:),&
	sm(:,:), sm_ini(:), epv(:,:), epv_ini(:), gw_dis(:,:), sw_dis(:,:), sm_dis(:,:), Qin(:,:), Qout(:,:), Qdiff(:,:)
	real(8) :: n, n_gw, m, alpha, beta, max_p_gw, pPlus, p_gw, p_sw, p_sm, pwp, gw_to_sm, L, excess_gw_vol,&
	sw_et_deficit, sm_deficit, sm_to_gw, sw_to_gw, gw_to_et, sw_th
	integer :: i, j
	integer(1), allocatable :: chd(:)
	
	open(9, file="args.ip")
	read(9, *) args
	elems = int(args(1))
	ts = int(args(2))
	ts_size = int(args(3))
	n = args(4)
	n_gw = args(5)
	m = args(6)
	beta = args(7)
	alpha = args(8)
	sw_th = args(9)
	
	allocate(gok(elems), k(elems), gws_ini(elems), sws_ini(elems), sm_ini(elems), epv_ini(elems))
	allocate(gws(elems,ts), sws(elems,ts), sm(elems,ts), epv(elems,ts), gw_dis(elems,ts), sw_dis(elems,ts), sm_dis(elems,ts),&
	Qin(elems,ts), Qout(elems,ts), Qdiff(elems,ts), chd(elems))
	allocate(p(ts), et(ts))

	open(10, file="gok.ip", form='unformatted')
	read(10) gok
	open(11, file="k.ip", form='unformatted')
	read(11) k
	open(12, file="p.ip", form='unformatted')
	read(12) p
	open(13, file="et.ip", form='unformatted')
	read(13) et
	open(14, file="gws_ini.ip", form='unformatted')
	read(14) gws_ini
	open(15, file="sws_ini.ip", form='unformatted')
	read(15) sws_ini
	open(16, file="epv_ini.ip", form='unformatted')
	read(16) epv_ini
	open(17, file="sm_ini.ip", form='unformatted')
	read(17) sm_ini
	open(18, file="chd.ip", form='unformatted')
	read(18) chd
	
	gws(:,1) = gws_ini(:)
	sws(:,1) = sws_ini(:)
	epv(:,1) = epv_ini(:)
	sm(:,1) = sm_ini(:)
	gw_dis(:,1) = 0
	sw_dis(:,1) = 0
	sm_dis(:,1) = 0
	Qin(:,1) = 0
	Qout(:,1) = 0
	
	do i = 1, elems
		if(chd(i) == 0) then
			do j = 2, ts
				L = gok(i) - gws(i,j-1)
				if(L<0 .OR. L==0) then
					excess_gw_vol = -L*n_gw + sm(i,j-1)
					gws(i,j) = gok(i)
					sm(i,j) = 0
					epv(i,j) = 0
					
					sws(i,j) = sws(i,j-1) + excess_gw_vol + p(j)*ts_size
					if (sws(i,j)>et(j)*ts_size) then
						sws(i,j) = sws(i,j) - et(j)*ts_size
					else
						sw_et_deficit = et(j)*ts_size - sws(i,j)
						sws(i,j) = 0
						gws(i,j) = gws(i,j) - (sw_et_deficit/n)
						epv(i,j) = (gok(i) - gws(i,j))*n
					end if
					
					gw_dis(i,j) = (gws(i,j) - gws(i,j-1))*n
					sm_dis(i,j) = (sm(i,j)) - sm(i,j-1)
					sw_dis(i,j) = sws(i,j) - sws(i,j-1)
					Qin(i,j) = p(j)*ts_size - et(j)*ts_size
					Qout(i,j) = gw_dis(i,j) + sw_dis(i,j) + sm_dis(i,j)
					sw_et_deficit = 0
				else
					epv(i,j) = L*n
					pwp = epv(i,j)*m
					
					!P dist
					pPlus = p(j)*ts_size
					p_sm = pPlus * (1 - (sm(i,j-1)/epv(i,j))**beta)
					if(p_sm<0) then
						p_sm = 0
					end if
					p_gw = min((min(k(i)*ts_size,pPlus)*(sm(i,j-1)/epv(i,j))**alpha), pPlus-p_sm, k(i)*ts_size)
					p_sw = pPlus - p_sm - p_gw
					
					!SW Push
					max_p_gw = k(i)*ts_size*((sm(i,j-1)/epv(i,j))**alpha)
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
					if(sws(i,j)>et(j)*ts_size+sw_th) then
						sws(i,j) = sws(i,j) - et(j)*ts_size
					else if(sm(i,j)>et(j)*ts_size) then
						sm(i,j) = sm(i,j) - et(j)*ts_size
					else
						gw_to_et = et(j)*ts_size
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
					Qin(i,j) = p(j)*ts_size - et(j)*ts_size
					Qout(i,j) = gw_dis(i,j) + sw_dis(i,j) + sm_dis(i,j)
					gw_to_et = 0
					sw_to_gw = 0
				end if
			end do
		else
			do j = 2, ts
				excess_gw_vol = sm(i,j-1)
				gws(i,j) = gws(i,j-1)
				sm(i,j) = 0
				epv(i,j) = 0
				sws(i,j) = sws(i,j-1) + p(j)*ts_size - et(j)*ts_size + excess_gw_vol		
				gw_dis(i,j) = 0
				sw_dis(i,j) = sws(i,j) - sws(i,j-1)
				sm_dis(i,j) = 0
				Qin(i,j) = p(j)*ts_size - et(j)*ts_size
				Qout(i,j) = gw_dis(i,j) + sw_dis(i,j) + sm_dis(i,j)
			end do
		end if	
	end do
	Qdiff = Qin - Qout
	open(20, file="gws.op", form="unformatted")
	write(20) gws
	open(21, file="sws.op", form="unformatted")
	write(21) sws
	open(22, file="sm.op", form="unformatted")
	write(22) sm
	open(23, file="epv.op", form="unformatted")
	write(23) epv
	open(24, file="gw_dis.op", form="unformatted")
	write(24) gw_dis
	open(25, file="sw_dis.op", form="unformatted")
	write(25) sw_dis
	open(26, file="sm_dis.op", form="unformatted")
	write(26) sm_dis
	open(27, file="Qin.op", form="unformatted")
	write(27) Qin
	open(28, file="Qout.op", form="unformatted")
	write(28) Qout
	open(29, file="Qdiff.op", form="unformatted")
	write(29) Qdiff
end program