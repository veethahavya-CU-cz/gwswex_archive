subroutine run()
	USE helpers, only: kSM, kGW, vanGI_fgsl

	implicit none

	real*8 :: L, sw_et_deficit, excess_gw_vol, sm_eq, k_inf, inf, excess_p, inf_deficit, sw_inf, &
		k_inf_gw, inf_gw, et_deficit, sw_et, start, finish
	integer :: e, t

	character(255) :: gws_file, gws_ini_file, sws_file, sws_ini_file, sm_file, sm_ini_file, epv_file, epv_ini_file, &
		gw_dis_file, sw_dis_file, sm_dis_file, Qin_file, Qout_file, Qdiff_file

	gws_ini_file = trim(input_path)//"gws_ini.ip"
	sws_ini_file = trim(input_path)//"sws_ini.ip"
	sm_ini_file = trim(input_path)//"sm_ini.ip"
	epv_ini_file = trim(input_path)//"epv_ini.ip"

	write(lu,*) "run entered"

	open(tu, file=gws_ini_file, form='unformatted', action='read')
	read(tu) gws(:,1)
	close(tu, status='keep')
	open(tu, file=sws_ini_file, form='unformatted', action='read')
	read(tu) sws(:,1)
	close(tu, status='keep')
	open(tu, file=sm_ini_file, form='unformatted', action='read')
	read(tu) sm(:,1)
	close(tu, status='keep')
	open(tu, file=epv_ini_file, form='unformatted', action='read')
	read(tu) epv(:,1)
	close(tu, status='keep')

	do t = 2, nts
		write(lu,*) "outer loop entered. ts ", t-1
		!$OMP PARALLEL DO SHARED(gws, sws, sm, epv) &
		!$OMP PRIVATE(L, sw_et_deficit, excess_gw_vol, sm_eq, k_inf, inf, excess_p, inf_deficit, sw_inf, k_inf_gw, inf_gw) &
		!$OMP PRIVATE(et_deficit, sw_et)
		!!! PUT ALL CODE THAT NEEDS TO BE SEQUENTIAL WITHIN FUNCTIONS AND CALL THEM ACCORDINGLY
		do e = 1, elems
			write(lu,*) "inner loop entered. elem", e
			write(lu,*) "gok", gok(e)
			write(lu,*) "bot", bot(e)
			if(.NOT. chd(e)) then
				L = gok(e) - gws(e,t-1) !prev. GW depth
				if(L<0 .OR. L==0) then !NO UZ case
					write(lu,*) "noUZ entered"
					!excess GW correction
					write(lu,*) "gws is ", gws(e,t-1)
					write(lu,*) "sws is ", sws(e,t-1)
					write(lu,*) "sm is ", sm(e,t-1)
					excess_gw_vol = -L*n(e) + sm(e,t-1)
					gws(e,t) = gok(e)
					sm(e,t) = 0
					epv(e,t) = 0
					sws(e,t) = sws(e,t-1) + excess_gw_vol + p(e,t)*dt
					write(lu,*) "excess_gw_vol ", excess_gw_vol
					write(lu,*) "gws after +p ", gws(e,t)
					write(lu,*) "sws after +p", sws(e,t)
					write(lu,*) "sm after +p", sm(e,t)
					!ET extraction
					if (sws(e,t)>et(e,t)*dt) then
						sws(e,t) = sws(e,t) - et(e,t)*dt
					else
						sw_et_deficit = et(e,t)*dt - sws(e,t)
						sws(e,t) = 0
						gws(e,t) = gws(e,t) - (sw_et_deficit/n(e))
						epv(e,t) = (gok(e) - gws(e,t))*n(e)
					end if
					write(lu,*) "gws after -et ", gws(e,t)
					write(lu,*) "sws after -et", sws(e,t)
					write(lu,*) "sm after -et", sm(e,t)
					!calc storage discharges
					gw_dis(e,t) = (gws(e,t) - gws(e,t-1))*n(e)
					sm_dis(e,t) = (sm(e,t)) - sm(e,t-1)
					sw_dis(e,t) = sws(e,t) - sws(e,t-1)
					Qin(e,t) = p(e,t)*dt - et(e,t)*dt
					Qout(e,t) = gw_dis(e,t) + sw_dis(e,t) + sm_dis(e,t)
					sw_et_deficit = 0
				else
					write(lu,*) "UZ entered"
					!P dist and SW push
					write(lu,*) "L is", L
					write(lu,*) "P is", p(e,t)*dt
					write(lu,*) "sm is", sm(e,t-1)
					write(lu,*) "epv is", epv(e,t-1)
					write(lu,*) "sm/epv", sm(e,t-1)/epv(e,t-1)
					k_inf = kSM(min(sm(e,t-1)/epv(e,t-1), 1.0)*n(e), k(e)) !calc K from wetness at the begining of this dt i.e. end of last dt
					write(lu,*) "got k", k_inf
					inf = min(k_inf*dt, p(e,t)*dt)
					write(lu,*) "inf aka p_sm is ", inf
					excess_p = p(e,t)*dt - inf
					write(lu,*) "excess p aka p_sw is ", excess_p
					write(lu,*) "sws is", sws(e,t-1)
					write(lu,*) "ET is", et(e,t)*dt
					sw_et = min(sws(e,t-1)+excess_p, et(e,t)*dt)
					inf_deficit = k_inf*dt - inf
					write(lu,*) "inf_deficit", inf_deficit
					sw_inf = min(inf_deficit, sws(e,t-1)+excess_p-sw_et)
					write(lu,*) "sw_inf", sw_inf
					sws(e,t) = sws(e,t-1) - sw_inf + excess_p - sw_et
					et_deficit = et(e,t)*dt - sw_et
					if(gws(e,t-1) <= bot(e)) then
						et_deficit = 0
					end if
					write(lu,*) "sw et removed", sw_et
					write(lu,*) "sws calcd", sws(e,t)
					sm(e,t) = sm(e,t-1) + inf + sw_inf - et_deficit
					write(lu,*) "sm et removed", et_deficit
					write(lu,*) "sm calcd", sm(e,t)
					call cpu_time(start)
					sm_eq = vanGI_fgsl(L)
					call cpu_time(finish)
					write(lu,*) "vanGI_fgsl time ", finish-start
					write(lu,*) "gws is ", gws(e,t-1)
					write(lu,*) "vanGI_fgsl called. sm_eq is ", sm_eq
					k_inf_gw = kGW(min(sm(e,t)/epv(e,t-1), 1.0)*n(e), k(e)) !calc K from current wetness (after P and SW inf)
					inf_gw = min(sm(e,t)-sm_eq, k_inf_gw*dt) !if sm<sm_eq, inf_gw is -ve ...
					if(gws(e,t-1) + inf_gw/n(e) < bot(e)) then
						inf_gw = - min(abs((gws(e,t-1) - bot(e)))*n(e), abs(k_inf_gw*dt))
					end if
					write(lu,*) "k_inf_gw is", k_inf_gw
					write(lu,*) "inf_gw is", inf_gw
					sm(e,t) = sm(e,t) - inf_gw !... deficit sm gets added to sm from gw
					write(lu,*) "sm recalcd ", sm(e,t)
					gws(e,t) = gws(e,t-1) + inf_gw/n(e) !... and subtracted from gw
					write(lu,*) "gws calcd", gws(e,t)
					if(gws(e,t)>gok(e)) then
						excess_gw_vol = (gws(e,t)-gok(e))*n(e) + sm(e,t)
						gws(e,t) = gok(e)
						sm(e,t) = 0
						sws(e,t) = sws(e,t) + excess_gw_vol
						write(lu,*) "gws recalcd", gws(e,t)
					end if
					epv(e,t) = (gok(e) - gws(e,t))*n(e)
					if(sm(e,t)>epv(e,t)) then
						sws(e,t) = sws(e,t) + (sm(e,t)-epv(e,t))
						sm(e,t) = epv(e,t)
					end if
					L = gok(e) - gws(e,t)
					sm_eq = vanGI_fgsl(L) !!!gw-sm balancing: consider adding a convergence criteria here
					write(lu,*) "new sm_eq", sm_eq
					k_inf_gw = kGW(min(sm(e,t)/epv(e,t), 1.0)*n(e), k(e))*dt - max(inf_gw, 0.00) !subtract k_inf_gw already utilized and allow freely capilary rise beyond k_inf_gw
					write(lu,*) "k_inf_gw remaining", k_inf_gw
					inf_gw = min(sm(e,t)-sm_eq, max(k_inf_gw*dt,0.0))
					if(gws(e,t) + inf_gw/n(e) < bot(e)) then
						inf_gw = - min(abs((gws(e,t) - bot(e)))*n(e), k_inf_gw*dt)
						if(sm(e,t)<0) then
							sm(e,t) = 0
						end if
					end if
					write(lu,*) "addnl inf_gw", inf_gw
					sm(e,t) = sm(e,t) - inf_gw
					gws(e,t) = gws(e,t) + inf_gw/n(e)
					write(lu,*) "sm-gw balanced", sm(e,t), gws(e,t)

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

	gws_file = trim(output_path)//"gws.op"
	gw_dis_file = trim(output_path)//"gw_dis.op"
	sws_file = trim(output_path)//"sws.op"
	sw_dis_file = trim(output_path)//"sw_dis.op"
	sm_file = trim(output_path)//"sm.op"
	sm_dis_file = trim(output_path)//"sm_dis.op"
	epv_file = trim(output_path)//"epv.op"
	Qin_file = trim(output_path)//"Qin.op"
	Qout_file = trim(output_path)//"Qout.op"
	Qdiff_file = trim(output_path)//"Qdiff.op"

	open(tu, file=gws_file, form='unformatted', action='write')
	write(tu) gws
	close(tu)
	open(tu, file=gw_dis_file, form='unformatted', action='write')
	write(tu) gw_dis
	close(tu)
	open(tu, file=sws_file, form='unformatted', action='write')
	write(tu) sws
	close(tu)
	open(tu, file=sw_dis_file, form='unformatted', action='write')
	write(tu) sw_dis
	close(tu)
	open(tu, file=sm_file, form='unformatted', action='write')
	write(tu) sm
	close(tu)
	open(tu, file=sm_dis_file, form='unformatted', action='write')
	write(tu) sm_dis
	close(tu)
	open(tu, file=epv_file, form='unformatted', action='write')
	write(tu) epv
	close(tu)
	open(tu, file=Qin_file, form='unformatted', action='write')
	write(tu) Qin
	close(tu)
	open(tu, file=Qout_file, form='unformatted', action='write')
	write(tu) Qout
	close(tu)
	open(tu, file=Qdiff_file, form='unformatted', action='write')
	write(tu) Qdiff
	close(tu)

end subroutine