function kGW(s, ks)
	! add preferential flow and implement kGW(s,d)
		real*8, intent(in) :: s, ks
		real*8 :: kGW, theta_r, theta_s, n, m, sat
		theta_r = vanG_pars(1)
		theta_s = vanG_pars(2)
		n = vanG_pars(4)
		m = (1-(1/n))
		if(s<0.1) then
			sat = ((0.1-theta_r)/(theta_s-theta_r))
		else
			sat = ((s-theta_r)/(theta_s-theta_r))
		end if
		kGW = ks*sat*((1-(1-(sat)**(1/m))**m)**2)
end function kGW