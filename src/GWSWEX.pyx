import numpy as np
from scipy.integrate import quad

cdef float kSM(float s, float ks, vanG_pars):
    # add Horton_inf regime effects. only capilary/vG-M considered here
    cdef double theta_r, theta_s, n, m, sat
    theta_r = vanG_pars[0]
    theta_s = vanG_pars[1]
    n = vanG_pars[3]
    m = (1-(1/n))
    if(s<0.1):
        sat = ((0.1-theta_r)/(theta_s-theta_r))
    else:
        sat = ((s-theta_r)/(theta_s-theta_r))
    kSM = ks*sat*((1-(1-(sat)**(1/m))**m)**2)
    return kSM

cdef float kGW(s, ks, vanG_pars):
    # add preferential flow and implement kGW(s,d)
    cdef double kGW, theta_r, theta_s, n, m, sat
    theta_r = vanG_pars[0]
    theta_s = vanG_pars[1]
    n = vanG_pars[3]
    m = (1-(1/n))
    if(s<0.1):
        sat = ((0.1-theta_r)/(theta_s-theta_r))
    else:
        sat = ((s-theta_r)/(theta_s-theta_r))
    kGW = ks*sat*((1-(1-(sat)**(1/m))**m)**2)
    return kGW

cdef double vanGI(d, vanG_pars):
    d = d/100
    def theta(h_c):
        theta_r = vanG_pars[0]
        theta_s = vanG_pars[1]
        alpha = vanG_pars[2]
        n = vanG_pars[3]
        m = (1-(1/n))
        return np.float64(theta_r + ((theta_s - theta_r)/((1+(alpha*(abs(h_c)))**n))**m))
    return np.float64(quad(theta,d,0)[0])*100

cdef tuple crun(gws_ini, sws_ini, sm_ini, epv_ini, nts, elems, n, dt, k, bot, chd, p, et, gok, vanG_pars):
    cdef int e, t
    cdef double L, sw_et_deficit, excess_gw_vol, sm_eq, k_inf, inf, excess_p, inf_deficit, sw_inf, k_inf_gw, inf_gw, et_deficit, sw_et
    
    cdef double[:,:] gws = np.empty((elems,nts), dtype=np.double, order="c")
    cdef double[:,:] sws = np.empty((elems,nts), dtype=np.double, order="c")
    cdef double[:,:] sm = np.empty((elems,nts), dtype=np.double, order="c")
    cdef double[:,:] epv = np.empty((elems,nts), dtype=np.double, order="c")
    cdef double[:,:] gw_dis = np.empty((elems,nts), dtype=np.double, order="c")
    cdef double[:,:] sw_dis = np.empty((elems,nts), dtype=np.double, order="c")
    cdef double[:,:] sm_dis = np.empty((elems,nts), dtype=np.double, order="c")
    cdef double[:,:] Qin = np.empty((elems,nts), dtype=np.double, order="c")
    cdef double[:,:] Qout = np.empty((elems,nts), dtype=np.double, order="c")
    cdef double[:,:] Qdiff = np.empty((elems,nts), dtype=np.double, order="c")

    gws.base[:,0] = np.array(gws_ini, dtype=np.double, order="c")
    sws.base[:,0] = np.array(sws_ini, dtype=np.double, order="c")
    sm.base[:,0] = np.array(sm_ini, dtype=np.double, order="c")
    epv.base[:,0] = np.array(epv_ini, dtype=np.double, order="c")

    for t in np.arange(1, nts):
        for e in range(elems):
            if(not chd[e]):
                L = gok[e] - gws[e,t-1] #prev. GW depth
                if(L<0 or L==0): #NO UZ case
                    #excess GW correction
                    excess_gw_vol = -L*n[e] + sm[e,t-1]
                    gws[e][t] = gok[e]
                    sm[e][t] = 0
                    epv[e][t] = 0
                    sws[e][t] = sws[e,t-1] + excess_gw_vol + p[e][t]*dt

                    #ET extraction
                    if (sws[e][t]>et[e][t]*dt):
                        sws[e][t] = sws[e][t] - et[e][t]*dt
                    else:
                        sw_et_deficit = et[e][t]*dt - sws[e][t]
                        sws[e][t] = 0
                        gws[e][t] = gws[e][t] - (sw_et_deficit/n[e])
                        epv[e][t] = (gok[e] - gws[e][t])*n[e]

                    #calc storage discharges
                    gw_dis[e][t] = (gws[e][t] - gws[e,t-1])*n[e]
                    sm_dis[e][t] = (sm[e][t]) - sm[e,t-1]
                    sw_dis[e][t] = sws[e][t] - sws[e,t-1]
                    Qin[e][t] = p[e][t]*dt - et[e][t]*dt
                    Qout[e][t] = gw_dis[e][t] + sw_dis[e][t] + sm_dis[e][t]
                    sw_et_deficit = 0

                else:
                    #P dist and SW push
                    k_inf = kSM(min(sm[e,t-1]/epv[e,t-1], 1.0)*n[e], k[e], vanG_pars) #calc K from wetness at the begining of this dt i.e. end of last dt
                    inf = min(k_inf*dt, p[e][t]*dt)
                    excess_p = p[e][t]*dt - inf
                    sw_et = min(sws[e,t-1]+excess_p, et[e][t]*dt)
                    inf_deficit = k_inf*dt - inf
                    sw_inf = min(inf_deficit, sws[e,t-1]+excess_p-sw_et)
                    sws[e][t] = sws[e,t-1] - sw_inf + excess_p - sw_et
                    et_deficit = et[e][t]*dt - sw_et
                    if(gws[e,t-1] <= bot[e]):
                        et_deficit = 0
                    sm[e][t] = sm[e,t-1] + inf + sw_inf - et_deficit
                    sm_eq = vanGI(L,vanG_pars)
                    k_inf_gw = kGW(min(sm[e][t]/epv[e,t-1], 1.0)*n[e], k[e], vanG_pars) #calc K from current wetness (after P and SW inf)
                    inf_gw = min(sm[e][t]-sm_eq, k_inf_gw*dt) #if sm<sm_eq, inf_gw is -ve ...
                    if(gws[e,t-1] + inf_gw/n[e] < bot[e]):
                        inf_gw = - min(abs((gws[e,t-1] - bot[e]))*n[e], abs(k_inf_gw*dt))
                    sm[e][t] = sm[e][t] - inf_gw #... deficit sm gets added to sm from gw
                    gws[e][t] = gws[e,t-1] + inf_gw/n[e] #... and subtracted from gw
                    if(gws[e][t]>gok[e]):
                        excess_gw_vol = (gws[e][t]-gok[e])*n[e] + sm[e][t]
                        gws[e][t] = gok[e]
                        sm[e][t] = 0
                        sws[e][t] = sws[e][t] + excess_gw_vol
                    epv[e][t] = (gok[e] - gws[e][t])*n[e]
                    if(sm[e][t]>epv[e][t]):
                        sws[e][t] = sws[e][t] + (sm[e][t]-epv[e][t])
                        sm[e][t] = epv[e][t]
                    L = gok[e] - gws[e][t]
                    sm_eq = vanGI(L,vanG_pars) ###gw-sm balancing: consider adding a convergence criteria here
                    k_inf_gw = kGW(min(sm[e][t]/epv[e][t], 1.0)*n[e], k[e], vanG_pars)*dt - max(inf_gw, 0.00) #subtract k_inf_gw already utilized and allow freely capilary rise beyond k_inf_gw
                    inf_gw = min(sm[e][t]-sm_eq, max(k_inf_gw*dt,0.0))
                    if(gws[e][t] + inf_gw/n[e] < bot[e]):
                        inf_gw = - min(abs((gws[e][t] - bot[e]))*n[e], abs(k_inf_gw*dt))
                        if(sm[e][t]<0):
                            sm[e][t] = 0
                    sm[e][t] = sm[e][t] - inf_gw
                    gws[e][t] = gws[e][t] + inf_gw/n[e]
                    
                    epv[e][t] = (gok[e] - gws[e][t])*n[e]
                    gw_dis[e][t] = (gws[e][t] - gws[e,t-1])*n[e]
                    sw_dis[e][t] = sws[e][t] - (sws[e,t-1])
                    sm_dis[e][t] = sm[e][t] - sm[e,t-1]
                    Qin[e][t] = p[e][t]*dt - et[e][t]*dt
                    Qout[e][t] = gw_dis[e][t] + sw_dis[e][t] + sm_dis[e][t]
            else:
                excess_gw_vol = sm[e,t-1]
                gws[e][t] = gws[e,t-1]
                sm[e][t] = 0
                epv[e][t] = 0
                sws[e][t] = sws[e,t-1] + p[e][t]*dt - et[e][t]*dt + excess_gw_vol		
                gw_dis[e][t] = 0
                sw_dis[e][t] = sws[e][t] - sws[e,t-1]
                sm_dis[e][t] = 0
                Qin[e][t] = p[e][t]*dt - et[e][t]*dt
                Qout[e][t] = gw_dis[e][t] + sw_dis[e][t] + sm_dis[e][t]
    Qdiff = np.array(Qin) - np.array(Qout)
    return (gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff)

def run(gws_ini, sws_ini, sm_ini, epv_ini, nts, elems, n, dt, k, bot, chd, p, et, gok, vanG_pars):
    return crun(gws_ini, sws_ini, sm_ini, epv_ini, nts, elems, n, dt, k, bot, chd, p, et, gok, vanG_pars)