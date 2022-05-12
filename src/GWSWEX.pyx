import numpy as np
import os
from scipy.integrate import quad

with open("cyGWSWEX.log", "w"):
    pass

log = open("cyGWSWEX.log", "a")
cdef logit(str msg):
    log.write(msg + "\n")

cdef float kSM(float s, float ks, vanG_pars):
    # add Horton_inf regime effects. only capilary/vG-M considered here
    cdef double theta_r, theta_s, n, m, sat
    theta_r = vanG_pars[0]
    theta_s = vanG_pars[1]
    n = vanG_pars[3]
    m = (1-(1/n))
    if(s<0.1):
        sat = ((0.1-theta_r)/(theta_s-theta_r))
    elif(s<theta_s):
        sat = ((s-theta_r)/(theta_s-theta_r))
    else:
        return ks
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
    elif(s<theta_s):
        sat = ((s-theta_r)/(theta_s-theta_r))
    else:
        return ks
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

cpdef tuple run(gws_ini, sws_ini, sm_ini, epv_ini, nts, elems, n, dt, k, bot, chd, p, et, gok, vanG_pars):
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
        logit("\n\n*** ts - {} ***".format(t+1))
        for e in range(elems):
            logit("\n* elem - {} ***".format(e+1))
            if(not chd[e]):
                L = gok[e] - gws[e,t-1] #prev. GW depth
                logit("L is "+str(L))
                if(L<0 or L==0): #NO UZ case
                    #excess GW correction
                    logit("noUZ entered")
                    logit("sws was "+str(sws[e][t-1]))
                    logit("gws was "+str(gws[e][t-1]))
                    logit("sm was "+str(sm[e][t-1]))
                    logit("p is "+str(p[e][t]))
                    logit("et is "+str(et[e][t]))
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
                    logit("sws is "+str(sws[e][t]))
                    logit("gws is "+str(gws[e][t]))
                    logit("sm is "+str(sm[e][t]))

                    #calc storage discharges
                    gw_dis[e][t] = (gws[e][t] - gws[e,t-1])*n[e]
                    sm_dis[e][t] = (sm[e][t]) - sm[e,t-1]
                    sw_dis[e][t] = sws[e][t] - sws[e,t-1]
                    Qin[e][t] = p[e][t]*dt - et[e][t]*dt
                    Qout[e][t] = gw_dis[e][t] + sw_dis[e][t] + sm_dis[e][t]
                    sw_et_deficit = 0

                else:
                    #P dist and SW push
                    logit("UZ entered")
                    logit("sws was "+str(sws[e][t-1]))
                    logit("gws was "+str(gws[e][t-1]))
                    logit("sm was "+str(sm[e][t-1]))
                    logit("sat is "+str(min(sm[e,t-1]/epv[e,t-1], 1.0)*n[e]))
                    k_inf = kSM(min(sm[e,t-1]/epv[e,t-1], 1.0)*n[e], k[e], vanG_pars) #calc K from wetness at the begining of this dt i.e. end of last dt
                    logit("Ksm is "+str(k_inf))
                    inf = min(k_inf*dt, p[e][t]*dt)
                    logit("inf is "+str(inf))
                    excess_p = p[e][t]*dt - inf
                    logit("excess_p is "+str(excess_p))
                    sw_et = min(sws[e,t-1]+excess_p, et[e][t]*dt)
                    logit("sw_et is "+str(sw_et))
                    inf_deficit = k_inf*dt - inf
                    logit("inf_deficit is "+str(inf_deficit))
                    sw_inf = min(inf_deficit, sws[e,t-1]+excess_p-sw_et)
                    logit("sw_inf is "+str(sw_inf))
                    sws[e][t] = sws[e,t-1] - sw_inf + excess_p - sw_et
                    logit("sws is "+str(sws[e][t]))
                    et_deficit = et[e][t]*dt - sw_et
                    logit("et_deficit is "+str(et_deficit))
                    if(gws[e,t-1] <= bot[e]):
                        et_deficit = 0
                    sm[e][t] = sm[e,t-1] + inf + sw_inf - et_deficit
                    logit("sm is "+str(sm[e][t]))
                    sm_eq = vanGI(L,vanG_pars)
                    logit("sm_eq is "+str(sm_eq))
                    logit("sat is "+str(min(sm[e,t]/epv[e,t-1], 1.0)*n[e]))
                    k_inf_gw = kGW(min(sm[e][t]/epv[e,t-1], 1.0)*n[e], k[e], vanG_pars) #calc K from current wetness (after P and SW inf)
                    logit("k_inf_gw is "+str(k_inf_gw))
                    inf_gw = min(sm[e][t]-sm_eq, k_inf_gw*dt) #if sm<sm_eq, inf_gw is -ve ...
                    logit("inf_gw calcd "+str(inf_gw))
                    if(gws[e,t-1] + inf_gw/n[e] < bot[e]):
                        inf_gw = - min(abs((gws[e,t-1] - bot[e]))*n[e], abs(k_inf_gw*dt))
                    logit("inf_gw is "+str(inf_gw))
                    sm[e][t] = sm[e][t] - inf_gw #... deficit sm gets added to sm from gw
                    logit("sm is "+str(sm[e][t]))
                    gws[e][t] = gws[e,t-1] + inf_gw/n[e] #... and subtracted from gw
                    logit("gws is "+str(gws[e][t]))
                    if(gws[e][t]>gok[e]):
                        excess_gw_vol = (gws[e][t]-gok[e])*n[e] + sm[e][t]
                        gws[e][t] = gok[e]
                        sm[e][t] = 0
                        sws[e][t] = sws[e][t] + excess_gw_vol
                    logit("gws is "+str(gws[e][t]))
                    epv[e][t] = (gok[e] - gws[e][t])*n[e]
                    logit("epv is "+str(epv[e][t]))
                    if(sm[e][t]>epv[e][t]):
                        sws[e][t] = sws[e][t] + (sm[e][t]-epv[e][t])
                        sm[e][t] = epv[e][t]
                    logit("sm is "+str(sm[e][t]))
                    L = gok[e] - gws[e][t]
                    sm_eq = vanGI(L,vanG_pars) ###gw-sm balancing: consider adding a convergence criteria here
                    logit("sm_eq is "+str(sm_eq))
                    if not epv[e][t] == 0:
                        logit("sat is "+str(min(sm[e][t]/epv[e][t], 1.0)*n[e])) #ZERO DIVISION HERE WHEN SM=ePV=0!!!!!!!!!!!!!!!!!!!!!!!!!
                        k_inf_gw = kGW(min(sm[e][t]/epv[e][t], 1.0)*n[e], k[e], vanG_pars)*dt - max(inf_gw, 0.00) #subtract k_inf_gw already utilized and allow freely capilary rise beyond k_inf_gw
                    else:
                        k_inf_gw = k[e]
                    logit("k_inf_gw is "+str(k_inf_gw))
                    inf_gw = min(sm[e][t]-sm_eq, max(k_inf_gw*dt,0.0))
                    logit("inf_gw is "+str(inf_gw))
                    if(gws[e][t] + inf_gw/n[e] < bot[e]):
                        inf_gw = - min(abs((gws[e][t] - bot[e]))*n[e], abs(k_inf_gw*dt))
                        if(sm[e][t]<0):
                            sm[e][t] = 0
                    logit("inf_gw is "+str(inf_gw))
                    sm[e][t] = sm[e][t] - inf_gw
                    gws[e][t] = gws[e][t] + inf_gw/n[e]
                    logit("gws is "+str(gws[e][t]))
                    logit("sm is "+str(sm[e][t]))

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
    log.close()
    return (gws, sws, sm, epv, gw_dis, sw_dis, sm_dis, Qin, Qout, Qdiff)