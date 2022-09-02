""" Generate spdata file used by GTC and Orbit.

---------------
    Author: Xishuo Wei. (xishuow@uci.edu, weixishuo@gmail.com)
    
"""

import numpy as np
from scipy import interpolate
from geqdsk import read
from fortran import spline_function_f
import argparse
import logging
import pickle
import os


# import sys
# path = __file__
# this_dir = os.path.dirname(path)
# sys.path.append(os.path.join(this_dir,"RBF/"))
# import pyximport
# pyximport.install()
# from rbf.interpolate import RBFInterpolant


def mapping(filename,lsp,lst, psimax_ratio, figs,nR=200):
    import matplotlib
    if figs=="save":
        matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    npsi = lsp - 1
    ntheta = lst - 1 # max(400,lst) - 1
    nR = np.max((nR,npsi,ntheta))    # grids used for inverse interpolation

    logging.info(f"Reading {filename}")
    data = read(filename)
    R2d = data["r"]
    Z2d = data["z"]

    psi2d = data["psi"]-data["simagx"]
    R2d_full = data["r"]
    Z2d_full = data["z"]
    psi_full = data["psi"]
    
    logging.info("construting psi(R,Z) RBF interpolants...")
    psi_RZ_f = interpolate.Rbf(R2d,Z2d,psi2d,kind="cubic")# interpolate.interp2d(R2d,Z2d,psi2d.T,kind="cubic")
    R_flat = R2d.flatten()
    Z_flat = Z2d.flatten()
    RZ_input = np.concatenate((R_flat[:,None],Z_flat[:,None]),axis=1)
    psi_input = psi2d.flatten()
    #psi_RBF_rep = RBFInterpolant(RZ_input, psi_input)
    logging.info("done!")

    # levels = [data["sibdry"]-data["simagx"]]
    # fig, ax = plt.subplots(figsize=(7,7),dpi=100)
    # CS = ax.contour(R2d,Z2d,psi2ds,levels,colors="cyan",linewidths=2)
    # CS = ax.contour(R_sd2d,Z_sd2d,psi_superdense,levels,colors="black",linestyles="dashed",linewidths=1)
    # ax.plot(data["rbdry"], data["zbdry"],'r--')
    # ax.set_aspect('equal', 'box')
    # ax.plot(data["xlim"], data["ylim"],'k-',linewidth=0.5)

    R_ctr = data['rmagx']
    Z_ctr = data['zmagx']
    theta_reps = []
    psimin = 0
    psimax = (data["sibdry"]-data["simagx"])*psimax_ratio#psi_RZ_f(max(data["rbdry"]),Z_ctr).squeeze()

    delp = (psimax - psimin) / npsi
    delt = np.pi*2 / ntheta
    bdry_r = np.sqrt((data["rbdry"]-R_ctr)**2+(data["zbdry"]-Z_ctr)**2)
    boundary_theta = np.arccos((data["rbdry"]-R_ctr)/bdry_r)
    lower_points = data["zbdry"] < Z_ctr                                 
    boundary_theta[lower_points] = np.pi*2.0-boundary_theta[lower_points]
    for jt in range(len(data["zbdry"])//4,0,-1):                         
        if boundary_theta[jt]<boundary_theta[jt-1]:                      
            boundary_theta[jt-1] -= np.pi*2                              
    for jt in range(len(data["zbdry"])//4*3,len(boundary_theta)-1):      
        if boundary_theta[jt]>boundary_theta[jt+1]:                      
            boundary_theta[jt+1] += np.pi*2                              

    trace_ind = np.arange(nR)
    psimesh = np.linspace(psimin, psimax, npsi+1)
    Ronpt = np.zeros((npsi+1,ntheta+1)) + R_ctr
    Zonpt = np.zeros((npsi+1,ntheta+1)) + Z_ctr
    jt = 0
    tmp_psi1d = np.zeros(nR)

    logging.info("Inverse interpolation to get R(psi,theta0) and Z(psi,theta0)")
    for theta in np.linspace(0, np.pi*2, ntheta+1):

        nearest_bdry_point = np.argmin(np.abs(boundary_theta-theta))
        longr = np.linspace(0,bdry_r[nearest_bdry_point]*2,nR*2)
        longR = longr * np.cos(theta) + R_ctr
        longZ = longr * np.sin(theta) + Z_ctr

        psi_prev=-1000.0
        for iilongR in range(nR*2):
            R_loc = longR[iilongR]
            Z_loc = longZ[iilongR]
            psi_loc = psi_RZ_f(R_loc,Z_loc)
            if psi_loc >= data["sibdry"]-data["simagx"] or psi_loc <= psi_prev :
                break
            psi_prev = psi_loc
            
        tmpr = np.linspace(0,longr[iilongR],nR)
        tmpR = tmpr * np.cos(theta) + R_ctr
        tmpZ = tmpr * np.sin(theta) + Z_ctr
    #     tmp_psi = psi_RZ_f(tmpR, tmpZ)
        for ind in range(nR):
            tmp_psi1d[ind] = psi_RZ_f(tmpR[ind],tmpZ[ind])
        for ind in range(nR-2,-1,-1):
            if tmp_psi1d[ind] > tmp_psi1d[ind+1]:
                tmp_psi1d[ind] = tmp_psi1d[ind+1]-2e-6
        inv_interpR = interpolate.CubicSpline(tmp_psi1d, tmpR)
        inv_interpZ = interpolate.CubicSpline(tmp_psi1d, tmpZ)

        tmpR_onpt = inv_interpR(psimesh[1:])
        tmpZ_onpt = inv_interpZ(psimesh[1:])

        Ronpt[1:, jt] = tmpR_onpt
        Zonpt[1:, jt] = tmpZ_onpt
        jt +=1


    tmesh = np.linspace(0, np.pi*2, ntheta+1)
    for ip in range(1,npsi+1):
        Ront = Ronpt[ip,:]
        Zont = Zonpt[ip,:]
        Rrep = interpolate.splrep(tmesh, Ront, per=1, s=0)
        Zrep = interpolate.splrep(tmesh, Zont, per=1, s=0)
        theta_reps.append((Rrep,Zrep))
    logging.info("done!")

    BR = np.zeros((npsi+1,ntheta+1))
    BZ = np.zeros((npsi+1,ntheta+1))
    bp_dir = np.zeros((ntheta+1,2))
    logging.info("Taking local psip to get Bp...")
    # altenative method, use local psi to get Bp
    for ip in range(1,npsi+1):
        bp_dir[1:ntheta,:] = np.array((Ronpt[ip, 2:] - Ronpt[ip, :ntheta-1], Zonpt[ip,2:] - Zonpt[ip, :ntheta-1])).T
        bp_dir[0] = np.array([Ronpt[ip, 1] - Ronpt[ip, -2], Zonpt[ip, 1] - Zonpt[ip, -2]])
        bp_dir[-1] = bp_dir[0]
        nbp = bp_dir / np.sqrt(bp_dir[:,0]**2+bp_dir[:,1]**2)[:,None]
        ndl = nbp[:,[1,0]]
        ndl[:,0] *= -1
        length=0.01
        dl = length * ndl
        R_inner = Ronpt[ip,:] - dl[:,0] / 2
        R_outer = Ronpt[ip,:] + dl[:,0] / 2
        Z_inner = Zonpt[ip,:] - dl[:,1] / 2
        Z_outer = Zonpt[ip,:] + dl[:,1] / 2

        psi_inner = psi_RZ_f(R_inner,Z_inner)
        psi_outer = psi_RZ_f(R_outer,Z_outer)
        Bp = np.abs(psi_outer-psi_inner) / (Ronpt[ip,:]*length)
        BR[ip,:] = Bp*nbp[:,0]
        BZ[ip,:] = Bp*nbp[:,1]
        BR[ip,ntheta] = BR[ip,0]
        BZ[ip,ntheta] = BZ[ip,0]
    # for ip in range(1,npsi+1):
    #     for jt in range(0, ntheta):
    #         if jt==0:
    #             bp_dir = np.array([Ronpt[ip, jt+1] - Ronpt[ip, -2], Zonpt[ip, jt+1] - Zonpt[ip, -2]])
    #         else:
    #             bp_dir = np.array([Ronpt[ip, jt+1] - Ronpt[ip, jt-1], Zonpt[ip, jt+1] - Zonpt[ip, jt-1]])
    #         nbp = bp_dir / np.sqrt(bp_dir[0]**2+bp_dir[1]**2)
    #         ndl = np.array((-nbp[1], nbp[0]))
    #         length = 0.01
    #         dl = length*ndl
    #         R_inner = Ronpt[ip, jt]-dl[0]/2
    #         R_outer = Ronpt[ip, jt]+dl[0]/2
    #         Z_inner = Zonpt[ip, jt]-dl[1]/2
    #         Z_outer = Zonpt[ip, jt]+dl[1]/2

    #         psi_inner = psi_RZ_f(R_inner,Z_inner)
    #         psi_outer = psi_RZ_f(R_outer,Z_outer)
    #         Bp = np.abs(psi_outer-psi_inner)/(Ronpt[ip,jt]*length)
    #         BR[ip,jt] = Bp*nbp[0]
    #         BZ[ip,jt] = Bp*nbp[1]
    #     BR[ip,ntheta] = BR[ip,0]
    #     BZ[ip,ntheta] = BZ[ip,0]

    logging.info("Taking derivative of psi to get Bp...")
    #BR_from_psi = np.zeros((npsi+1,ntheta+1))
    #BZ_from_psi = np.zeros((npsi+1,ntheta+1))
    #
    #for ip in range(1,npsi+1):
    #    for jt in range(0, ntheta):
    #        R_loc = Ronpt[ip,jt]
    #        Z_loc = Zonpt[ip,jt]
    #        dpsi_dR = psi_RBF_rep([[R_loc, Z_loc]], diff=(1,0))
    #        dpsi_dZ = psi_RBF_rep([[R_loc, Z_loc]], diff=(0,1))
    #        # dpsi_dR = interpolate.bisplev(R_loc, Z_loc, psi_RZ_rep, dx=1)
    #        # dpsi_dZ = interpolate.bisplev(R_loc, Z_loc, psi_RZ_rep, dy=1)
    #        BZ_from_psi[ip,jt] = dpsi_dR/R_loc
    #        BR_from_psi[ip,jt] = -dpsi_dZ/R_loc
    #    BR_from_psi[ip,ntheta] = BR_from_psi[ip,0]
    #    BZ_from_psi[ip,ntheta] = BZ_from_psi[ip,0]

    logging.info("done!")

    # Construct Bt from fpol
    BT = np.zeros((npsi+1,ntheta+1))
    gcurrent = data["fpol"]
    dpgrids = np.linspace(0, data["sibdry"]-data["simagx"],len(gcurrent))

    if np.sign(gcurrent[0]) < 0:
        logging.warning("Negative g detected. Enforce it to be positive...")
    fit_g = interpolate.CubicSpline(dpgrids, gcurrent)
    smooth_g = abs(fit_g(psimesh))

    fit_q = interpolate.CubicSpline(dpgrids, data["qpsi"])
    smooth_q = fit_q(psimesh)

    fit_p = interpolate.CubicSpline(dpgrids, data["pressure"])
    smooth_p = fit_p(psimesh)

    for ip in range(0,npsi+1):
        for jt in range(ntheta+1):
            tmpR = Ronpt[ip,jt]
            BT[ip,jt] = smooth_g[ip] / tmpR

    # Trace field lines
    theta_trace = []
    phi_trace = []
    q_from_trace = [smooth_q[0]]

    logging.info("Tracing field lines...")
    for ip in range(1, npsi+1):
        R_fit = interpolate.CubicSpline(tmesh, Ronpt[ip,:])
        Z_fit = interpolate.CubicSpline(tmesh, Zonpt[ip,:])
    #     Bp_ont = np.sqrt(BR_from_psi[ip,:]**2 + BZ_from_psi[ip,:]**2)
    #    BR_fit = interpolate.CubicSpline(tmesh, BR_from_psi[ip,:])
    #    BZ_fit = interpolate.CubicSpline(tmesh, BZ_from_psi[ip,:])
        BR_fit = interpolate.CubicSpline(tmesh, BR[ip,:])
        BZ_fit = interpolate.CubicSpline(tmesh, BZ[ip,:])
        Bt_fit = interpolate.CubicSpline(tmesh, BT[ip,:])
        # start from theta=0
        theta = 0
        phi = 0
        dphi = 0.01 * smooth_q[ip]
        theta_this = [theta]
        phi_this = [phi]
        max_iter = 2500
        iiter = 0
        while abs(theta) < np.pi*2 and iiter < max_iter:
            iiter += 1
            tt = np.mod(theta,np.pi*2)
    #         Bp_loc = Bp_fit(tt)
            BR_loc = BR_fit(tt)
            BZ_loc = BZ_fit(tt)
            Bp_loc = np.sqrt(BR_loc**2+BZ_loc**2)
            Bt_loc = Bt_fit(tt)
            B_loc = np.sqrt(Bp_loc**2+Bt_loc**2)
            r_loc = np.sqrt((R_fit(tt)-R_ctr)**2+(Z_fit(tt)-Z_ctr)**2)
            R_loc = R_fit(tt)
            phi_mid = phi+dphi / 2
            theta_mid = theta+abs(dphi * (R_loc*Bp_loc)/(r_loc*Bt_loc) * (-BR_loc*np.sin(tt)+BZ_loc*np.cos(tt))/Bp_loc) / 2
            tt = np.mod(theta_mid,np.pi*2)
            Bt_loc = Bt_fit(tt)
            B_loc = np.sqrt(Bp_loc**2+Bt_loc**2)
            r_loc = np.sqrt((R_fit(tt)-R_ctr)**2+(Z_fit(tt)-Z_ctr)**2)
            R_loc = R_fit(tt)
            phi += dphi
            # note that Bp does not parallel to grad_\theta
            theta += abs(dphi * (R_loc*Bp_loc)/(r_loc*Bt_loc) * (-BR_loc*np.sin(tt)+BZ_loc*np.cos(tt))/Bp_loc)


            theta_this.append(theta)
            phi_this.append(phi)
        if iiter >= max_iter - 1:
            logging.error("unfinished iter at", ip)
        phi_this[-1] = (phi_this[-1] - phi_this[-2])/(theta_this[-1] - theta_this[-2]) * (np.pi*2-theta_this[-2]) + phi_this[-2]
        theta_this[-1] = np.pi*2


        theta_trace.append(theta_this)
        phi_trace.append(phi_this)
        q_from_trace.append(phi_this[-1]/theta_this[-1])

    logging.info("done!")

    phi_end = [smooth_q[0]]
    for ip in range(npsi):
        phi_end.append(phi_trace[ip][-1]/np.pi/2)

    #plt.plot(psimesh,phi_end)
    #plt.plot(psimesh,smooth_q,'--')

    # generate functions on theta_prime
    tpmesh = np.linspace(0, np.pi*2, ntheta+1)
    Ronptp = np.zeros((npsi+1,ntheta+1)) + R_ctr
    Zonptp = np.zeros((npsi+1,ntheta+1)) + Z_ctr
    t0_interp_end = []
    BR_ontp = np.zeros((npsi+1,ntheta+1))
    BZ_ontp = np.zeros((npsi+1,ntheta+1))
    t0onptp = np.zeros((npsi+1,ntheta+1))
    t0onptp[0,:] = tpmesh
    for ip in range(1,npsi+1):
        # construct t on tprime using trace data
        tp_this = np.array(phi_trace[ip-1]) / ((phi_trace[ip-1][-1])/theta_trace[ip-1][-1])
        t0_this = np.array(theta_trace[ip-1])
        t0_fit = interpolate.CubicSpline(tp_this, t0_this)
        t0_interp = t0_fit(tpmesh)
        t0onptp[ip] = t0_interp
        t0_interp_end.append(t0_interp[-1])
        Ronptp[ip, :] = interpolate.splev(t0_interp, theta_reps[ip-1][0])
        Zonptp[ip, :] = interpolate.splev(t0_interp, theta_reps[ip-1][1])
        Ronptp[ip, -1] = Ronptp[ip, 0]
        Zonptp[ip, -1] = Zonptp[ip, 0]
    #    BR_fit = interpolate.CubicSpline(tmesh, BR_from_psi[ip,:],bc_type='periodic')
    #    BZ_fit = interpolate.CubicSpline(tmesh, BZ_from_psi[ip,:],bc_type='periodic')
        BR_fit = interpolate.CubicSpline(tmesh, BR[ip,:],bc_type='periodic')
        BZ_fit = interpolate.CubicSpline(tmesh, BZ[ip,:],bc_type='periodic')
        BR_ontp[ip,:] = BR_fit(t0_interp)
        BZ_ontp[ip,:] = BZ_fit(t0_interp)

        BR_ontp[ip,-1] = BR_ontp[ip,0]
        BZ_ontp[ip,-1] = BZ_ontp[ip,0]
        Ronptp[ip,-1] = Ronptp[ip,0]
        Zonptp[ip,-1] = Zonptp[ip,0]

    # calculate nu function
    logging.info("calculating cpsi and nu ...")
    cpsi = np.zeros(npsi+1)
    delt = np.pi*2/ntheta
    orig_ri = np.zeros((npsi+1,ntheta+1))
    for ip in range(1,npsi+1):
    #     Bp_ont = np.sqrt(BR_from_psi[ip,:]**2 + BZ_from_psi[ip,:]**2)
        BR_fit = interpolate.CubicSpline(tpmesh, BR_ontp[ip,:],bc_type="periodic")
        BZ_fit = interpolate.CubicSpline(tpmesh, BZ_ontp[ip,:],bc_type="periodic")
        R_tp_rep = interpolate.splrep(tpmesh, Ronptp[ip,:], per=1)
        Z_tp_rep = interpolate.splrep(tpmesh, Zonptp[ip,:], per=1)
        ri = 0
        for jt in range(0,ntheta):
            tt = tpmesh[jt]+ delt / 2
            R_loc = interpolate.splev(tt, R_tp_rep)
            Z_loc = interpolate.splev(tt, Z_tp_rep)
            dR = interpolate.splev(tt, R_tp_rep, der=1)
            dZ = interpolate.splev(tt, Z_tp_rep, der=1)
            ri_loc = (BR_fit(tt)*dR+BZ_fit(tt)*dZ)
            ri += ri_loc * delt / (np.pi*2)
            orig_ri[ip,jt]=ri_loc
        cpsi[ip] = ri
        orig_ri[ip,-1] = orig_ri[ip,0]

    # Construct nu func for theta_b and zeta_b
    nu_onptp = np.zeros((npsi+1,ntheta+1))
    for ip in range(1,npsi+1):
    #     Bp_ont = np.sqrt(BR_from_psi[ip,:]**2 + BZ_from_psi[ip,:]**2)
        BR_fit = interpolate.CubicSpline(tpmesh, BR_ontp[ip,:],bc_type="periodic")
        BZ_fit = interpolate.CubicSpline(tpmesh, BZ_ontp[ip,:],bc_type="periodic")
        R_tp_rep = interpolate.splrep(tpmesh, Ronptp[ip,:], per=1)
        Z_tp_rep = interpolate.splrep(tpmesh, Zonptp[ip,:], per=1)
        nu_onptp[ip,0] = 0
        for jt in range(0,ntheta):
            tt = tpmesh[jt]+ delt / 2
            R_loc = interpolate.splev(tt, R_tp_rep)
            Z_loc = interpolate.splev(tt, Z_tp_rep)
            dR = interpolate.splev(tt, R_tp_rep, der=1)
            dZ = interpolate.splev(tt, Z_tp_rep, der=1)
            ri_loc = (BR_fit(tt)*dR+BZ_fit(tt)*dZ)
            nu_onptp[ip,jt+1] = nu_onptp[ip,jt] + (ri_loc - cpsi[ip]) * delt
        nu_onptp[ip,-1] = nu_onptp[ip,0]
        gqi = smooth_g[ip]+cpsi[ip]/smooth_q[ip]
        nu_onptp[ip,:] /= gqi
    
    dnudtp = np.zeros((npsi+1,ntheta+1))
    for ip in range(1,npsi+1):
        nu_tp_rep = interpolate.splrep(tpmesh, nu_onptp[ip,:], per=1)
        dnudtp[ip,:] = interpolate.splev(tpmesh,nu_tp_rep,der=1)

    logging.info("done!")

    logging.info("Interpolate values on theta_b...")
    BR_onptb = np.zeros((npsi+1,ntheta+1))
    BZ_onptb = np.zeros((npsi+1,ntheta+1))
    BT_onptb = np.zeros((npsi+1,ntheta+1)) + BT[0,0]
    tb_mesh = np.linspace(0,np.pi*2,ntheta+1)
    Ronptb = np.zeros((npsi+1,ntheta+1)) + R_ctr
    Zonptb = np.zeros((npsi+1,ntheta+1)) + Z_ctr
    nu_onptb = np.zeros((npsi+1,ntheta+1))
    delt = np.pi*2/(ntheta)
    for ip in range(1,npsi+1):
        #splines, a lot of
        tb_on_tp = tpmesh + nu_onptp[ip, :] / smooth_q[ip]
        R_fit = interpolate.CubicSpline(tb_on_tp, Ronptp[ip,:], bc_type="periodic")
        Z_fit = interpolate.CubicSpline(tb_on_tp, Zonptp[ip,:], bc_type="periodic")
        BR_fit = interpolate.CubicSpline(tb_on_tp, BR_ontp[ip,:], bc_type="periodic")
        BZ_fit = interpolate.CubicSpline(tb_on_tp, BZ_ontp[ip,:], bc_type="periodic")
        nu_fit = interpolate.CubicSpline(tb_on_tp, nu_onptp[ip,:], bc_type="periodic")

        Ronptb[ip,:] = R_fit(tb_mesh)
        Zonptb[ip,:] = Z_fit(tb_mesh)
        BR_onptb[ip,:] = BR_fit(tb_mesh)
        BZ_onptb[ip,:] = BZ_fit(tb_mesh)
        nu_onptb[ip,:] = nu_fit(tb_mesh)

        BT_onptb[ip,:] = smooth_g[ip] / Ronptb[ip,:]
    logging.info("done!")

    # logging.info("mapping from dense theta grid to lst theta grid")
    # BR_onptb_lst = np.zeros((lsp, lst))
    # BT_onptb_lst = np.zeros((lsp, lst)) + smooth_g[1] / R_ctr
    # BZ_onptb_lst = np.zeros((lsp, lst))
    # tb_mesh_lst = np.linspace(0,np.pi*2,lst)
    # Ronptb_lst = np.zeros((lsp, lst)) + R_ctr
    # Zonptb_lst = np.zeros((lsp, lst)) + Z_ctr

    # for ip in range(1, lsp):
    #     BR_fit = interpolate.CubicSpline(tb_mesh, BR_onptb[ip,:], bc_type='periodic')
    #     BZ_fit = interpolate.CubicSpline(tb_mesh, BZ_onptb[ip,:], bc_type='periodic')
    #     R_fit = interpolate.CubicSpline(tb_mesh, Ronptb[ip,:], bc_type='periodic')
    #     Z_fit = interpolate.CubicSpline(tb_mesh, Zonptb[ip,:], bc_type='periodic')

    #     BR_onptb_lst[ip,:] = BR_fit(tb_mesh_lst)
    #     BZ_onptb_lst[ip,:] = BZ_fit(tb_mesh_lst)
    #     Ronptb_lst[ip,:] = R_fit(tb_mesh_lst)
    #     Zonptb_lst[ip,:] = Z_fit(tb_mesh_lst)
    #     BT_onptb_lst[ip,:] = smooth_g / Ronptb_lst[ip,:]

    map_data = {
        "BR_onptb": BR_onptb,
        "BT_onptb": BT_onptb,
        "BZ_onptb": BZ_onptb,
        "psimesh": psimesh,
        "tb_mesh": tb_mesh,
        "nu_neg": -nu_onptb,
        "Ronptb": Ronptb,
        "Zonptb": Zonptb,
        "qpsi": smooth_q,
        "gpsi": smooth_g,
        "cpsi": cpsi,
        "ppsi": smooth_p,
        "R_ctr": R_ctr,
        "Z_ctr": Z_ctr,
        "lsp": lsp,
        "lst": lst
    }

    # check consistency
    # 1. q = B*nabla phi / (B* nabla tp) = BT*r/(R*dtpdt*(BZ*cos(tt)-BR*sin(tt)))
    q_2d_onptp = np.zeros((npsi+1,ntheta+1)) + smooth_q[0]
    for ip in range(1,npsi+1):
        tp_this = np.array(phi_trace[ip-1]) / q_from_trace[ip]
        t0_this = np.array(theta_trace[ip-1])
        tp_fit_t0 = interpolate.splrep(t0_this, tp_this, s=0)
        dtpdt0 = interpolate.splev(t0onptp[ip,:], tp_fit_t0, der=1)
        rr = np.sqrt((Ronptp[ip,:] - R_ctr)**2+(Zonptp[ip,:]-Z_ctr)**2)
        q_2d_onptp[ip, :] = smooth_g[ip]/Ronptp[ip,:]*rr/(Ronptp[ip,:]*dtpdt0*(BZ_ontp[ip,:]*np.cos(t0onptp[ip,:])-BR_ontp[ip,:]*np.sin(t0onptp[ip,:])))

    mean_q2d = np.mean(q_2d_onptp,axis=1)
    var_q2d = np.max(q_2d_onptp,axis=1) - np.min(q_2d_onptp,axis=1)
    fig = plt.figure(figsize=(8,6))
    plt.plot(psimesh,smooth_q,linewidth=2,color="tab:blue")
    plt.errorbar(psimesh,mean_q2d,var_q2d,color="tab:orange",ecolor="r",capsize=3, errorevery=5)
    plt.legend([r"$q$", r"$B\cdot{\nabla\zeta}/B\cdot{\nabla\theta}$"], fontsize=20)
    plt.xlabel(r"$\psi_p$")
    if figs=="save":
        plt.savefig("q_check.png")

    # 2. check whter \bar{I} is the covariant coefficient
    ri_2d_onptb = np.zeros((npsi+1,ntheta+1))
    dnudtb = np.zeros((npsi+1,ntheta+1))
    for ip in range(1,npsi+1):
    #     Bp_ont = np.sqrt(BR_from_psi[ip,:]**2 + BZ_from_psi[ip,:]**2)
        BR_fit = interpolate.CubicSpline(tb_mesh, BR_onptb[ip,:],bc_type="periodic")
        BZ_fit = interpolate.CubicSpline(tb_mesh, BZ_onptb[ip,:],bc_type="periodic")
        R_tb_rep = interpolate.splrep(tb_mesh, Ronptb[ip,:], per=1, s=0)
        Z_tb_rep = interpolate.splrep(tb_mesh, Zonptb[ip,:], per=1, s=0)
        nu_tb_rep = interpolate.splrep(tb_mesh, nu_onptb[ip,:], per=1, s=0)
        ri = 0
        for jt in range(0,ntheta):
            tt = tb_mesh[jt]
            R_loc = interpolate.splev(tt, R_tb_rep)
            Z_loc = interpolate.splev(tt, Z_tb_rep)
            dR = interpolate.splev(tt, R_tb_rep, der=1)
            dZ = interpolate.splev(tt, Z_tb_rep, der=1)
            dnu = interpolate.splev(tt, nu_tb_rep, der=1)
            dnudtb[ip,jt] = dnu
            ri_2d_onptb[ip,jt] = (BR_onptb[ip,jt]*dR+BZ_onptb[ip,jt]*dZ-BT_onptb[ip,jt]*dnu*Ronptb[ip,jt])
        ri_2d_onptb[ip,-1] = ri_2d_onptb[ip,0]

    ri_2d_onptb[ip,-1] = ri_2d_onptb[ip,0]
    mean_ri2d = np.mean(ri_2d_onptb,axis=1)
    var_ri2d = [-np.min(ri_2d_onptb,axis=1)+mean_ri2d, np.max(ri_2d_onptb,axis=1)-mean_ri2d]
    fig = plt.figure(figsize=(8,6))
    plt.plot(psimesh,cpsi,linewidth=2,color="tab:blue")
    plt.errorbar(psimesh,mean_ri2d,var_ri2d,color="tab:orange",ecolor="r",capsize=3, errorevery=5)
    plt.legend([r"$\bar{I}$", r"$B\cdot e_\theta$"], fontsize=20)
    plt.xlabel(r"$\psi_p$")
    if figs=="save":
        plt.savefig("I_check.png")

    # 3. Jacobian * B^2 .vs. gq+I
    dRdp = np.zeros((npsi+1,ntheta+1))
    dRdtb = np.zeros((npsi+1,ntheta+1))
    dZdp = np.zeros((npsi+1,ntheta+1))
    dZdtb = np.zeros((npsi+1,ntheta+1))
    dnudp = np.zeros((npsi+1,ntheta+1))
    dnudtb = np.zeros((npsi+1,ntheta+1))
    upper = np.zeros(ntheta+1)
    lower = np.zeros(ntheta+1)
    outer = np.zeros(npsi+1)
    delta = np.zeros((npsi+1,ntheta+1))
    inner = np.zeros(npsi+1)

    for ip in range(1,npsi+1):
        R_fit = interpolate.CubicSpline(tb_mesh, Ronptb[ip,:], bc_type="periodic")
        Z_fit = interpolate.CubicSpline(tb_mesh, Zonptb[ip,:], bc_type="periodic")
        nu_fit = interpolate.CubicSpline(tb_mesh, nu_onptb[ip,:], bc_type="periodic")
        dRdtb[ip,:] = R_fit(tb_mesh, nu=1)
        dZdtb[ip,:] = Z_fit(tb_mesh, nu=1)
        dnudtb[ip,:] = nu_fit(tb_mesh, nu=1)

    for jt in range(1,ntheta+1):
        # R_fit = interpolate.Akima1DInterpolator(psimesh, Ronptb[:,jt])
        # Z_fit = interpolate.Akima1DInterpolator(psimesh, Zonptb[:,jt])
        # nu_fit = interpolate.Akima1DInterpolator(psimesh, nu_onptb[:,jt])

        R_fit = interpolate.PchipInterpolator(psimesh[:-1], Ronptb[:-1,jt],extrapolate=True)
        Z_fit = interpolate.PchipInterpolator(psimesh[:-1], Zonptb[:-1,jt],extrapolate=True)
        nu_fit = interpolate.PchipInterpolator(psimesh[:-1], nu_onptb[:-1,jt],extrapolate=True)

        dp = (psimesh[1]-psimesh[0])/200
        outer = R_fit(psimesh[1:]+dp)
        inner = R_fit(psimesh[1:]-dp)
        dRdp[1:,jt] = (outer - inner) / (dp*2)

        outer = Z_fit(psimesh[1:]+dp)
        inner = Z_fit(psimesh[1:]-dp)
        dZdp[1:,jt] = (outer - inner) / (dp*2)

        outer = nu_fit(psimesh[1:]+dp)
        inner = nu_fit(psimesh[1:]-dp)
        dnudp[1:,jt] = (outer - inner) / (dp*2)

        delta[1:,jt] = BR_onptb[1:,jt]*dRdp[1:,jt] + BZ_onptb[1:,jt]*dZdp[1:,jt]-BT_onptb[1:,jt]*dnudp[1:,jt]*Ronptb[ip,jt]

    map_data["delta"] = delta

    dRdp[:,0] = dRdp[:,-1]
    dZdp[:,0] = dZdp[:,-1]
    dnudp[:,0] = dnudp[:,-1]
    Jac_fromcoord = np.zeros((npsi+1,ntheta+1))
    Jac_fromcoord = Ronptb*(dRdp*dZdtb-dRdtb*dZdp)
    Bmag2 = BZ_onptb**2+BR_onptb**2+BT_onptb**2
    Jac_B2 = Jac_fromcoord * Bmag2
    mean_Jacb2 = np.mean(Jac_B2, axis=1)
    var_Jacb2 = [(-np.min(Jac_B2, axis=1)+mean_Jacb2)[1:], (np.max(Jac_B2,axis=1)-mean_Jacb2)[1:]]
    fig = plt.figure(figsize=(8,6))
    plt.plot(psimesh, smooth_g*smooth_q+cpsi,color="tab:blue")
    plt.errorbar(psimesh[1:],mean_Jacb2[1:],var_Jacb2,color="tab:orange",ecolor="r",capsize=3,errorevery=5)
    plt.legend(["(gq+I)", "J*B^2"],fontsize=20)
    plt.xlabel(r"$\psi_p$")
    if figs=="save":
        plt.savefig("J_check.png")


    fig,ax = plt.subplots(dpi=120)
    ax.set_aspect('equal', 'box')
    # ax.plot(data["xlim"], data["ylim"],'k-',linewidth=0.5)
    psi_end = npsi+1
    for jt in range(0,ntheta+1,5):
        plt.plot(Ronptb[:psi_end,jt],Zonptb[:psi_end,jt],'k-')
    for jt in range(0,ntheta+1,5):
        plt.plot(Ronptp[:psi_end,jt],Zonptp[:psi_end,jt],'b-')
    for ip in range(0,psi_end,10):
        plt.plot(Ronptb[ip,:],Zonptb[ip,:],'k-')
    ax.plot(data["rbdry"], data["zbdry"],'r--')
    if figs=="save":
        plt.savefig("Boozer.png")

    fig,ax = plt.subplots(dpi=120)
    ax.set_aspect('equal', 'box')
    ax.pcolormesh(Ronptb,Zonptb,nu_onptb)
    ax.set_title(r"$\nu$")
    if figs=="save":
        plt.savefig("nu.png")
    
    fig,ax = plt.subplots(dpi=120)
    ax.set_aspect('equal', 'box')
    ax.pcolormesh(Ronptb,Zonptb,delta)
    ax.set_title(r"$\delta$")
    if figs=="save":
        plt.savefig("delta.png")

    if figs=="show":
        plt.show()

#    neg_q = None
#    if np.sign(gcurrent[0]) < 0:
#        while neg_q not in ("","y","n"):
#            neg_q = input("Negative g detected in this gfile and it is enforced positive in spdata. Do you want to make q profile negative instead? y/n default:n")
#        if neg_q == "y":
#            map_data["qpsi"] = -map_data["qpsi"]

    return map_data

def write_spdata(map_data, output):
    lsp = map_data["lsp"]
    lst = map_data["lst"]

    logging.info("construt bsp")
    bsp = np.zeros((9,lsp,lst), dtype="float32")
    bsp[0,:,:] = np.sqrt(map_data["BR_onptb"][:,:]**2 + map_data["BZ_onptb"][:,:]**2 + map_data["BT_onptb"][:,:]**2)

    spdpsi = map_data["psimesh"][1]-map_data["psimesh"][0]
    spdtheta = map_data["tb_mesh"][1]-map_data["tb_mesh"][0]
    bsp = spline_function_f.spline_function.construct_spline2d(spdpsi,spdtheta,bsp,1,2,lsp,lst)

    logging.info("construt xsp")
    xsp = np.zeros((9,lsp,lst), dtype="float32")
    xsp[0,:,:] = map_data["Ronptb"][:,:]
    # xsp[0,:,:] = Ronptp[:,:]
    xsp = spline_function_f.spline_function.construct_spline2d(spdpsi,spdtheta,xsp,1,2,lsp,lst)

    logging.info("construt zsp")
    zsp = np.zeros((9,lsp,lst), dtype="float32")
    zsp[0,:,:] = map_data["Zonptb"][:,:] - map_data["Z_ctr"]
    # zsp[0,:,:] = Zonptp[:,:]
    zsp = spline_function_f.spline_function.construct_spline2d(spdpsi,spdtheta,zsp,1,2,lsp,lst)

    logging.info("construct delta")
    delsp = np.zeros((9,lsp,lst), dtype="float32")
    delsp[0,:,:] = map_data["delta"][:,:]
    delsp = spline_function_f.spline_function.construct_spline2d(spdpsi,spdtheta,delsp,0,2,lsp,lst)

    logging.info("construct fsp")
    fsp = np.zeros((9,lsp,lst), dtype="float32")
    fsp[0,:,:] = map_data["nu_neg"][:,:]
    fsp = spline_function_f.spline_function.construct_spline2d(spdpsi,spdtheta,fsp,0,2,lsp,lst)

    logging.info("construt jacsp")
    jacsp = np.zeros((9,lsp,lst), dtype="float32")
    for jt in range(0,lst):
        jacsp[0,:,jt] = (map_data["gpsi"]* map_data["qpsi"] + map_data["cpsi"]) / bsp[0,:,jt]

    jacsp = spline_function_f.spline_function.construct_spline2d(spdpsi,spdtheta,jacsp,1,2,lsp,lst)

    logging.info("constuct rd2d")
    rd2d = np.zeros((9,lsp,lst), dtype="float32")
    for jt in range(0,lst):
        rd2d[0,:,jt] = map_data["cpsi"]
    rd2d = spline_function_f.spline_function.construct_spline2d(spdpsi,spdtheta,rd2d,0,2,lsp,lst)

    logging.info("construt qpsi")
    qpsi = np.zeros((3,lsp), dtype="float32")
    qpsi[0,:] = map_data["qpsi"]
    qpsi = spline_function_f.spline_function.construct_spline1d(spdpsi,qpsi,0,lsp)

    logging.info("construct gpsi")
    gpsi = np.zeros((3,lsp), dtype="float32")
    gpsi[0,:] = map_data["gpsi"]
    gpsi = spline_function_f.spline_function.construct_spline1d(spdpsi,gpsi,0,lsp)

    logging.info("construct ppsi")
    ppsi = np.zeros((3,lsp), dtype="float32")
    ppsi[0,:] = map_data["ppsi"]
    ppsi = spline_function_f.spline_function.construct_spline1d(spdpsi,ppsi,0,lsp)

    logging.info("construct rpsi")
    rpsi = np.zeros((3,lsp), dtype="float32")
    rpsi[0,:] = (map_data["Ronptb"][:,0]-map_data["R_ctr"]) / map_data["R_ctr"]
    rpsi = spline_function_f.spline_function.construct_spline1d(spdpsi,rpsi,1,lsp)

    logging.info("construct torpsi")
    torpsi = np.zeros((3,lsp), dtype="float32")
    for ip in range(1,lsp):
        torpsi[0, ip] = torpsi[0, ip-1] + qpsi[0, ip-1] * spdpsi + 0.5*qpsi[1, ip-1] * spdpsi**2\
                        + 1/3*qpsi[2,ip-1] * spdpsi**3

    torpsi = spline_function_f.spline_function.construct_spline1d(spdpsi,torpsi,0,lsp)

    # write toroidal psi for profile.dat
    with open("torpsispickle","wb") as fb:
        pickle.dump(map_data["psimesh"],fb)
        pickle.dump(torpsi,fb)

    # write to spdata

    from fortran import spdata_interface_f

    output_file = output
    logging.info(f"Writing spdata to {output_file}...")

    spdata_interface_f.spdataio.filename(output_file)
    spdata_interface_f.spdataio.allocate(lsp,lst)

    spdata_interface_f.spdataio.bsp = bsp
    spdata_interface_f.spdataio.xsp = xsp
    spdata_interface_f.spdataio.zsp = zsp
    spdata_interface_f.spdataio.gsp = jacsp
    spdata_interface_f.spdataio.qpsi = qpsi
    spdata_interface_f.spdataio.gpsi = gpsi
    spdata_interface_f.spdataio.rd = rd2d
    spdata_interface_f.spdataio.ppsi = ppsi
    spdata_interface_f.spdataio.rpsi = rpsi
    spdata_interface_f.spdataio.fsp = fsp
    spdata_interface_f.spdataio.delsp = delsp
    spdata_interface_f.spdataio.torpsi = torpsi
    spdata_interface_f.spdataio.psiw = map_data["psimesh"][-1]
    spdata_interface_f.spdataio.ped = map_data["psimesh"][-1]

    spdata_interface_f.spdataio.write()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ini", help="A file containing parameters."
                        "Not necessary if parameters inputs in command lines.")
    parser.add_argument("--input", help="Path of input gfile")
    parser.add_argument("--output", help="Path of output spdata file")
    parser.add_argument("--lsp", help="Grid# in psi direction", type=int)
    parser.add_argument("--lst", help="Grid# in theta direction. Must be an even number", type=int)
    parser.add_argument("--psimax",help="ratio of psimax to LCFS, default=0.995",type=float,default=0.995)
    parser.add_argument("--detail-print", help="Set to False if no debug info is wanted. default: False", default=False, type=bool)
    parser.add_argument("--figs",help="Set to 'show' for popup plots. Set to 'save' to save figures. default: show",default="show")
    args = parser.parse_args()

    if args.ini:
        import json,sys
        with open(args.ini, "r") as fj:
            paras = json.load(fj)
            filenames = paras["gfilename"]
            for ind, inputfile in enumerate(filenames):
                data = mapping(inputfile, paras["npsi"], paras["ntheta"],paras["figs"])
                write_spdata(data)
            sys.exit()

    paras = args
    if paras.detail_print:
        logging.basicConfig(level=logging.INFO)
    if not paras.input:
        filename = input("Please input path of gfile:\n")
    else:
        filename = paras.input

    if not paras.lsp:
        lsp = input("Please input npsi (=lsp in GTC and spdata), default: 91\n")
        if not lsp:
            lsp = 91
        else:
            lsp = int(lsp)
    else:
        lsp = paras.lsp

    if not paras.lst:
        lst = input("Please input ntheta (=lst in GTC, =lst+1 in spdata), default: 122\n")
        if not lst:
            lst = 122
        else:
            lst = int(lst)
    else:
        lst = paras.lst

    psimax=paras.psimax
    print(f"Ratio of last mapping flux surface to LCFS: {psimax}. Use --psimax=[value] to specify other value.")

    if not paras.output:
        output = input("Please enter the output file path, default: './spdata.dat'\n")
        if not output:
            output = "./spdata.dat"
    else:
        output = paras.output

    map_data = mapping(filename, lsp, lst, psimax, args.figs)
    write_spdata(map_data, output)
