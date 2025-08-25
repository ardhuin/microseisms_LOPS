import numpy as np
import os
import datetime
import netCDF4
from netCDF4 import Dataset
import xarray as xr
from scipy.interpolate import interp1d
from datetime import datetime, timedelta


def seismic_define_C2(readc, CgR):
    """
    Defines seismic response function (C2, nC2, c).
    
    Parameters
    ----------
    readc : int
        0 = analytic fit (Kedar et al. 2008 / Longuet-Higgins 1950)
        1 = read from Rayleigh_source.txt
    CgR : float
        Group speed. If < 0, read from Rayleigh_Cg.txt
    
    Returns
    -------
    C2 : ndarray
        Seismic response function (c^2)
    nC2 : int
        Number of discretized points
    c : ndarray
        
    """
    # --- Default version: analytic fit ---
    nC2 = 601
    x0 = np.linspace(0, 6, nC2)  # sig*h/beta, discretized
    
    a, b, c = 1.0, 10.0, 0.4
    c1 = 0.6 * (a + c * x0) / (a + b * (x0 - 0.85)**2) \
         + 0.28 / (2.0 + 1.0 * x0) \
         - 0.008 * (x0 - 2.5)

    a, b, c, d = 6.0, 30.0, 0.5, 10.0
    c2 = (0.30 * (a + c * x0) / (a + b * (x0 - 2.77)**2)
           + 0.40 / (1.0 + x0)
           + 0.006 * (x0 - 5)) * np.maximum(np.tanh((x0 - 1) * d), 0) \
         + 0.02 * np.exp(-(1.5 * (x0 - 3.6))**2) \
         - np.maximum(0.02 * (x0 - 5), 0)

    C2 = c1**2 + c2**2
    c0 = np.sqrt(c1**2)

    # --- Version 2: read from file ---
    if readc == 1:
        nmodes = 1
        all_c = np.loadtxt("Rayleigh_source.txt")
        nC, nmodein = all_c.shape[0], all_c.shape[1] - 1

        if nmodes > nmodein:
            nmodes = nmodein

        sighoverbetain = np.pi * all_c[:, 0]
        nC2 = 1001
        x0 = np.linspace(0, 10, nC2)
        C2 = np.zeros_like(x0)

        for i in range(nmodes):
            interp_func = interp1d(sighoverbetain, all_c[:, i+1],
                                   bounds_error=False, fill_value=np.nan)
            c = interp_func(x0)
            mask = np.isfinite(c)
            C2[mask] += c[mask]**2

        C2[0] = C2[1]
        c0 = np.sqrt(C2)

    # --- Group speed definition ---
    if CgR < 0:
        all_Cg = np.loadtxt("Rayleigh_Cg.txt")
        nCg, nmodein = all_Cg.shape[0], all_Cg.shape[1] - 1

        sighoverbetain = np.pi * all_Cg[:, 0]
        nC2 = 1001
        x0 = np.linspace(0, 10, nC2)

        # First mode (reference)
        interp_Cg = interp1d(sighoverbetain, all_Cg[:, 1],
                             bounds_error=False, fill_value=np.nan)
        CgRcmax = interp_Cg(x0)

        interp_c0 = interp1d(sighoverbetain, all_c[:, 1],
                             bounds_error=False, fill_value=np.nan)
        c0 = interp_c0(x0)

        # Loop over additional modes
        for i in range(1, nmodein):
            interp_c = interp1d(sighoverbetain, all_c[:, i+1],
                                bounds_error=False, fill_value=np.nan)
            interp_Cgi = interp1d(sighoverbetain, all_Cg[:, i+1],
                                  bounds_error=False, fill_value=np.nan)

            c = interp_c(x0)
            Cgi = interp_Cgi(x0)

            mask = np.isfinite(c) & np.isfinite(c0) & (c > c0)
            CgRcmax[mask] = Cgi[mask]
            c0[mask] = c[mask]

        # Fill NaNs with previous value
        isnan_mask = ~np.isfinite(CgRcmax)
        if np.any(isnan_mask):
            first_valid = np.argmax(np.isfinite(CgRcmax))
            CgRcmax[isnan_mask] = CgRcmax[first_valid]

    else:
        CgRcmax = np.zeros(nC2) + CgR

    return C2, nC2, c




def dispNewtonTH(f, dep):
    """
    Inverts the linear dispersion relation (2*pi*f)^2 = g*k*tanh(k*dep)
    to get k from f and dep.
    
    Parameters
    ----------
    f : array_like
        Frequency in Hz
    dep : float
        Water depth (m)
        
    Returns
    -------
    k : ndarray
        Wavenumber corresponding to frequency and depth
    """
    g = 9.81
    eps = 1e-6

    f = np.atleast_1d(f)
    sig = 2 * np.pi * f
    Y = dep * sig**2 / g   # squared dimensionless frequency

    X = np.sqrt(Y)         # initial guess that fits deep water
    
    # mask: only apply Newton iteration where needed
    mask = X < 8

    if np.any(mask):
        F = np.ones_like(X[mask])
        while np.max(np.abs(F)) > eps:
            H = np.tanh(X[mask])
            F = Y[mask] - X[mask] * H
            FD = -H - X[mask] / np.cosh(X[mask])**2
            X[mask] = X[mask] - F / FD

    return X / dep


def dist_sphere(lon1, lon2, lat1, lat2):
    # Haversine formula for spherical distance in degrees
    R = 6371.0
    lat1, lat2 = np.radians(lat1), np.radians(lat2)
    lon1, lon2 = np.radians(lon1), np.radians(lon2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    c = 2*np.arcsin(np.sqrt(a))
    return np.degrees(c)


###########################################################################
# --- Main function ---
def seismic_response_primary_sandwaves(
        wave_spectrum,
        bottom_topography_spectrum,
        dpt, lat, lon, CgR, Q, #date1, date2,
        lono, lato, dA, statname):
    
    
    g    = 9.81
    rhow = 1026
    rhos = 2600
    betas= 2800
    R_E  = 4e7/(2*np.pi)

    # --- read bottom topo spectrum ---
    with open(bottom_topography_spectrum, "r") as f:
        A = np.fromfile(f, sep=" ")

    # first two = integers, next two = reals
    nkbx = int(A[0])
    nkby = int(A[1])
    dkbx = float(A[2])
    dkby = float(A[3])

    nvals = nkbx * nkby
    flat = A[4:4 + nvals]
    if flat.size != nvals:
        raise ValueError(f"Expected {nvals} values for botspec, found {flat.size}")
    botspec = flat.reshape((nkby, nkbx), order="F").T  # -> shape (nkbx, nkby)
    

    kbxmax = ((nkbx - nkbx%2)/2.0) * dkbx
    kbymax = ((nkby - nkby%2)/2.0) * dkby
    kbotx = np.linspace(-kbxmax, kbxmax, nkbx)
    kboty = np.linspace(-kbymax - (1-nkby%2)*dkby, kbymax, nkby)

    # --- seismic coefficients ---
    C2, nC2, c = seismic_define_C2(0, CgR)
    factor1 = 2*np.pi * (1/rhos)**2 / (betas**5 * R_E)

    # --- time axis ---
    #nt = int((date2-date1)/dt) + 1
    #times = np.linspace(date1, date2, nt)

    # --- read wave spectrum file ---
    ds = Dataset(wave_spectrum, 'r')
    freqp = np.array(ds.variables['frequency'][:])
    theta = np.array(ds.variables['direction'][:])
    times = ds.variables['time'][:]
    lonp  = ds.variables['longitude'][0,0] # assumes fixed position
    latp  = ds.variables['latitude'][0,0]
    dptpall = ds.variables['dpt'][:]
    efthall = ds.variables['efth']
        
    nf = len(freqp)
    nd = len(theta)
    nt = len(times)
    
    xfr = np.exp(np.log(freqp[nf - 1] / freqp[0]) / (nf - 1))  # geometric progression factor
    df = freqp * 0.5 * (xfr - 1.0 / xfr)                       # frequency intervals in wave model


    
    omega = 2*np.pi*freqp

    alpha = dist_sphere(lono, lonp, lato, latp) * (np.pi/180)
    print('alpha:',alpha,xfr)
    dtor = np.pi/180
    k = dispNewtonTH(freqp, dpt)


# 5. Pre-computes Efth spectrum to pressure transformation 

    Eftop1f = np.zeros((nd, nf))
    botspeci = np.zeros((nd, nf))

    dtor = np.pi / 180.0

    for i in range(nf):
        for j in range(nd):
            # compute wave-number components in x and y
            kbx = k[i] * np.sin(theta[j] * dtor)
            kby = k[i] * np.cos(theta[j] * dtor)

            # wavenumber in grid units 
            kbotxi = (kbxmax + kbx) / dkbx
            kbotyi = (kbymax + kby) / dkby
            # clamp indices (Python 0-based vs MATLAB 1-based)
            ibk = max(min(int(np.floor(kbotxi)), nkbx - 2), 0)
            jbk = max(min(int(np.floor(kbotyi)), nkby - 2), 0)

            # fractional part for bilinear interpolation
            xbk = kbotxi - np.floor(kbotxi)
            ybk = kbotyi - np.floor(kbotyi)

            # bilinear interpolation of bottom spectrum onto wave wavenumber 
            botspeci[j, i] = (
                (botspec[ibk, jbk] * (1 - ybk) + botspec[ibk, jbk + 1] * ybk) * (1 - xbk)
                + (botspec[ibk + 1, jbk] * (1 - ybk) + botspec[ibk + 1, jbk + 1] * ybk) * xbk
            )
            #print('Testing bottom spectrum:',i,j,ibk,jbk,freqp[i],k[i],dkbx,dkby,kbxmax,kbx,kby,botspeci[j, i])



    # outputs
    F_delta = np.zeros((nt, nf))
    Fx_delta= np.zeros((nt, nf))
    Fy_delta= np.zeros((nt, nf))


    # intermediate transfer functions
    Eftop1f = np.zeros((nd,nf))
    alphas = np.zeros(nf)

    coeff = np.ones(nf)
    coeff_all = np.ones(nf)
    attenuation = np.ones(nf)
    coeff_love=np.ones(nf)
    coeff_all_love = np.ones(nf)
    attenuation_love=np.ones(nf)
        
    for i in range(nf):
            om = omega[i]
            k0 = k[i]
            h0 = dpt #dptpall[0,it]
            kh = k0*dpt

            if kh < 7:
                sinh2kh=np.sinh(2*kh)
                dkdD=-2*k0**2/(2*kh+ sinh2kh)
                dkDdD=k0*np.sinh(2*kh)/(2*kh+ sinh2kh)
                Cg=om/k0*0.5*(1+(2*kh)/sinh2kh)
                dCgdD=om/k0*(dkDdD/sinh2kh-2*kh*dkDdD*np.cosh(2*kh)/(sinh2kh**2))-dkdD*om/k0**2*0.5*(1+(2*kh)/sinh2kh)
                alpha1=-(k0+h0*dkdD)*np.tanh(k0*h0)/k0
                alpha2=-0.5*dCgdD/(Cg*k0)
                alpha3=-dkdD/(k0**2)
                alphas[i] = alpha1+alpha2+alpha3
            else:
                alphas[i] = -1.0

            Eftop1f[:,i] = botspeci[:,i] * (rhow*g*k0*alphas[i]/np.cosh(kh))**2 * (2*np.pi/nd)

            #print('TEST Eftop1f:',i,freqp[i],dpt,omega[i]**2,g*k0*np.tanh(kh),(rhow*g*k0*alphas[i]/np.cosh(kh))**2 * (2*np.pi/nd),k0,alphas[i],np.sum(botspeci[:,i]),np.sum(Eftop1f[:,i]) )

            omehoverbeta=0*omega[i]*dpt/betas;   # set to zero because interaction in shallow water
            # gets the nearest discretized omega*h/beta
            ind = (100.0 * omehoverbeta).astype(int)
            if ind > nC2-1: 
                ind= nC2-1
            coeff[i] = (factor1 * omega[i] * C2[ind]) / np.sin(alpha)
  
            # copy to coeff_love (placeholder)
            coeff_love[i] = coeff[i]  # to be corrected...
            
            
            # attenuation for one orbit
            b = np.exp(-omega[i] * (2.0 * np.pi) * (R_E / (abs(CgR) * Q[i])))

            # terms for shorter arc
            attenuation[i] = np.exp(-omega[i] * alpha * (R_E / (abs(CgR) * Q[i]))) / (1.0 - b)

            # terms for longer arc
            attenuation[i] += np.exp(-omega[i] * (2.0 * np.pi - alpha) * (R_E / (abs(CgR) * Q[i]))) / (1.0 - b)

            # copy to Love-wave attenuation (placeholder)
            attenuation_love[i] = attenuation[i]  # to be corrected...


            coeff_all[i] = coeff[i] * attenuation[i] * dA
            coeff_all_love[i] = coeff_love[i] * attenuation_love[i] * dA
            #print('TEST COEF:',i,C2[ind],coeff[i],'##',attenuation[i],'##',dA,'##',coeff_all[i])
    

    # main loop over time steps 
    for it in range(nt):
        #itf = int(round((times[it]-times[0])/dt))

        efth = efthall[it,0,:,:].T
        p1f = np.sum(Eftop1f*efth, axis=0)
        chk1 = np.sum(Eftop1f, axis=0)
        chk2 = np.sum(efth, axis=0)
        #for i in range(nf):
        #    print('I ...',it,i,p1f[i],chk1[i],chk2[i],coeff_all[i])
        t1xf= np.sum(Eftop1f*efth*(np.sin(theta*dtor)[:,None]**2)/(alphas[None,:]**2), axis=0)
        t1yf= np.sum(Eftop1f*efth*(np.cos(theta*dtor)[:,None]**2)/(alphas[None,:]**2), axis=0)

        source = coeff_all * p1f
        sourcex= coeff_all_love * t1xf
        sourcey= coeff_all_love * t1yf

        F_delta[it,:] = source
        Fx_delta[it,:]= sourcex
        Fy_delta[it,:]= sourcey
        

    return F_delta, Fx_delta, Fy_delta, freqp, df, times


##########################################################
def station_info(freq, statname):
    network = 'undefined'
    startdate = 'undefined'
    if statname == 'ADK':
        lono = -176.6844
        lato = 51.8837
        network = 'IU'
        startdate = '2009/07/19'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1E-8         
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'AIS':
        lono = 77.57
        lato = -37.80
        network = 'G'
        startdate = '1993/12/25'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1E-8         
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'ALOHA':
        lato = 22 + 45 / 60.0
        lono = -158.0
        network = 'UH'
        startdate = '2007/02/14'
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 6
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'ANMO':
        lato = 34.946
        lono = -106.457
        network = 'IU'
        startdate = '2007/02/14'
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 1.5
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'ATD':
        lono = 42.85
        lato = 11.53  # Arta Cave - Arta, Republic of Djibouti 
        network = 'G'
        startdate = '2010/01/01'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1E-8         
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'BBSR':
        lato = 32.3712
        lono = -64.6962
        network = 'IU'
        startdate = '2002/07/19'
        thresh1 = 1E-14  # SAME AS CMLA
        thresh2 = 10     # SAME AS CMLA
        thresh3 = 1E-9   # SAME AS CMLA
        maxd = 8        # SAME AS CMLA
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'BFO':
        station_full_name = 'Black Forest Observatory, Schiltach, Germany'  
        lato = 48.33
        lono = 8.33         
        network = 'II'
        startdate = '1986/01/01'
        thresh1 = 1E-15
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'BKS':
        lato = 37.8762
        lono = -122.2356
        network = 'BK'
        startdate = '1991/05/01'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        maxd = 5
        r = np.maximum(0, (1.4 - 8 * freq))
    elif statname == 'BORG':
        lato = 64.7474
        lono = -21.3268   
        network = 'IU'
        startdate = '1994/07/30'
        thresh1 = 3E-12
        thresh2 = 8
        thresh3 = 1.E-8
        maxd = 48
        r = np.maximum(0, (0.8 - 5 * freq))
        picture = 'http://www.iris.edu/hq/gallery/photo/443'
    elif statname == 'CAN':
        lono = 149.00
        lato = -35.32
        network = 'G'
        startdate = '1987/11/27'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1E-9
        maxd = 8
        r = np.maximum(0, (0.8 - 5 * freq))  # to be adjusted ... 
    elif statname == 'CMLA':
        lono = -25.5243
        lato = 37.7637
        network = 'II'
        startdate = '1996/03/10'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1E-9
        maxd = 8
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'CCM':
        lono = -91.245
        lato = 38.056
        network = 'IU'
        startdate = '1989/01/01'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1E-9
        maxd = 2
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'COR':
        lato = 44.586
        lono = -123.305
        network = 'IU'
        startdate = '1992/01/01'
        thresh1 = 1E-14
        thresh2 = 8
        thresh3 = 1.E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq)) / 6
    elif statname == 'COYC':
        lato = -45.57
        lono = -72.08
        network = 'G'
        startdate = '2004/12/17'
        thresh1 = 1E-14
        thresh2 = 8
        thresh3 = 1.E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq)) / 6
    elif statname == 'CRZF':
        lato = -46.43
        lono = 51.86
        network = 'G'
        startdate = '1986/02/01'
        thresh1 = 1E-13
        thresh2 = 4
        thresh3 = 1.E-8
        maxd = 10
    elif statname == 'DBG':
        station_full_name = 'Daneborg, Greenland'
        lato = 74.31
        lono = -20.22
        startdate = '2010/08/11'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        network = 'II'
        maxd = 5
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'DRV':
        lato = -66.66
        lono = 140.00
        network = 'G'
        startdate = '1986/02/01'
        thresh1 = 1E-14
        thresh2 = 8
        thresh3 = 1.E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq)) / 6
    elif statname == 'DWPF':
        lato = 28.11
        lono = -81.43
        network = 'IU'
        startdate = '2010/09/24'
        thresh1 = 1E-14
        thresh2 = 8
        thresh3 = 1.E-8
        maxd = 2
        r = np.maximum(0, (0.8 - 5 * freq)) / 6
    elif statname == 'EDA':
        lono = 10.153427
        lato = 3.778868
        network = 'G'
        startdate = '2019/05/04'
        thresh1 = 1E-14
        thresh2 = 8
        thresh3 = 1.E-8
        maxd = 2
        maxdp = 1
        r = np.maximum(0, (0.8 - 5 * freq))/6
    elif statname == 'ELSH':
        lono = 1.136600
        lato = 51.147600
        network = 'GB'
        startdate = '2008/06/21'
        thresh1 = 1E-15
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 10
        maxdp = 0.3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'ERM':
        lono = 143.1571
        lato = 42.0150
        network = 'IU'
        startdate = '1990/05/21'
        thresh1 = float('nan')
        thresh2 = float('nan')
        thresh3 = float('nan')
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'ESK':
        lato = 55.317
        lono = -3.205
        network = 'II'
        startdate = '1987/11/13'
        thresh1 = 1E-14
        thresh2 = 30
        thresh3 = 5.E-8
        maxd = 5
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'FDF':
        lato = 14.7349711
        lono = -61.146311
        network = 'G'
        startdate = '2008/09/01'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1.E-9
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'FOMA':
        lato = -24.98
        lono = 46.98
        network = 'G'
        startdate = '2008/09/01'
        thresh1 = 1E-14
        thresh2 = 30
        thresh3 = 5.E-8
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'GRFO':
        station_full_name = 'Grafenberg, Germany'
        lono = 49.69
        lato = 11.22
        network = 'IU'
        startdate = '1994/01/26'
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 5
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'H2O':
        lono = -141.991736
        lato = 27.881910
        network = 'H2'
        startdate = '1999/10/04'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'HRV':
        lato = 42.5064
        lono = -71.5583
        network = 'IU'
        startdate = '1988/01/01'
        thresh1 = 0.5E-14
        thresh2 = 3
        thresh3 = 1.E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'ICC':
        network = 'G'
        startdate = '1996/04/01'
        lato = -20.28
        lono = -70.03
        thresh1 = 1E-14
        thresh2 = 8
        thresh3 = 1.E-9
        maxd = 4
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'IDGL':
        lono = -7.51
        lato = 55.07
        network = 'EI'
        startdate = '2011/01/01'
        thresh1 = 1e-13
        thresh2 = 4
        thresh3 = 1E-11
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'IGLA':
        lono = -9.375000
        lato = 53.419500
        network = 'EI'
        startdate = '2011/01/01'
        thresh1 = 3e-15
        thresh2 = 4
        thresh3 = 1E-11
        maxd = 10
        maxdp = 0.6
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'INU':
        network = 'G'
        lato = 35.350
        lono = 137.029  # Inuyama, Japon
        thresh1 = 1E-14
        thresh2 = 8
        thresh3 = 1.E-9
        maxd = 4
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'JOHN':
        lono = -169.5292
        lato = 16.7329
        network = 'IU'
        startdate = '1998/08/04'
        thresh1 = 1E-13  # same as KIP
        thresh2 = 8  # same as KIP
        thresh3 = 1.E-11  # same as KIP
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'JTS':
        lono = -84.953
        lato = 10.291
        network = 'IU'
        thresh1 = 1E-15
        thresh2 = 10
        thresh3 = 5.E-8
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'KBS':
        station_full_name = 'Ny-Alesund, Spitzbergen, Norway'
        lono = 11.94
        lato = 78.92
        network = 'IU'
        startdate = '2010/06/08'
        thresh1 = 1E-13  # same as KIP
        thresh2 = 8  # same as KIP
        thresh3 = 1.E-11  # same as KIP
        maxd = 5
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'KEKH':
        lato = 21.98
        lono = -159.71
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 6
        r = np.maximum(0, (0.8 - 5 * freq))
        station_full_name = 'Kekaha, Kauai, Hawaii'
        network = 'PT'
    elif statname == 'KDAK':
        lato = 57.7828
        lono = -152.5835
        network = 'II'
        startdate = '1997/06/09'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'KIP':
        lato = 21.4233
        lono = -158.015
        network = 'G'
        startdate = '1986/04/17'
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 6
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'KIEV':
        lato = 50.70
        lono = 29.22
        network = 'IU'
        startdate = '1995/01/30'
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'KNTN':
        lato = -2.77
        lono = -171.72
        network = 'IU'
        startdate = '2007/12/04'
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 6
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'KONO':
        lato = 59.6491
        lono = 9.5982
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'KWAJ':
        lono = 167.6130
        lato = 8.8019
        thresh1 = 1e-13  # SAME AS WAKE
        thresh2 = 4  # SAME AS WAKE
        thresh3 = 1E-11  # SAME AS WAKE
        maxd = 10  # SAME AS WAKE
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'LCO':
        lato = -29.01
        lono = -70.70  # Las Campanas Astronomical Observatory, Chile
        network = 'IU'
        startdate = '2009/07/30'
        thresh1 = 1E-13  # 'to be defined';
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 6
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'LVC':
        lato = -22.61
        lono = -68.91
        network = 'IU'
        startdate = '2009/04/10'
        thresh1 = 1E-13  # 'to be defined';
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 2
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'MAJO':
        lono = 138.2070
        lato = 36.5427
        thresh1 = float('nan')
        thresh2 = float('nan')
        thresh3 = float('nan')
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'MBO':
        station_full_name = 'M Bour, Senegal'
        lato = 14.39
        lono = -16.96
        network = 'G'
        startdate = '1985/09/01'
        thresh1 = 1E-13  # 'to be defined';
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 0.5
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'MBWA':
        station_full_name = 'Marble Bar, Western Australia'
        lato = -21.159
        lono = 119.731
        network = 'IU'
        startdate = '2001/09/01'
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 4
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'MIDW':
        lono = -177.3697
        lato = 28.2157
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 6
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'NACB':
        lato = -24.1738
        lono = 121.5947
        network = 'TW'
        thresh1 = 3e-15
        thresh2 = 4
        thresh3 = 1E-11
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'NE71':
        lato = 31.68973
        lono = -115.90526
        thresh1 = 1e-14
        thresh2 = 5
        thresh3 = 1
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'NE75':
        lato = 27.29334
        lono = -112.85649
        thresh1 = 1E-14
        thresh2 = 5
        thresh3 = 1
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'NE79':
        lato = 23.11937
        lono = -109.75611
        thresh1 = 1E-14
        thresh2 = 5
        thresh3 = 1
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'NNA':
        lato = -11.988
        lono = -76.842
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        network = 'II'
        maxd = 1
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'NOR':
        station_full_name = 'Station Nord, Greenland'
        lato = 81.60
        lono = -16.66
        startdate = '2010/07/06'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        network = 'II'
        maxd = 2
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'NWAO':
        station_full_name = 'Narrogin, Australia'
        lato = -32.93
        lono = 117.24
        network = 'IU'
        startdate = '1991/11/25'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        maxd = 4
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'OTAV':
        lato = 0.24
        lono = -78.45
        network = 'IU'
        startdate = '2009/04/01'
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1.E-9
        maxd = 0.5
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'PFO':
        lato = 33.6092
        lono = -116.4553
        thresh1 = 1E-13
        thresh2 = 4
        thresh3 = 1.E-8
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'POHA':
        lato = 19.7575
        lono = -155.5325
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 9
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'PAB':
        lato = 39.5458
        lono = -4.3483
        thresh1 = 0.5E-14
        thresh2 = 10
        thresh3 = 1.E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'PAF':
        lato = -49.35
        lono = 70.21
        network = 'G'
        thresh1 = 1E-13
        thresh2 = 4
        thresh3 = 1.E-8
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'PAS':
        lato = 34.15
        lono = -118.15
        thresh1 = 1E-13
        thresh2 = 4
        thresh3 = 1.E-8
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'PAYG':
        network = 'IU'
        startdate = '2010/05/05'
        lato = -0.67
        lono = -90.29
        thresh1 = 1E-14
        thresh2 = 8
        thresh3 = 1.E-9
        maxd = 4
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'PEL':
        network = 'G'
        startdate = '1995/10/04'
        lato = -33.14
        lono = -70.67  # Peldehue, Chile
        thresh1 = 1E-14  # TO BE ADJUSTED
        thresh2 = 8
        thresh3 = 1.E-9
        maxd = 2
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'PET':
        lono = 158.6531
        lato = 53.0239
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1E-8
        maxd = 4
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'PTCN':
        network = 'IU'
        startdate = '1996/12/29'
        lato = -25.07
        lono = -130.10  # Pitcairn Island, South Pacific
        thresh1 = 1E-14  # TO BE ADJUSTED
        thresh2 = 8
        thresh3 = 1.E-9
        maxd = 4
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'PPT':
        lato = -17.569
        lono = -149.576
        network = 'G'
        startdate = '1986/05/31'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq)) * 0.2
    elif statname == 'PPTF':
        lato = -17.569
        lono = -149.576
        network = 'G'
        startdate = '2009/05/27'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        maxd = 10
        r = np.maximum(0, (1.4 - 8 * freq))
    elif statname == 'NOR':
        station_full_name = 'Station Nord, Greenland'
        lato = 81.60
        lono = -16.66
        startdate = '2010/07/06'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        network = 'II'
        maxd = 2
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'PVAQ':
        lato = 37.40
        lono = -7.72
        network = 'PM'
        startdate = '2006/12/23'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1.E-9
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'RAR':
        lato = -21.21
        lono = -159.77
        network = 'IU'
        startdate = '1992/03/07'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1.E-9
        maxd = 30
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'RER':
        station_full_name = 'Riviere de l Est - Sainte Rose - La Reunion island, France'
        lato = -21.17
        lono = 55.74
        network = 'G'
        startdate = '1986/02/10'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1.E-9
        maxd = 6
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'ROCAM':
        lato = -19.76
        lono = 63.37
        network = 'G'
        startdate = '2012/12/15'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1.E-9
        maxd = 50
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'RODM':
        lato = -19.70  # Rodrigues Island, Republic of Mauritius
        lono = 63.44
        network = 'G'
        startdate = '2012/12/15'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1.E-9
        maxd = 50
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'ROSA':
        lato = 38.72
        lono = -28.25
        network = 'PM'
        startdate = '2008/03/06'
        thresh1 = 1E-15
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'RPN':
        lato = -27.13
        lono = -109.33  # Rapanui, Easter Island, Chile
        network = 'II'
        thresh1 = 1E-13  # TO BE ADJUSTED...
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 1.5
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'RSSD':
        lato = 44.121
        lono = -104.036
        network = 'IU'
        startdate = '1993/01/01'
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 1.5
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'SACV':
        station_full_name = 'Santiago Island, Cape Verde'
        lato = 14.97
        lono = -23.61
        network = 'II'
        startdate = '2000/05/29'
        thresh1 = 1E-15
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'SCO':  # FAIT FA 08/12/2013
        station_full_name = 'Ittoqqortoormiit, Greenland'
        lato = 70.49
        lono = -21.95
        startdate = '2010/08/05'
        thresh1 = 1E-12
        thresh2 = 10
        thresh3 = 1.E-6
        network = 'II'
        maxd = 5
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'SCZ':
        lato = 36.59798
        lono = -121.40481
        thresh1 = 1E-14
        thresh2 = 8
        thresh3 = 1.E-8
        maxd = 5
        network = 'G'
        r = np.maximum(0, (0.8 - 5 * freq))

    elif statname == 'SFJD':  # FAIT MO 09/06/11
        lono = -50.6215
        lato = 66.9960
        thresh1 = 5E-15
        thresh2 = 10
        thresh3 = 1.E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'SOEG':
        station_full_name = 'Sodalen, Greenland'
        lato = 68.20
        lono = -31.38
        startdate = '2011/07/27'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        network = 'II'
        maxd = 5
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'SPB':
        station_full_name = 'Sao Paulo, Brazil'
        lato = -23.59
        lono = -47.43
        network = 'G'
        startdate = '1996/06/17'
        thresh1 = 1E-15
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 5
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'SSB':
        lato = 45.279
        lono = 4.542
        network = 'G'
        startdate = '1982/05/02'
        thresh1 = 1E-15
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 3
        maxdp = 0.2
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'SSPA':  # FAIT MO 09/06/11
        lato = 40.6358
        lono = -77.8880
        thresh1 = 0.5E-14
        thresh2 = 10
        thresh3 = 1.E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'SUR':
        lato = -32.38
        lono = 20.81
        network = 'II'
        startdate = '1990/10/30'
        thresh1 = 1.E-13
        thresh2 = 8
        thresh3 = 1.E-8
        maxd = 3
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'PET':
        lono = 158.6531
        lato = 53.0239
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'TAM':
        lono = 5.53  # Tamanrasset, Algeria
        network = 'G'
        lato = 22.79
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1E-8
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
        startdate = '1983/11/16'
        thresh1 = 1E-15
        thresh2 = 10
        thresh3 = 5.E-9
        maxd = 0.5
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'TAOE':
        lato = -8.85
        lono = -140.15
        station_full_name = 'Taiohae - Marquesas islands, France'
        network = 'G'
        startdate = '2004/11/01'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        r = np.maximum(0, (0.8 - 5 * freq))
        maxd = 10
        
    elif statname == 'TEIG':
        lato = 20.23
        lono = -88.28  # Tamanrasset, Algeria
        network = 'IU'
        startdate = '2009/09/18'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        maxd = 1
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'TRIS':
        lato = -37.07
        lono = -12.32
        network = 'G'
        startdate = '2004/03/03'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        maxd = 4
        r = np.maximum(0, (0.8 - 5 * freq))
    elif statname == 'TRQA':
        station_full_name = 'Tornquist, Argentina'
        lato = -38.06
        lono = -61.98
        network = 'IU'
        startdate = '2000/10/28'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        maxd = 4
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'UCC':
        lato = 50.279
        lono = 1.5
        network = 'G'
        startdate = '1982/05/02'
        thresh1 = 1E-15
        thresh2 = 8
        thresh3 = 1.E-11
        maxd = 5
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'UNM':
        lato = 19.329662
        lono = -99.178065
        network = 'G'
        startdate = '1990/06/06'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 5.E-8
        maxd = 4
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'WAKE':
        lono = 166.6536
        lato = 19.2833
        thresh1 = 1e-13
        thresh2 = 4
        thresh3 = 1E-11
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
        
    elif statname == 'WKER1':
        thresh1 = 1e-13
        thresh2 = 4
        thresh3 = 1E-11
        maxd = 10
        network = 'IUEM'
        r = np.maximum(0, (0.8 - 5 * freq))
        lato = -46.634166
        lono = 60.1303
        
    elif statname == 'XMAS':
        thresh1 = 1e-13
        thresh2 = 4
        thresh3 = 1E-11
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
        lono = -153.910
        lato = 0.000
        
    elif statname == 'YSS':
        lono = 142.7550
        lato = 46.9539
        thresh1 = np.nan
        thresh2 = np.nan
        thresh3 = np.nan
        maxd = 10
        r = np.maximum(0, (0.8 - 5 * freq))
        
    else:
        print('STATION NOT FOUND')
        lato = np.nan
        lono = np.nan
        thresh1 = np.nan
        thresh2 = np.nan
        thresh3 = np.nan
        maxd = np.nan
        r = np.nan * np.ones(len(freq))
            
    return lato, lono, thresh1, thresh2, thresh3, maxd, r, network, startdate


############################################

import numpy as np

def microseisms_filterseisms(date0, frq0, spectre0, thresh1, thresh2, thresh3):
    # SUMMARY: removes seisms from seismic data
    # The following criteria are used:
    #    - (the relative change from one time step to the other at low
    #    frequencies should be less than thresh2 OR the ground displacement variance at low frequency 
    #    should be less than thresh1) AND the ground displacement variance at low frequency 
    #    should be less than thresh3
    
    nts = len(date0)
    frq2 = np.tile(frq0, (nts, 1))  # replicate frq0 along axis 0 (rows)
    df0 = frq0[1] - frq0[0]  # frequency step
    
    # Checks low frequency data
    Indexf = np.where((frq0 > 0.04) & (frq0 < 0.08))[0]
    Esum = np.sum(spectre0[:, Indexf], axis=1) * df0
    Er1 = np.ones(nts)
    Er2 = np.ones(nts)
    Er3 = np.ones(nts)
    
    # Calculate Er1
    Er1[1:nts-1] = np.exp(np.abs(np.log(Esum[1:nts-1]**2 / (Esum[0:nts-2] * Esum[2:nts]))))
    
    # Calculate Er2 and Er3
    Er2[0:nts-1] = np.exp(np.abs(np.log(Esum[0:nts-1] / Esum[1:nts]  ))) 
    Er3[0:nts-1] = np.exp(np.abs(np.log(Esum[1:nts]   / Esum[0:nts-1]))) 
    
    Er = Er1
    # Adjust Er based on conditions for finite values
    I = np.where(np.isfinite(Er1) == 0 & np.isfinite(Er2) == 1)[0]
    Er[I] = Er2[I]
    I = np.where(np.isfinite(Er) == 0 & np.isfinite(Er3) == 1)[0]
    Er[I] = Er3[I]
    I = np.where(np.isfinite(Er1) == 1 & np.isfinite(Er2) == 1)[0]
    Er[I] = np.maximum(Er2[I], Er[I])
    I = np.where(np.isfinite(Er1) == 1 & np.isfinite(Er3) == 1)[0]
    Er[I] = np.maximum(Er3[I], Er[I])
    
    # Shift Esum to also remove abnormally low noise (e.g., for KIP)
    Esum1 = Esum[0:nts-2]
    Esum2 = Esum[2:nts]
    Esum3 = Esum[1:nts-1]
    Ind1 = np.where(np.isfinite(Esum3) & np.isfinite(Esum1))[0]
    Esum[Ind1+1] = np.maximum(Esum1[Ind1], Esum3[Ind1])
    Ind1 = np.where(np.isfinite(Esum2) & np.isfinite(Esum3))[0]
    Esum[Ind1+1] = np.maximum(Esum1[Ind1], Esum3[Ind1])
    
    # Detection of seisms: either a very strong signal: Esum > thresh3
    Indbad = np.where((Er > np.abs(thresh2) & Esum > thresh1) | Esum > thresh3)[0]
    
    # Remove detected seisms from the spectrum
    spectre1 = spectre0.copy()
    spectre1[Indbad, :] = np.nan
    Esumout = np.sum(spectre1[:, Indexf], axis=1) * df0
    Erout = Er
    
    # Checks high frequency
    Indexf = np.where((frq0 > 0.4) & (frq0 < 0.48))[0]
    Esum = np.sum(spectre0[:, Indexf], axis=1) * df0
    Er1 = np.ones(nts)
    Er2 = np.ones(nts)
    Er3 = np.ones(nts)
    
    # Calculate Er1 for high frequency range
    Er1[1:nts-1] = Esum[1:nts-1]**2 / (Esum[0:nts-2] * Esum[2:nts])
    
    # Calculate Er2 and Er3 for high frequency range
    Er2[0:nts-1] = np.exp(np.abs(np.log(Esum[0:nts-1] / Esum[1:nts])))
    Er3[0:nts-1] = np.exp(np.abs(np.log(Esum[1:nts] / Esum[0:nts-1])))
    
    Er = Er1
    # Adjust Er based on conditions for finite values
    I = np.where(np.isfinite(Er1) == 0 & np.isfinite(Er2) == 1)[0]
    Er[I] = Er2[I]
    I = np.where(np.isfinite(Er) == 0 & np.isfinite(Er3) == 1)[0]
    Er[I] = Er3[I]
    I = np.where(np.isfinite(Er1) == 1 & np.isfinite(Er2) == 1)[0]
    Er[I] = np.maximum(Er2[I], Er[I])
    I = np.where(np.isfinite(Er1) == 1 & np.isfinite(Er3) == 1)[0]
    Er[I] = np.maximum(Er3[I], Er[I])
    
    # Detection of seisms: either a very strong signal: Esum > thresh3
    Indbad = np.where((Er > np.abs(thresh2) & Esum > 1E-5 * thresh1) |
                      Esum > thresh3 | np.isfinite(Er) == 0)[0]
    
    # Remove detected seisms from the spectrum
    spectre1[Indbad, :] = np.nan
    
    # Return the results
    return spectre1, Erout, Er2, Er3, Esumout


def read_spectrograms(mydatapath, statname, iydata):
    ynamed = f'{iydata:04d}'
    seismicdatafile = os.path.join(mydatapath, str(iydata), statname, f'{statname}.LHZ.{ynamed}_03h_depl.mat')
    seismicncfile = os.path.join(mydatapath, str(iydata), statname, f'{statname[:3]}-LHZ-03h_{ynamed}.nc')
    
    if os.path.exists(seismicncfile):
       S = netCDF4.Dataset(seismicncfile, 'r')
       spectre0 = np.squeeze(np.ma.getdata(S.variables['edzf'][:]) ).T
       spectre0_qc = np.squeeze(np.ma.getdata(S.variables['edzf_qc_level'][:])).T
       frequency = np.ma.getdata(S.variables['frequency'][:])
       date0 = np.ma.getdata(S.variables['time'][:]) 
    else:
        date0 = np.linspace(date1.toordinal(), date2.toordinal(), int((date2 - date1).days * 8))
        spectre0 = np.squeeze(np.zeros((len(date0), 48))) + 1E-18
        spectre0_qc = np.zeros_like(spectre0)
        
    return date0, frequency, spectre0, spectre0_qc



def read_modeled_spectrograms(file):
    S = netCDF4.Dataset(file, 'r')
    # This function will read the model results from the netCDF files
    Ef = np.squeeze(np.ma.getdata(S.variables['edzf'][:]))
    freqs = np.ma.getdata(S.variables['frequency'][:])
    f1 = np.ma.getdata(S.variables['frequency1'][:])
    f2 = np.ma.getdata(S.variables['frequency2'][:])
    dfs = f2-f1
    times = np.ma.getdata(S.variables['time'][:]) 
    return Ef, freqs, dfs, times


# Compute delta (seismic displacement)
def compute_delta_model(Ef_noref, Ef_ref, dfs, rs, ifmin, ifmax):
    [nt,nf]=Ef_ref.shape
    df2 = np.tile(dfs, (Ef_noref.shape[0], 1))
    print('df2:',df2.shape)
    r2 = np.tile(rs, (Ef_noref.shape[0], 1))
    delta_noref = 1E6 * np.sqrt(np.sum(Ef_noref[:,ifmin:ifmax] * df2[:,ifmin:ifmax], axis=1))
    delta_ref = 1E6 * np.sqrt(np.sum((Ef_noref[:,ifmin:ifmax] + r2[:,ifmin:ifmax] * (Ef_ref[:,ifmin:ifmax] - Ef_noref[:,ifmin:ifmax])) * df2[:,ifmin:ifmax], axis=1))
    return delta_noref, delta_ref
    
    
    # Compute delta (seismic displacement)
def compute_delta_obs(spectre0, frq0, freqs,ifmin, ifmax,date0,times,smooth=1,percentile=50):
# Handling dates and frequency data
    datesi = date0  # Assuming date0 is already in the correct format
    nts = len(date0)
    frq2 = np.tile(frq0, (nts, 1))  # Equivalent to `repmat` in MATLAB
    dfsi = frq0[1] - frq0[0]

    # Index for frequencies within the given range
    Indexf = np.where((frq0 > freqs[ifmin] * 0.95) & (frq0 < freqs[ifmax-1] * 1.05))[0]

    # Summing over the specified frequency range
    Esi = np.sum(spectre0[:, Indexf], axis=1) * dfsi

    # Filtering out based on QC values
    #II = np.where(spectre0_qc < 4)
    #spectre1[II] = np.nan

    # Summing for spectre1
    Esi1 = np.sum(spectre0[:, Indexf], axis=1) * dfsi
  
    # Calculating dsi1
    dsi1 = np.sqrt(Esi1) * 1e6

# Create an interpolation function
    #interpolate = interp1d(datesi, dsi1, kind='linear', fill_value="extrapolate")
    #delta_obs=interpolate(times)
    delta_obs=np.zeros(len(times))
    dt=times[1]-times[0]
    for it in range(len(times)):
        inds=np.where((datesi >= times[it]-dt*0.5*smooth) & (datesi < times[it]+dt*0.5*smooth))[0]
        if len(inds) > 0: 
           if percentile==50:
              delta_obs[it]=np.nanmedian(dsi1[inds])
           else: 
              delta_obs[it]=np.nanmin(dsi1[inds])
    return delta_obs
