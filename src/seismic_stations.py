import numpy as np

def seismic_stations(freq, statname):
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
        r = np.maximum(0, 0.8 - 5 * np.array(freq))

    elif statname == 'AIS':
        lono = 77.57
        lato = -37.80
        network = 'G'
        startdate = '1993/12/25'
        thresh1 = 1E-14
        thresh2 = 10
        thresh3 = 1E-8
        maxd = 10
        r = np.maximum(0, 0.8 - 5 * np.array(freq))

    elif statname == 'ALOHA':
        lato = 22 + 45 / 60.0
        lono = -158.0
        network = 'UH'
        startdate = '2007/02/14'
        thresh1 = 1E-13
        thresh2 = 8
        thresh3 = 1E-11
        maxd = 6
        r = np.maximum(0, 0.8 - 5 * np.array(freq))
    elif statname == 'EDA' 
        lato=	3.778868 	
        lono=	10.153427
        network='G'
        startdate='2019/05/04'
        thresh1=1E-14
        thresh2=8
        thresh3=1.E-8;
        maxd=2;
        maxdp=1;           
        r=max(0, (0.8-5*freq))./6;

    else:
        print('STATION NOT FOUND')
        lato = np.nan
        lono = np.nan
        thresh1 = np.nan
        thresh2 = np.nan
        thresh3 = np.nan
        maxd = np.nan
        r = np.full_like(freq, np.nan, dtype=float)

    maxdp = maxd  # Assuming maxdp = maxd like in MATLAB

    return lato, lono, thresh1, thresh2, thresh3, maxd, maxdp, r, network, startdate

