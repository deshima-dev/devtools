# -*- coding utf-8 -*-
#
# functions.py: Define functions for merge script
#
# Author : Tetsutaro Ueda
# Created: 2017/11/02
#-------------------------------- PUBLIC ITEMS
__all__ = [
    'create_bintablehdu',
    'load_obsinst',
    'get_maskid_corresp'
    'calibrate_to_power',
    'convert_asciitime',
    'convert_timestamp'
]

#-------------------------------- IMPORT MODULES
#-------- Standard Modules
from datetime import datetime
#-------- Dependent Packages
import numpy as np
import scipy.interpolate
from astropy.io import fits

#-------------------------------- FUNCTIONS
#---------------- Create Binary Table HDU from 'hdu_dict'
def create_bintablehdu(hd):
    header = fits.Header()
    for (i, j) in zip(hd['hdr_vals'].items(), hd['hdr_coms'].items()):
        header[i[0]] = i[1], j[1]

    columns = [
        fits.Column(name=i[0], format=j[1], array=i[1], unit=k[1])
        for (i, j, k) in zip(
            hd['col_vals'].items(),
            hd['col_form'].items(),
            hd['col_unit'].items()
        )
    ]

    hdu = fits.BinTableHDU.from_columns(columns, header)

    for i in hd['hdr_coms'].items():
        hdu.header.comments[i[0]] = i[1]

    return hdu


#---------------- Get data for 'OBSINFO'
def load_obsinst(obsinst):
    if not '.obs' in obsinst:
        raise ValueError('The input file must be an observational instruction!!')

    with open(obsinst, 'r') as f:
        for line in f:
            if '% OBSERVER' in line:
                observer = line.rstrip().split('=')[1]                          # Get 'OBSERVER'
            elif '% SRC_NAME' in line:
                obs_object = line.rstrip().split('=')[1]                        # Get 'OBJECT'
            elif '% EPOCH' in line:
                equinox = line.rstrip().split('=')[1].strip('J').strip('B')     # Get 'EQUINOX'

    return observer, obs_object, equinox


#---------------- Get Correspondance of 'master' and 'kid'
def get_maskid_corresp(des, filt):
#-------- Correspondance between 'masterid' and 'attribute'
    mas_vs_attr = [[x, y] for (x, y) in zip(des['masterid'], des['attribute'])]
    del mas_vs_attr[13:17]
    mas_vs_attr.append([-1, np.nan])
#-------- Correspondance between 'masterid' and 'kidid'
    mas_vs_kid = [[x, y, z[0]] for (x, y, z) in zip(filt['masterid'], filt['kidid'], filt['F_filter, dF_filter'])]
#-------- Correspondance between 'kidid' and 'attribute'
    kid_vs_attr = [
        [x[0], x[1], y[1], x[2]]
         for (x, y) in zip(
            sorted(mas_vs_kid, key=lambda x: x[0]),
            sorted(mas_vs_attr, key=lambda x: x[0])
        )
    ]
    kid_vs_attr = sorted(kid_vs_attr, key=lambda x: x[1])   # Sorted by 'kidid'
    kid_vs_attr = list(map(list, zip(*kid_vs_attr)))        # Transpose
#-------- Get Data
    masterids = np.array([kid_vs_attr[0]])
    kidids    = np.array([kid_vs_attr[1]])
    kidattr   = kid_vs_attr[2]
    kidfreqs  = np.array([kid_vs_attr[3]])
#---- Get kidtypes
    kidtypes = []
    for x in kidattr:
        if x=='wideband':   attr = 0
        elif x=='filter':   attr = 1
        elif x=='blind':    attr = 2
        else:               attr = np.nan
        kidtypes.append(attr)

    return masterids, kidids, np.array([kidtypes]), kidfreqs


#---------------- Calibrate 'amplitude' and 'phase' to 'power'
#-------- ax^(1/2) + b
def sqrt_func(x, p0, p1):
    return p0*np.sqrt(x) + p1

#-------- Calibrate
def calibrate_to_power(pixelid, Troom, rhdus, ddb):
    nkid = rhdus['READOUT'].header['NKID%d' %pixelid]

    kiddict = {}
    for (i, j) in zip(ddb['KIDFILT'].data['kidid'], ddb['KIDFILT'].data['masterid']):
        kiddict[i] = j

    linphase = np.transpose([rhdus['READOUT'].data['Amp, Ph, linPh %d' %i].T[2] for i in range(nkid)])
    linyfc   = rhdus['KIDSINFO'].data['yfc, linyfc'].T[1]
    Qr       = rhdus['KIDSINFO'].data['Qr, dQr (300K)'].T[0]

    fshift = np.array((linphase - linyfc)/ (4.*Qr)).T

#---- Responsivity curve
    (p0, dp0, p1, dp1) = ddb['KIDRESP'].data['fit params'].T
    Tload = ddb['KIDRESP'].data['Raw Tload']
    Pload = ddb['KIDRESP'].data['Raw Pload']

    y = sqrt_func(Tload.T, p0, p1) - sqrt_func(Troom, p0, p1)
    y = y.T

#---- Interpolate
    Tsignal = []
    for i in range(nkid):
        masterid = kiddict[i]
        if masterid < 0:
            Tsignal.append([np.nan for i in range(len(fshift[i]))])
            continue

        tck = scipy.interpolate.splrep(y[i], Tload[i], s=0)
        Tsignal.append(scipy.interpolate.splev(fshift[i], tck, der=0))

#---- Convert to power
    Psignal = []
    for i in range(nkid):
        tck = scipy.interpolate.splrep(Tload[i], Pload[i], s=0)
        Psignal.append(scipy.interpolate.splev(Tsignal[i], tck, der=0))


    return np.array(Tsignal).T, np.array([Psignal]).T


#---------------- Convert times to FITS times
#-------- Ascii time
def convert_asciitime(asciitime, form_fitstime):
    asciitime = [datetime.strptime('%14.6f' %t, '%Y%m%d%H%M%S.%f') for t in asciitime]
    asciitime = [datetime.strftime(t, form_fitstime) for t in asciitime]

    return np.array(asciitime)

#-------- Timestamp
def convert_timestamp(timestamp):
    timestamp = [datetime.utcfromtimestamp(t) for t in timestamp]
    timestamp = [datetime.strftime(t, FORM_FITSTIME_P) for t in timestamp]

    return np.array(timestamp)
