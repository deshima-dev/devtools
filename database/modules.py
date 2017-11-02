import os
from pathlib import Path

import numpy as np
from datetime import datetime
from astropy.io import fits
from astropy.io import ascii
from astropy import table
from astropy import coordinates
from astropy import units as u

D_ASTE    = (10.0* u.m).value                            # diameter of the ASTE
LON_ASTE  = coordinates.Angle('-67d42m11.89525s').deg    # longitude of the ASTE
LAT_ASTE  = coordinates.Angle('-22d58m17.69447s').deg    # latitude of the ASTE
PATH_DDBDICT = Path('/Users/tetsu/DESHIMA/devtools/merge/DB/DDB_dict.yaml').expanduser()
FORM_FITSTIME  = '%Y-%m-%dT%H:%M:%S'

def convert_fitstime(time_lis):
    if len(time_lis[0]) <= 14:
        return [datetime.strftime(datetime.strptime(time_lis[i], '%Y%m%d%H%M%S'),
                                  FORM_FITSTIME)
                for i in range(len(time_lis))]
    else:
        return [datetime.strftime(datetime.strptime(time_lis[i], '%Y%m%d%H%M%S.%f'),
                                  FORM_FITSTIME + '.%f')
                for i in range(len(time_lis))]
    

def addHeaderComments(hdr, key, com):
    hdr.comments[key] = com


def createBinTableHdu(data_dict):
    header = fits.Header()
    for (i, j, k) in zip(data_dict['hdr_key_lis'], data_dict['hdr_val_lis'], data_dict['hdr_com_lis']):
        header[i] = j, k

    columns = [fits.Column(name=data_dict['cols_key_lis'][i],
                           format=data_dict['tform'][i],
                           array=data_dict['cols_data_lis'][i],
                           unit=data_dict['tunit'][i])
               for i in range(len(data_dict['cols_key_lis']))]

    hdu = fits.BinTableHDU.from_columns(columns, header)
    [addHeaderComments(hdu.header, list(hdu.header)[i], data_dict['hdr_com_lis'][i])
     for i in range(- (len(data_dict['hdr_com_lis']) - data_dict['hdr_com_lis'].index('label for field 0')), 0)]

    return hdu

##### dict for KID design
def design_dict():
    hdr_key_lis = ['EXTNAME', 'FILENAME', ]
    hdr_val_lis = ['KIDDES', None, ]
    hdr_com_lis = ['name of binary data',
                   'input filename',
                   'label for field 0', 'data format of field 0',
                   'label for field 1', 'data format of field 1',
                   'label for field 2', 'data format of field 2',
                   'label for field 3', 'data format of field 3', 'data unit of field 3',
                   'label for field 4', 'data format of field 4', 'data unit of field 4']
    cols_key_lis = ['pixelid', 'masterid', 'attribute', 'fr_design', 'F_filter_design']
    cols_data_lis = []
    tform = ['I', 'I', '8A', 'E', 'E']
    tunit = [None, None, None, 'GHz', 'GHz']

    d_dict = {'hdr_key_lis': hdr_key_lis,
              'hdr_val_lis': hdr_val_lis,
              'hdr_com_lis': hdr_com_lis,
              'cols_key_lis': cols_key_lis,
              'cols_data_lis': cols_data_lis,
              'tform': tform,
              'tunit': tunit}
    return d_dict

##### dict for measured filter
def filter_dict():
    hdr_key_lis = ['EXTNAME', 'FILENAME', 'JSONNAME',]
    hdr_val_lis = ['KIDFILT', None, None,]
    hdr_com_lis = ['name of binary data',
                   'input data filename',
                   'input json filename',
                   'label for field 0', 'data format of field 0',
                   'label for field 1', 'data format of field 1',
                   'label for field 2', 'data format of field 2',
                   'label for field 3', 'data format of field 3',
                   'label for field 4', 'data format of field 4', 'data unit of field 4',
                   'label for field 5', 'data format of field 5',
                   '(A, dA, F0, dF0, W, dW) for Lorentzian', 'data format of field 6',
                   'label for field 7', 'data format of field 7', 'data unit of field 7',
                   'label for field 8', 'data format of field 8']
    cols_key_lis = ['pixelid', 'kidid', 'masterid', 'runid, framelen, refid',
                    'F_filter, dF_filter', 'Q_filter, dQ_filter', 'fit params', 'Raw Toptica F', 'Raw df resp.']
    cols_data_lis = []
    tform = ['I', 'I', 'I', '3I',
             '2E', '2E', '6E', None, None]
    tunit = [None, None, None, None,
             'GHz', None, None, 'GHz', None]

    f_dict = {'hdr_key_lis': hdr_key_lis,
              'hdr_val_lis': hdr_val_lis,
              'hdr_com_lis': hdr_com_lis,
              'cols_key_lis': cols_key_lis,
              'cols_data_lis': cols_data_lis,
              'tform': tform,
              'tunit': tunit}
    return f_dict

##### dict for responsivity
def resp_dict():
    hdr_key_lis = ['EXTNAME', 'FILENAME', 'JSONNAME',]
    hdr_val_lis = ['KIDRESP', None, None,]
    hdr_com_lis = ['name of binary data',
                   'input data filename',
                   'input json filename',
                   'label for field 0', 'data format of field 0',
                   'label for field 1', 'data format of field 1',
                   'label for field 2', 'data format of field 2',
                   'label for field 3', 'data format of field 3',
                   '(p0, dp0, p1, dp1): p0*sqrt(x) + p1', 'data format of field 4', 
                   'label for field 5', 'data format of field 5', 'data unit of field 5',
                   'label for field 6', 'data format of field 6', 'data unit of field 6',
                   'label for field 7', 'data format of field 7',
                   'label for field 8', 'data format of field 8']
    cols_key_lis = ['pixelid', 'kidid', 'masterid', 'runid, framelen',
                    'fit params', 'Raw Tload', 'Raw Pload', 'Raw df resp.', 'Raw dferr']
    cols_data_lis = []
    tform = ['I', 'I', 'I', '2I',
             '4E', None, None, None, None]
    tunit = [None, None, None, None,
             None, 'K', 'W', None, None]

    f_dict = {'hdr_key_lis': hdr_key_lis,
              'hdr_val_lis': hdr_val_lis,
              'hdr_com_lis': hdr_com_lis,
              'cols_key_lis': cols_key_lis,
              'cols_data_lis': cols_data_lis,
              'tform': tform,
              'tunit': tunit}
    return f_dict


##### dict for rcp
def rcp_dict():
    hdr_key_lis = ['EXTNAME', 'FILENAME']
    hdr_val_lis = ['KIDRCP', None]
    hdr_com_lis = ['name of binary table', 'filename??',
                   'label for field 0', 'data format of field 0',
                   'label for field 1', 'data format of field 1',
                   'label for field 2', 'data format of field 2', 'data unit of field 2',
                   'label for field 3', 'data format of field 3', 'data unit of field 3']
    cols_key_lis = ['pixelid', 'offsetaz', 'offsetel', 'gain']
    cols_data_lis = None
    tform = ['K', 'D', 'D', 'D']
    tunit = [None, 'deg', 'deg', '1']

    r_dict = {'hdr_key_lis': hdr_key_lis,
              'hdr_val_lis': hdr_val_lis,
              'hdr_com_lis': hdr_com_lis,
              'cols_key_lis': cols_key_lis,
              'cols_data_lis': cols_data_lis,
              'tform': tform,
              'tunit': tunit}

    return r_dict

##### calibration function
def calibrate_to_power(pixelid, Troom, rhdus, ddb):
    nkid = rhdus['READOUT'].header['NKID%d' %pixelid]
    
    kiddict = {}
    for i,j in zip(ddb['KIDFILT'].data['kidid'],ddb['KIDFILT'].data['masterid']):
        kiddict[i] = j
    #print( len(kiddict), kiddict )

    linphase = np.transpose( [rhdus['READOUT'].data['Amp, Ph, linPh %d' %i].T[2] for i in range(nkid)] )
    linyfc   = rhdus['KIDSINFO'].data['yfc, linyfc'].T[1]
    Qr       = rhdus['KIDSINFO'].data['Qr, dQr (300K)'].T[0]
    
    fshift = np.array( (linphase-linyfc)/4./Qr ).T
    #print(np.shape(fshift), fshift[0])
        
    ##### responsivity curve
    (p0, dp0, p1, dp1) = ddb['KIDRESP'].data['fit params'].T
    Tload = ddb['KIDRESP'].data['Raw Tload']
    Pload = ddb['KIDRESP'].data['Raw Pload']
    
    def func(x, p0, p1):
        return p0*np.sqrt(x) + p1

    y = func(Tload.T, p0, p1) - func(Troom, p0, p1)
    y = y.T
    #print(np.shape(y), y[0])
    
    #####
    ##### interpolate
    import scipy.interpolate
    Tsignal = []
    for i in range(nkid):
        masterid = kiddict[i]
        if masterid<0:
            Tsignal.append( [np.nan for i in range( len(fshift[i]) )] )
            continue
            
        #print(i, y[i], Tload[i])
        tck = scipy.interpolate.splrep(y[i], Tload[i], s=0)
        Tsignal.append( scipy.interpolate.splev(fshift[i], tck, der=0) )
    
    #####
    ##### convert to power
    Psignal = []
    for i in range(nkid):
        tck = scipy.interpolate.splrep(Tload[i], Pload[i], s=0)
        Psignal.append( scipy.interpolate.splev(Tsignal[i], tck, der=0) )
    
    return np.array(Tsignal), np.array(Psignal)

