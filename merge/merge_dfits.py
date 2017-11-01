#  -*- coding: utf-8 -*-
#
#  merge_dfits.py: Merge All HDUs and create 'DFITS'
#
#  Author : Tetsutaro Ueda
#  Created: 2017/10/24
#-------------------------------- IMPORT MODULES
#-------- Standard Modules
import os
from pathlib import Path
#-------- Dependent Packages
import yaml
import numpy as np
from datetime import datetime
from astropy.io import fits
from astropy.io import ascii
from astropy import table
from astropy import coordinates
from astropy import units


#-------------------------------- CONSTANTS
TELESCOP        = 'ASTE'
D_ASTE          = (10.0* units.m).value                        # Diameter  of the ASTE
LON_ASTE        = coordinates.Angle('-67d42m11.89525s').deg    # Longitude of the ASTE
LAT_ASTE        = coordinates.Angle('-22d58m17.69447s').deg    # Latitude  of the ASTE
FORM_FITSTIME   = '%Y-%m-%dT%H:%M:%S'                          # YYYY-mm-ddTHH:MM:SS
FORM_FITSTIME_P = '%Y-%m-%dT%H:%M:%S.%f'                       # YYYY-mm-ddTHH:MM:SS.ss
PATH_DFITSDICT  = str(Path('~/DESHIMA/devtools/merge/DFITS_dict.yaml').expanduser())

#-------------------------------- PUBLIC ITEMS
#__all__ = ['fromaste']


#-------------------------------- FUNCTIONS
#------------------------ Main
def dfits_fromaste(ddb_fits, obsinst, antennalog, rout_data, weatherlog=None):
    """Read logging data of ASTE and merge them into a FITS object.

    Args.:
        ddb_fits     (str): File name of the DESHIMA Database FITS.
        obsinst      (str): File name of the Observation Instruction.
        antennalog   (str): File name of the Antennalogging.
        readout_data (str): File name of the Readout Data.
        weatherlog   (str): File name of the Weatherlogging.
                            Default is None.

    Returns:
        hdus (HDUlist): HDU list containing the merged data.
    """

#---------------- Path
    ddb_fits   = str(Path(ddb_fits).expanduser())
    obsinst    = str(Path(obsinst).expanduser())
    antennalog = str(Path(antennalog).expanduser())
    rout_data  = str(Path(rout_data).expanduser())
    if weatherlog is not None:
        weatherlog = str(Path(weatherlog).expanduser())

#---------------- Create DFITS HDUs
    hdus = fits.HDUList()

    hdus.append(fits.PrimaryHDU())                              # PRIMARY
    hdus.append(make_obsinfo(ddb_fits, obsinst, antennalog))    # OBSINFO
    hdus.append(make_antenna(antennalog))                       # ANTENNA
    hdus.append(make_readout(ddb_fits, rout_data))                        # READOUT

    if weatherlog is not None:
        hdus.append(make_weather(weatherlog))                   # WEATHER

    return hdus


#------------------------ Create HDUs
#---------------- OBSINFO
def make_obsinfo(ddb_fits, obsinst, antennalog):
    #-------- Read 'obsinst'
    with open(obsinst, 'r') as f:
        for line in f:
            if '% OBSERVER' in line:
                observer   = line.rstrip().split('=')[1]                          # Get OBSERVER
            elif '% SRC_NAME' in line:
                obs_object = line.rstrip().split('=')[1]                          # Get OBJECT
            elif '% EPOCH' in line:
                equinox    = line.rstrip().split('=')[1].strip('J').strip('B')    # Get EQUINOX

    #-------- Read 'antennalog'
    antennalog = ascii.read(antennalog)
    #---- Get 'Date-obs' from 'antennalog'
    date_obs   = datetime.strptime(antennalog['time'][0].astype(np.str), '%Y%m%d%H%M%S.%f')
    date_obs   = datetime.strftime(date_obs, FORM_FITSTIME)
    #---- Get 'RA' and 'DEC'
    ra  = np.mean(antennalog['ra-prg'][:-1])
    dec = np.mean(antennalog['dec-prg'][:-1])

    dh = fits.open(ddb_fits)
    pixelid   = dh['KIDRCP'].data['pixelid']
    offsetaz  = dh['KIDRCP'].data['offsetaz']
    offsetel  = dh['KIDRCP'].data['offsetel']
    interval  = np.array([1/ 196])
    integtime = np.array([1/ 196])
    beamsize  = np.array([0.005])
    gain      = dh['KIDRCP'].data['gain']
    mas_vs_attr = [[x, y]
                   for (x, y) in zip(
                       list(dh['KIDDES'].data['masterid']),
                       list(dh['KIDDES'].data['attribute']))]
    del mas_vs_attr[13:17]
    mas_vs_attr.append([-1, np.nan])
    mas_vs_kid  = [[x, y, z[0]]
                   for (x, y, z) in zip(
                       dh['KIDFILT'].data['masterid'],
                       dh['KIDFILT'].data['kidid'],
                       dh['KIDFILT'].data['F_filter, dF_filter'])]
    mas_kid_corresp = [[x[0], x[1], y[1], x[2]]
                       for (x, y) in zip(
                           sorted(mas_vs_kid, key=lambda x: x[0]),
                           sorted(mas_vs_attr, key=lambda x: x[0]))]
    mas_kid_corresp = sorted(mas_kid_corresp, key=lambda x: x[1])

    ddbid = dh['PRIMARY'].header['DDB_ID']
    dh.close()
    mas_kid_corresp = list(map(list, zip(*mas_kid_corresp)))
    masterids = [mas_kid_corresp[0]]
    kidids    = [mas_kid_corresp[1]]
    kidattr   = mas_kid_corresp[2]
    kidfreqs  = [mas_kid_corresp[3]]

    kidtypes = []
    for x in kidattr:
        if x=='wideband':   attr = 0
        elif x=='filter':   attr = 1
        elif x=='blind':    attr = 2
        else:               attr = np.nan
        kidtypes.append(attr)
    kidtypes = [kidtypes]
    obsinfo_data_lis = [pixelid, offsetaz, offsetel, interval, integtime, beamsize, gain,
                        masterids, kidids, kidtypes, kidfreqs]
     #-------- Read the Dictionary 'obsinfo_dict'
    with open(PATH_DFITSDICT, 'r') as f:
        obsinfo_dict = yaml.load(f)['obsinfo_dict']

    #-------- Set Data to the Dictinary 'obsinfo_dict'
    obsinfo_dict['hdr_val_lis'][2:12] = [ddbid, TELESCOP, LON_ASTE, LAT_ASTE, date_obs,
                                         observer, obs_object, ra, dec, equinox]
    obsinfo_dict['cols_data_lis']     = obsinfo_data_lis

    return createBinTableHDU(obsinfo_dict)

#---------------- ANTENNA
def make_antenna(antennalog):
    filename = os.path.basename(antennalog)

#-------- Read 'antennalog'
    antlog_data = ascii.read(antennalog)
    antlog_len  = len(antlog_data) - 1
    antlog_data = antlog_data[:antlog_len]
    ant_time    = convert_asciitime(antlog_data['time'].astype(np.str))

#-------- Read the Dictionary 'antenna_dict'
    with open(PATH_DFITSDICT, 'r') as f:
        antenna_dict = yaml.load(f)['antenna_dict']

#-------- Set Data to the Dictinary 'obsinfo_dict'
    antenna_dict['hdr_val_lis'][1] = filename
    antenna_dict['cols_data_lis']  = [ant_time, antlog_data['type'],
                                      antlog_data['az-real'],
                                      antlog_data['el-real'],
                                      antlog_data['ra-prg'],
                                      antlog_data['dec-prg'],
                                      antlog_data['az-prog(center)'],
                                      antlog_data['el-prog(center)']]

#-------- Create 3rd. HDU 'ANTENNA'
    return createBinTableHDU(antenna_dict)


#---------------- READOUT (Not confirmed)
def make_readout(ddb_fits, rout_data):
    filename  = os.path.basename(rout_data)

#-------- Read 'ddb_fits' and 'rout_data'
    ddb   = fits.open(ddb_fits)
    rhdus = fits.open(rout_data)

    starttime = convert_timestamp(rhdus['READOUT'].data['timestamp'])
    pixelid_0 = rhdus['READOUT'].data['pixelid']

#---- Caliblate 'Amplitude' and 'Phase' to 'Power'
    Troom = 17. + 273.    # K, cabin temperature
    Tsignal, Psignal = calibrate_to_power(pixelid_0[0], Troom, rhdus, ddb)
    Tsignal = Tsignal.T
    Psignal = Psignal.T

#---- Close 'ddb_fits' and 'rout_data'
    ddb.close()
    rhdus.close()

#-------- Read the Dictionary 'readout_dict'
    with open(PATH_DFITSDICT, 'r') as f:
        readout_dict = yaml.load(f)['readout_dict']

#-------- Set Data to the Dictinary 'readout_dict'
    readout_dict['hdr_val_lis'][1] = os.path.basename(rout_data)
    readout_dict['cols_data_lis']  = [starttime, pixelid_0, Tsignal, Psignal]

#-------- Create 4th. HDU 'READOUT'
    return createBinTableHDU(readout_dict)


#/*--------------------------- Not confirmed ----------------------------*/
#---------------- FILTERS
def make_filters(ddb_fits):
    filename = os.path.basename(ddb_fits)

#-------- Read 'filters_data'
    with fits.open(ddb_fits) as ddb_hdus:
        filename      = ddb_hdus['KIDSINFO'].header['FILENAME']
        kidsinfo_data = ddb_hdus['KIDSINFO'].data

#-------- Read the Dictionary 'filters_dict'
    with open(PATH_DFITSDICT, 'r') as f:
        filters_dict = yaml.load(f)['filters_dict']

#-------- Set Data to the Dictinary 'readout_dict'
    crval3 = 320.0e+09
    crpix3 = 0.0
    cdelt3 = 10.0e+06

    filters_dict['hdr_val_lis'][1]  = filename
    filters_dict['hdr_val_lis'][3:] = [crval3, crpix3, cdelt3]
    filters_dict['cols_data_lis']   = [kidsinfo_data['pixelid'],
                                       kidsinfo_data['kidid'],
                                       kidsinfo_data['kidtypes'],
                                       kidsinfo_data['bandpass'],
                                       kidsinfo_data['lorentz']]
    filters_dict['tform'][3] = str(len(kidsinfo_data['bandpass'][0])) + 'D'

#-------- Create 5th. HDU 'FILTERS'
    return createBinTableHDU(filters_dict)
#/*----------------------------------------------------------------------*/


#---------------- WEATHER
def make_weather(weatherlog):
    filename = os.path.basename(weatherlog)

#-------- Read 'weatherlog'
    wealog_data = ascii.read(weatherlog)
    wea_time = convert_asciitime(wealog_data['time'].astype(np.str))

#-------- Set Data to the Dictinary 'weather_dict'
    with open(PATH_DFITSDICT, 'r') as f:
        weather_dict = yaml.load(f)['weather_dict']

    weather_dict['hdr_val_lis'][1] = filename
    weather_dict['cols_data_lis']  = [wea_time,
                                      wealog_data['tmperature'],
                                      wealog_data['presure'],
                                      wealog_data['vapor-pressure'],
                                      wealog_data['aux1'],
                                      wealog_data['aux2']]

#---- Create 6th. HDU 'WEATHER'
    return createBinTableHDU(weather_dict)


#------------------------ Others
#---------------- Create Binary Table HDU
def createBinTableHDU(data_dict):
#-------- Set Header and Comments
    header = fits.Header()
    for (i, j, k) in zip(data_dict['hdr_key_lis'], data_dict['hdr_val_lis'], data_dict['hdr_com_lis']):
        header[i] = j, k

#-------- Create Collumns of the Binary Table
    columns = [fits.Column(name=data_dict['cols_key_lis'][i],
                           format=data_dict['tform'][i],
                           array=data_dict['cols_data_lis'][i],
                           unit=data_dict['tunit'][i])
               for i in range(len(data_dict['cols_key_lis']))]

    hdu = fits.BinTableHDU.from_columns(columns, header)

#-------- Add comments
    [addHeaderComments(hdu.header, list(hdu.header)[i], data_dict['hdr_com_lis'][i])
     for i in range(-(len(data_dict['hdr_com_lis']) - data_dict['hdr_com_lis'].index('label for field 1')), 0)]

    return hdu


#---------------- Add Comments to Header
def addHeaderComments(hdr, key, com):
    hdr.comments[key] = com


#---------------- Chunk the List
def chunk_list(iterable, n):
    return [iterable[i: i+n] for i in range(0, len(iterable), n)]


#---------------- Convert timestamps Applied FITS Times
def convert_timestamp(timestamp):
    timestamp = [datetime.utcfromtimestamp(t) for t in timestamp]
    timestamp = [datetime.strftime(t, FORM_FITSTIME_P) for t in timestamp]
    return np.array(timestamp)


#-------- Convert ASCII Times Applied FITS Times
def convert_asciitime(asciitime):
    if len(asciitime[0]) <= 14:
        asciitime = [datetime.strptime(t, '%Y%m%d%H%M%S') for t in asciitime]
        form_fitstime = FORM_FITSTIME
    else:
        asciitime = [datetime.strptime(t, '%Y%m%d%H%M%S.%f') for t in asciitime]
        form_fitstime = FORM_FITSTIME_P

    asciitime = [datetime.strftime(t, form_fitstime) for t in asciitime]
    return np.array(asciitime)

def func(x, p0, p1):
    return p0*np.sqrt(x) + p1
    #return p0*(np.sqrt(x)-np.sqrt(x[0])) + p1

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
