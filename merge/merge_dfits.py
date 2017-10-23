#  -*- coding: utf-8 -*-
#
#  DFITS_merge.py: Merge All HDUs and create 'DFITS'
#
#  Author : Tetsutaro Ueda
#  Created: 2017/10/16
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
__all__ = ['fromaste']


#-------------------------------- FUNCTIONS
#---------------- Main
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

#-------- Path
    ddb_fits   = str(Path(ddb_fits).expanduser())
    obsinst    = str(Path(obsinst).expanduser())
    antennalog = str(Path(antennalog).expanduser())
    rout_data  = str(Path(rout_data).expanduser())
    if weatherlog is not None:
        weatherlog = str(Path(weatherlog).expanduser())

#-------- Create DFITS HDUs
    hdus = fits.HDUList()

    hdus.append(fits.PrimaryHDU())                              # PRIMARY
    hdus.append(make_obsinfo(ddb_fits, obsinst, antennalog))    # OBSINFO
    hdus.append(make_antenna(antennalog))                       # ANTENNA
    hdus.append(make_readout(rout_data))                        # READOUT
#    hdus.append(make_filters(ddb_fits))                         # FILTERS
    if weatherlog is not None:
        hdus.append(make_weather(weatherlog))                   # WEATHER

    return hdus


#---------------- Create HDUs
#-------- OBSINFO
def make_obsinfo(ddb_fits, obsinst, antennalog):
#---- Read 'ddb_fits'
    with fits.open(ddb_fits) as f:
        ddb_id        = f['PRIMARY'].header['DDB_ID']
        kidsinfo_data = f['KIDSINFO'].data
        rcp_data      = f['RCP'].data

#---- Read 'obsinst'
    with open(obsinst, 'r') as f:
        for line in f:
            if '% OBSERVER' in line:
                observer   = line.rstrip().split('=')[1]                          # Get OBSERVER
            elif '% SRC_NAME' in line:
                obs_object = line.rstrip().split('=')[1]                          # Get OBJECT
            elif '% EPOCH' in line:
                equinox    = line.rstrip().split('=')[1].strip('J').strip('B')    # Get EQUINOX

#---- Read 'antennalog'
    antennalog = ascii.read(antennalog)
    date_obs = datetime.strptime(antennalog['time'][0].astype(np.str), '%Y%m%d%H%M%S.%f')
    date_obs = datetime.strftime(date_obs, FORM_FITSTIME)

#/*--------------------------- Not confirmed ----------------------------*/
    ra  = 0
    dec = 0
    interval  = np.ones(1)* 1/ 196
    integtime = np.ones(1)* 1/ 196
#    beamsize  = np.ones(1)* 1* 10** -3/ D_ASTE  #D/λ
    beamsize  = np.ones(1)* 0.005    # 18 arcsec
#    kidfreqs  = np.arange(320.0e+09, 369.0e+09, 1.0e+09)
    kidfreqs  = np.arange(320.0e+09, 376.5e+09, 0.5e+09)
#/*----------------------------------------------------------------------*/
    kid_num   = len(kidsinfo_data['kidid'])
    kidfreqs  = kidfreqs + np.transpose(kidsinfo_data['lorentz'])[0]

    obsinfo_data_lis = [rcp_data['pixelid'],
                        rcp_data['offsetaz'],
                        rcp_data['offsetel'],
                        interval,
                        integtime,
                        beamsize,
                        rcp_data['gain'],
                        chunk_list(kidsinfo_data['kidid'], kid_num),
                        chunk_list(kidsinfo_data['kidtypes'], kid_num),
                        chunk_list(kidfreqs, kid_num)]

#---- Read the Dictionary 'obsinfo_dict'
    with open(PATH_DFITSDICT, 'r') as f:
        obsinfo_dict = yaml.load(f)['obsinfo_dict']

#---- Set Data to the Dictinary 'obsinfo_dict'
    obsinfo_dict['hdr_val_lis'][2:12] = [ddb_id, TELESCOP, LON_ASTE, LAT_ASTE, date_obs,
                                         observer, obs_object, ra, dec, equinox]
    obsinfo_dict['cols_data_lis']     = obsinfo_data_lis
    obsinfo_dict['tform'][7:10]       = [str(kid_num) + 'K',
                                         str(kid_num) + 'K',
                                         str(kid_num) + 'D']

#---- Create 2nd HDU 'OBSINFO'
    return createBinTableHDU(obsinfo_dict)


#-------- ANTENNA
def make_antenna(antennalog):
#---- Read 'antennalog'
    filename = os.path.basename(antennalog)

    antlog_data = ascii.read(antennalog)
    antlog_len  = len(antlog_data) - 1
    antlog_data = antlog_data[:antlog_len]
    ant_time    = convert_asciitime(antlog_data['time'].astype(np.str))

#---- Read the Dictionary 'antenna_dict'
    with open(PATH_DFITSDICT, 'r') as f:
        antenna_dict = yaml.load(f)['antenna_dict']

#---- Set Data to the Dictinary 'obsinfo_dict'
    antenna_dict['hdr_val_lis'][1] = filename
    antenna_dict['cols_data_lis']  = [ant_time, antlog_data['type'],
                                      antlog_data['az-real'],
                                      antlog_data['el-real'],
                                      antlog_data['ra-prg'],
                                      antlog_data['dec-prg'],
                                      antlog_data['az-prog(center)'],
                                      antlog_data['el-prog(center)']]

#---- Create 3rd. HDU 'ANTENNA'
    return createBinTableHDU(antenna_dict)


#-------------------- Not confirmed (for First Light) --------------------
#-------- READOUT (Not confirmed)
def make_readout(rout_data):
#---- Read 'readout_data'
    filename  = os.path.basename(rout_data)

    with fits.open(rout_data) as f:
        timestamp = f[1].data['timestamp']
        I = f[1].data['data'][:, 0::2]
        Q = f[1].data['data'][:, 1::2]
        phase = np.arctan(Q/I)

    starttime = convert_timestamp(timestamp)
    pixelid   = np.zeros(len(starttime))
    arraydata = phase

#---- Read the Dictionary 'readout_dict'
    with open(PATH_DFITSDICT, 'r') as f:
        readout_dict = yaml.load(f)['readout_dict']

#---- Set Data to the Dictinary 'readout_dict'
    readout_dict['hdr_val_lis'][1] = os.path.basename(rout_data)
    readout_dict['cols_data_lis']  = [starttime, pixelid, arraydata]
    readout_dict['tform'][2]       = str(len(arraydata[0])) + 'D'

#---- Create 4th. HDU 'READOUT'
    return createBinTableHDU(readout_dict)


#/*--------------------------- Not confirmed ----------------------------*/
#-------- FILTERS
def make_filters(ddb_fits):
#---- Read 'filters_data'
    with fits.open(ddb_fits) as ddb_hdus:
        filename      = ddb_hdus['KIDSINFO'].header['FILENAME']
        kidsinfo_data = ddb_hdus['KIDSINFO'].data

#---- Read the Dictionary 'filters_dict'
    with open(PATH_DFITSDICT, 'r') as f:
        filters_dict = yaml.load(f)['filters_dict']

#---- Set Data to the Dictinary 'readout_dict'
    crval3 = 320.0e+09
    crpix3 = 0.0
    cdelt3 = 10.0e+06

#    filters_dict['hdr_val_lis'][1] = os.path.basename(filters_data)
    filters_dict['hdr_val_lis'][1]  = filename
    filters_dict['hdr_val_lis'][3:] = [crval3, crpix3, cdelt3]
    filters_dict['cols_data_lis']   = [kidsinfo_data['pixelid'],
                                       kidsinfo_data['kidid'],
                                       kidsinfo_data['kidtypes'],
                                       kidsinfo_data['bandpass'],
                                       kidsinfo_data['lorentz']]
    filters_dict['tform'][3] = str(len(kidsinfo_data['bandpass'][0])) + 'D'

#---- Create 5th. HDU 'FILTERS'
    return createBinTableHDU(filters_dict)
#/*----------------------------------------------------------------------*/


#-------- WEATHER
def make_weather(weatherlog):
#---- Read 'weatherlog'
    filename = os.path.basename(weatherlog)

    wealog_data = ascii.read(weatherlog)
    wea_time = convert_asciitime(wealog_data['time'].astype(np.str))

#---- Set Data to the Dictinary 'weather_dict'
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


#---------------- Others
#-------- Create Binary Table HDU
def createBinTableHDU(data_dict):
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
     for i in range(-(len(data_dict['hdr_com_lis']) - data_dict['hdr_com_lis'].index('label for field 1')), 0)]

    return hdu


#-------- Add Comments to Header
def addHeaderComments(hdr, key, com):
    hdr.comments[key] = com


#-------- Chunk the List
def chunk_list(iterable, n):
    return [iterable[i: i+n] for i in range(0, len(iterable), n)]


#-------------------------------- Not confirmed ------------------------
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
