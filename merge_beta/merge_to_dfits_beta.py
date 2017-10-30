#  -*- coding: utf-8 -*-
#
#  merge_to_dfits_beta.py: Merge All HDUs and create 'DFITS' (Beta)
#
#  Author : Tetsutaro Ueda
#  Created: 2017/10/31
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
PATH_DFITSDICT  = str(Path('~/DESHIMA/devtools/merge_beta/DFITS_dict.yaml').expanduser())

#-------------------------------- PUBLIC ITEMS
__all__ = ['dfits_fromaste']


#-------------------------------- FUNCTIONS
#------------------------ Main
def dfits_fromaste(kid_mas_corresp, obsinst, antennalog, rout_data, weatherlog=None):
    """Read logging data of ASTE and merge them into a FITS object.

    Args.:
        kid_mas_corresp (str): File name of Correpondance table of 'kid_id' and 'master_id'.
        obsinst         (str): File name of the Observation Instruction.
        antennalog      (str): File name of the Antennalogging.
        readout_data    (str): File name of the Readout Data.
        weatherlog      (str): File name of the Weatherlogging.
                               Default is None.

    Returns:
        hdus (HDUlist): HDU list containing the merged data.
    """

#---------------- Path
    kid_mas_corresp = str(Path(kid_mas_corresp).expanduser())
    obsinst         = str(Path(obsinst).expanduser())
    antennalog      = str(Path(antennalog).expanduser())
    rout_data       = str(Path(rout_data).expanduser())
    if weatherlog is not None:
        weatherlog  = str(Path(weatherlog).expanduser())

#---------------- Create DFITS HDUs
    hdus = fits.HDUList()

    hdus.append(fits.PrimaryHDU())                                     # PRIMARY
    hdus.append(make_obsinfo(kid_mas_corresp, obsinst, antennalog))    # OBSINFO
    hdus.append(make_antenna(antennalog))                              # ANTENNA
    hdus.append(make_readout(rout_data))                               # READOUT
    if weatherlog is not None:
        hdus.append(make_weather(weatherlog))                          # WEATHER

    return hdus


#------------------------ Create HDUs
#---------------- OBSINFO
def make_obsinfo(kid_mas_corresp, obsinst, antennalog):
#-------- Read 'ddb_fits'
    '''
    with fits.open(ddb_fits) as f:
        ddb_id        = f['PRIMARY'].header['DDB_ID']
        kidsinfo_data = f['KIDSINFO'].data
        rcp_data      = f['RCP'].data
    '''
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

#-------- Read 'kid_master_corresp.csv'
    kid_vs_mas = ascii.read(kid_mas_corresp)
#---- Get kid inromation
    kidids   = np.array([list(kid_vs_mas['kid_id'])])
    kid_attr = list(kid_vs_mas['attribute'])
    kidtypes = []
    for x in kid_attr:
        if x=='wideband':   attr = 0
        elif x=='filter':   attr = 1
        elif x=='blind':    attr = 2
        else:               attr = np.nan
        kidtypes.append(attr)
    kidtypes = np.array([kidtypes])
    kidfreqs = np.array([list(kid_vs_mas['F_filter_meas'])])

#/*--------------------------- Not confirmed ----------------------------*/
    pixelid   = np.array([0])
    offsetaz  = np.array([0.0])
    offsetel  = np.array([0.0])
    interval  = np.ones(1)*(1/ 196)
    integtime = np.ones(1)*(1/ 196)
    beamsize  = np.ones(1)*0.005     # 18 arcsec
    gain      = np.array([1.0])
#/*----------------------------------------------------------------------*/

#-------- Read the Dictionary 'obsinfo_dict'
    with open(PATH_DFITSDICT, 'r') as f:
        obsinfo_dict = yaml.load(f)['obsinfo_dict']

#-------- Set Data to the Dictinary 'obsinfo_dict'
    obsinfo_dict['hdr_val_lis'][2:12] = [None, TELESCOP, LON_ASTE, LAT_ASTE, date_obs,
                                         observer, obs_object, ra, dec, equinox]
    obsinfo_dict['cols_data_lis']     = [pixelid, offsetaz, offsetel, interval,
                                         integtime, beamsize, gain,
                                         kidids, kidtypes, kidfreqs]

#-------- Create 2nd HDU 'OBSINFO'
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
def make_readout(rout_data):
    filename  = os.path.basename(rout_data)

#-------- Read 'readout_data'
    with fits.open(rout_data) as f:
        starttime = convert_timestamp(f[1].data['timestamp'])
        pixelid_0 = np.zeros(len(starttime))
        amplitude = np.transpose([f[1].data['Amp %d' %i] for i in range(0, 63)])
        phase     = np.transpose([f[1].data['Ph %d' %i] for i in range(0, 63)])

#-------- Read the Dictionary 'readout_dict'
    with open(PATH_DFITSDICT, 'r') as f:
        readout_dict = yaml.load(f)['readout_dict']

#-------- Set Data to the Dictinary 'readout_dict'
    readout_dict['hdr_val_lis'][1] = os.path.basename(rout_data)
    readout_dict['cols_data_lis']  = [starttime, pixelid_0, amplitude, phase]

#-------- Create 4th. HDU 'READOUT'
    return createBinTableHDU(readout_dict)


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
