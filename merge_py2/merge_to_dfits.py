# -*- coding utf-8 -*-
#
# merge_to_dfits.py: Read logging data and merge them into a FITS object
#
# Author : Tetsutaro Ueda
# Created: 2017/11/02
#-------------------------------- IMPORT MODULES
#-------- Standard Modules
import os
#-------- Dependent Packages
import yaml
import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy import coordinates
from astropy import units
import functions as fc

#-------------------------------- CONSTANTS
TELESCOP        = 'ASTE'
D_ASTE          = (10.0* units.m).value                        # Diameter  of the ASTE
LON_ASTE        = coordinates.Angle('-67d42m11.89525s').deg    # Longitude of the ASTE
LAT_ASTE        = coordinates.Angle('-22d58m17.69447s').deg    # Latitude  of the ASTE
FORM_FITSTIME   = '%Y-%m-%dT%H:%M:%S'                          # YYYY-mm-ddTHH:MM:SS
FORM_FITSTIME_P = '%Y-%m-%dT%H:%M:%S.%f'                       # YYYY-mm-ddTHH:MM:SS.ss
PATH_DFITSDICT  = os.path.expanduser('~/DESHIMA/devtools/merge_py2/dfits_dict.yaml')


class MergeToDfits:
    """Read logging data of ASTE and merge them into a FITS object

    Args.:
        ddbfits    (str): File name of the DESHIMA Database FITS.
        obsinst    (str): File name of the Observation Instruction.
        antennalog (str): File name of the Antennalogging.
        rout_data  (str): File name of the Readout Data (FITS).
        weatherlog (str): File name of the Weatherlogging.
                          Default is None.

    Example:
        "ddbfits"   : DDB_20171031.fits
        "obsinst"   : 20171024012916.obs
        "antennalog": 20171024012916.ant
        "rout_data" : reduced_20171024012916.fits
        "weatherlog": 20171024012916.wea

        Get HDUList of DFITS:
        >>> mtd = MergeToDfits(
                ddbfits=ddbfits,
                obsinst=obsinst,
                antennalog=antennalog,
                rout_data=rout_data
            )
        >>> dfits = mtd.dfits
        >>> dfits.indo()
        Filename: (No file associated with this HDUList)
        No.    Name         Type      Cards   Dimensions   Format
        0  PRIMARY     PrimaryHDU       4   ()
        1  OBSINFO     BinTableHDU     52   1R x 11C   ['K', 'D', 'D', 'D', 'D', 'D', 'D', '63K', '63K', '63K', '63D']
        2  ANTENNA     BinTableHDU     32   3854R x 8C   ['26A', '4A', 'D', 'D', 'D', 'D', 'D', 'D']
        3  READOUT     BinTableHDU     20   29234R x 4C   ['26A', 'K', '63D', '63D']

        Get each HDU:
        >>> mtd = MergeToDfits(
                ddbfits=ddbfits,
                obsinst=obsinst,
                antennalog=antennalog,
                rout_data=rout_data,
                weatherlog=weatherlog
            )
        >>> obsinfo = mtd.obsinfo
        >>> antenna = mtd.antenna
        >>> readout = mtd.readout
        >>> weather = mtd.weather
    """

    def __init__(self, ddbfits, obsinst, antennalog, rout_data, weatherlog=None):
#-------- Path
        self.ddbfits    = os.path.expanduser(ddbfits)
        self.obsinst    = os.path.expanduser(obsinst)
        self.antennalog = os.path.expanduser(antennalog)
        self.rout_data  = os.path.expanduser(rout_data)
        if weatherlog is None:
            self.wflag = 0
        else:
            self.weatherlog = os.path.expanduser(weatherlog)

#-------- DFITS Dictionary
        with open(PATH_DFITSDICT, 'r') as f:
            self.dfits_dict = yaml.load(f)

#-------- HDU
        self.antenna        # ANETNNA
        self.obsinfo        # OBSINFO
        self.readout        # READOUT
        if weatherlog is not None:
            self.weather    # WEATHER

#-------------------------------- METHODS
#---------------- Merge to DFITS
    @property
    def dfits(self):
        hdus = fits.HDUList()
        hdus.append(fits.PrimaryHDU())   # PRIMARY
        hdus.append(self.obsinfo)        # OBSINFO
        hdus.append(self.antenna)        # ANTENNA
        hdus.append(self.readout)        # READOUT
        if self.wflag != 0:
            hdus.append(self.weather)    # WEATHER

        return hdus

#---------------- OBSINFO
    @property
    def obsinfo(self):
#-------- Get the Dicitinary of 'OBSINFO': 'obsinfo_dict'
        od = self.dfits_dict['obsinfo_dict']
#-------- Get Header Values
        obsinst_data = fc.load_obsinst(self.obsinst)
        od['hdr_vals']['TELESCOP'] = TELESCOP
        od['hdr_vals']['SITELON']  = LON_ASTE
        od['hdr_vals']['SITELAT']  = LAT_ASTE
        od['hdr_vals']['DATE-OBS'] = self.ant_time[0][:19]
        od['hdr_vals']['OBSERVER'] = obsinst_data[0]
        od['hdr_vals']['OBJECT']   = obsinst_data[1]
        od['hdr_vals']['RA']       = np.mean(self.ant_ra)
        od['hdr_vals']['DEC']      = np.mean(self.ant_dec)
        od['hdr_vals']['EQUINOX']  = obsinst_data[2]
#-------- Get DDBID and Values for Columns
        with fits.open(self.ddbfits) as f:
            od['hdr_vals']['DDBID'] = f['PRIMARY'].header['DDB_ID']
            des  = f['KIDDES'].data
            filt = f['KIDFILT'].data
            rcp  = f['KIDRCP'].data

        od['col_vals']['pixelid']   = rcp['pixelid']
        od['col_vals']['offsetaz']  = rcp['offsetaz']
        od['col_vals']['offsetel']  = rcp['offsetel']
        od['col_vals']['gain']      = rcp['gain']
#/*--------------------------- Not confirmed ----------------------------*/
        od['col_vals']['interval']  = np.array([1/ 196])
        od['col_vals']['integtime'] = np.array([1/ 196])
        od['col_vals']['beamsize']  = np.array([0.005])     # 18 arcsec
#/*----------------------------------------------------------------------*/
        mas_kid_corresp = fc.get_maskid_corresp(des, filt)
        od['col_vals']['masterids'] = mas_kid_corresp[0]
        od['col_vals']['kidids']    = mas_kid_corresp[1]
        od['col_vals']['kidtypes']  = mas_kid_corresp[2]
        od['col_vals']['kidfreqs']  = mas_kid_corresp[3]

        return fc.create_bintablehdu(od)

#---------------- ANTENNA
    @property
    def antenna(self):
#-------- Get the Dicitinary of 'ANTENNA': 'antenna_dict'
        ad = self.dfits_dict['antenna_dict']
#---- Get Header Values
        ad['hdr_vals']['FILENAME'] = os.path.basename(self.antennalog)

#-------- Read 'antennalog'
        antlog_data = ascii.read(self.antennalog)[:-1]
        self.ant_time = fc.convert_asciitime(antlog_data['time'], FORM_FITSTIME_P)
        self.ant_ra   = antlog_data['ra-prg']
        self.ant_dec  = antlog_data['dec-prg']
#---- Get Values for Columns
        ad['col_vals']['time']      = self.ant_time
        ad['col_vals']['scantype']  = antlog_data['type']
        ad['col_vals']['az']        = antlog_data['az-real']
        ad['col_vals']['el']        = antlog_data['el-real']
        ad['col_vals']['ra']        = self.ant_ra
        ad['col_vals']['dec']       = self.ant_dec
        ad['col_vals']['az_center'] = antlog_data['az-prog(center)']
        ad['col_vals']['el_center'] = antlog_data['el-prog(center)']

        return fc.create_bintablehdu(ad)

#---------------- READOUT
    @property
    def readout(self):
#-------- Get the Dicitinary of 'READOUT': 'readout_dict'
        rd = self.dfits_dict['readout_dict']
#---- Get Header Values
        rd['hdr_vals']['FILENAME'] = os.path.basename(self.rout_data)
#-------- Open 'DDB' and 'rout_data'
        ddb   = fits.open(self.ddbfits)
        rhdus = fits.open(self.rout_data)
#---- Define Troom (Not confirmed)
        Troom = 17. + 273
#---- Get Values for Columns
        rd['col_vals']['starttime'] = fc.convert_timestamp(rhdus['READOUT'].data['timestamp'])
        rd['col_vals']['pixelid']   = rhdus['READOUT'].data['pixelid']
        rd['col_vals']['Tsignal'], rd['col_vals']['Psignal'] = fc.calibrate_to_power(0, Troom, rhdus, ddb)
#-------- Close 'DDB' and 'rout_data'
        ddb.close()
        rhdus.close()

        return fc.create_bintablehdu(rd)

#---------------- WEATHER
    @property
    def weather(self):
#-------- Error Handling: Case of 'weatherlog' is None
        if self.wflag == 0:
            raise ValueError('No "weatherlog" is inputed!!')
#-------- Get the Dicitinary of 'READOUT': 'readout_dict'
        wd = self.dfits_dict['weather_dict']
#---- Get Header Values
        wd['hdr_vals']['FILENAME'] = os.path.basename(self.weatherlog)
#-------- Read 'weatherlog'
        wlog_data = ascii.read(self.weatherlog)
#---- Get Values for Columns
        wd['col_vals']['time']           = fc.convert_asciitime(wlog_data['time'], FORM_FITSTIME)
        wd['col_vals']['temperature']    = wlog_data['tmperature']
        wd['col_vals']['pressure']       = wlog_data['presure']
        wd['col_vals']['vapor-pressure'] = wlog_data['vapor-pressure']
        wd['col_vals']['windspd']        = wlog_data['aux1']
        wd['col_vals']['winddir']        = wlog_data['aux2']

        return fc.create_bintablehdu(wd)
