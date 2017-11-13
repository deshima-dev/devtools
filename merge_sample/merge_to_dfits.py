# -*- coding utf-8 -*-
#
# merge_to_dfits.py: Read logging data and merge them into a FITS object
#
# Author : Tetsutaro Ueda
# Created: 2017/11/02
#-------------------------------- IMPORT MODULES
#-------- Standard Modules
from pathlib import Path
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
PATH_DFITSDICT  = Path('~/DESHIMA/devtools/merge_sample/dfits_dict.yaml').expanduser()


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
        self.ddbfits    = Path(ddbfits).expanduser()
        self.obsinst    = Path(obsinst).expanduser()
        self.antennalog = Path(antennalog).expanduser()
        self.rout_data  = Path(rout_data).expanduser()
        if weatherlog is None:
            self.wflag = 0
        else:
            self.weatherlog = Path(weatherlog).expanduser()

#-------- DFITS Dictionary
        with PATH_DFITSDICT.open() as f:
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
        """
            HDU list of DFITS
        """
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
        """
            HDU of 'OBSINFO'
        """
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
        """
            HDU of 'ANTENNA'
        """
#-------- Get the Dicitinary of 'ANTENNA': 'antenna_dict'
        ad = self.dfits_dict['antenna_dict']
#---- Get Header Values
        ad['hdr_vals']['FILENAME'] = self.antennalog.name

#-------- Read 'antennalog'
        antlog_data = ascii.read(self.antennalog)[:-1]
        self.ant_time = fc.convert_asciitime(antlog_data['time'])
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
        """
            HDU of 'READOUT'
        """
#-------- Get the Dicitinary of 'READOUT': 'readout_dict'
        rd = self.dfits_dict['readout_dict']
#---- Get Header Values
        rd['hdr_vals']['FILENAME'] = self.rout_data.name
#-------- Open 'DDB' and 'rout_data'
        ddb   = fits.open(self.ddbfits)
        rhdus = fits.open(self.rout_data)
#---- Define Troom (Not confirmed)
        Troom = 17. + 273
#---- Get Values for Columns
        nkid = rhdus['READOUT'].header['NKID0']
        reduce_data = np.transpose(
            [rhdus['READOUT'].data['Amp, Ph, linPh %d' %i].T for i in range(nkid)]
        )

        rd['col_vals']['starttime']  = fc.convert_timestamp(rhdus['READOUT'].data['timestamp'])
        rd['col_vals']['pixelid']    = rhdus['READOUT'].data['pixelid']
        rd['col_vals']['amplitude']  = reduce_data[:, 0]
        rd['col_vals']['phase']      = reduce_data[:, 1]
        rd['col_vals']['line_phase'] = reduce_data[:, 2]
        rd['col_vals']['Tsignal'], rd['col_vals']['Psignal'] = fc.calibrate_to_power(0, Troom, rhdus, ddb)
#-------- Close 'DDB' and 'rout_data'
        ddb.close()
        rhdus.close()

        return fc.create_bintablehdu(rd)

#---------------- WEATHER
    @property
    def weather(self):
        """
            HDU of 'WEATHER'
        """
#-------- Error Handling: Case of 'weatherlog' is None
        if self.wflag == 0:
            raise ValueError('No "weatherlog" is inputed!!')
#-------- Get the Dicitinary of 'READOUT': 'readout_dict'
        wd = self.dfits_dict['weather_dict']
#---- Get Header Values
        wd['hdr_vals']['FILENAME'] = self.weatherlog.name
#-------- Read 'weatherlog'
        wlog_data = ascii.read(self.weatherlog)
#---- Get Values for Columns
        wd['col_vals']['time']           = fc.convert_asciitime(wlog_data['time'])
        wd['col_vals']['temperature']    = wlog_data['tmperature']
        wd['col_vals']['pressure']       = wlog_data['presure']
        wd['col_vals']['vapor-pressure'] = wlog_data['vapor-pressure']
        wd['col_vals']['windspd']        = wlog_data['aux1']
        wd['col_vals']['winddir']        = wlog_data['aux2']

        return fc.create_bintablehdu(wd)
