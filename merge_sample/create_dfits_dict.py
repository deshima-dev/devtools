# -*- coding utf-8 -*-
#
# create_dfits_dict.py: Create the Dictionary for DFITS
#
# Author : Tetsutaro Ueda
# Created: 2017/11/02
#-------------------------------- IMPORT MODULES
#-------- Standard Modules
from pathlib import Path
from collections import OrderedDict
#-------- Dependent Packages
import yaml

#-------------------------------- DEFINE CONSTANTS
PATH_DFITSDICT = Path('~/DESHIMA/devtools/merge_sample/dfits_dict.yaml').expanduser()

#-------------------------------- CREATE DICTIONARIES
#---------------- OBSINFO
def obsinfo_dict():
    hdr_vals = OrderedDict([
        ('EXTNAME', 'OBSINFO'), ('FITSTYPE', 'DESHIMAv0'), ('DDBID', None),
        ('TELESCOP', None), ('SITELON', None), ('SITELAT', None),
        ('DATE-OBS', None), ('OBSERVER', None), ('OBJECT', None),
        ('RA', None), ('DEC', None), ('EQUINOX', None), ('KIDTYPE0', 'wideband'),
        ('KIDTYPE1', 'filter'), ('KIDTYPE2', 'blind')
    ])
    hdr_coms = {
        'EXTNAME': 'name of binary table',
        'FITSTYPE': 'declares DESHIMA FITS',
        'DDBID': 'ID of DESHIMA database',
        'TELESCOP': 'name of used telescope',
        'SITELON': 'site longitude in units of deg',
        'SITELAT': 'site latitude in units of deg',
        'DATE-OBS': 'YYYY-mm-ddTHH:MM:SS',
        'OBSERVER': 'name of observer',
        'OBJECT': 'name of observed object',
        'RA': 'right ascension of the object in units of dec',
        'DEC': 'declination of the object in units of deg',
        'EQUINOX': 'equinox of coordinates',
        'KIDTYPE0': '0 is the kid of "wideband"',
        'KIDTYPE1': '1 is the kid of "filter"',
        'KIDTYPE2': '2 is the kid of "blind"',
        'TTYPE1': 'label for field 1',
        'TFORM1': 'data format of field 1',
        'TTYPE2': 'label for field 2',
        'TFORM2': 'data format of field 2',
        'TUNIT2': 'data unit of field 2',
        'TTYPE3': 'label for field 3',
        'TFORM3': 'data format of field 3',
        'TUNIT3': 'data unit of field 3',
        'TTYPE4': 'label for field 4',
        'TFORM4': 'data format of field 4',
        'TUNIT4': 'data unit of field 4',
        'TTYPE5': 'label for field 5',
        'TFORM5': 'data format of field 5',
        'TUNIT5': 'data unit of field 5',
        'TTYPE6': 'label for field 6',
        'TFORM6': 'data format of field 6',
        'TUNIT6': 'data unit of field 6',
        'TTYPE7': 'label for field 7',
        'TFORM7': 'data format of field 7',
        'TUNIT7': 'data unit of field 7',
        'TTYPE8': 'label for field 8',
        'TFORM8': 'data format of field 8',
        'TTYPE9': 'label for field 9',
        'TFORM9': 'data format of field 9',
        'TTYPE10': 'label for field 10',
        'TFORM10': 'data format of field 10',
        'TTYPE11': 'label for field 11',
        'TFORM11': 'data format of field 11',
        'TUNIT11': 'data unit of field 11'
    }
    col_vals = OrderedDict([
        ('pixelid', None), ('offsetaz', None), ('offsetel', None),
        ('interval', None), ('integtime', None), ('beamsize', None),
        ('gain', None), ('masterids', None), ('kidids', None),
        ('kidtypes', None), ('kidfreqs', None)
    ])
    col_form = OrderedDict([
        ('pixelid', 'K'), ('offsetaz', 'D'), ('offsetel', 'D'),
        ('interval', 'D'), ('integtime', 'D'), ('beamsize', 'D'),
        ('gain', 'D'), ('masterids', '63K'), ('kidids', '63K'),
        ('kidtypes', '63K'), ('kidfreqs', '63D')
    ])
    col_unit = OrderedDict([
        ('pixelid', None), ('offsetaz', 'deg'), ('offsetel', 'deg'),
        ('interval', 's'), ('integtime', 's'), ('beamsize', 'deg'),
        ('gain', '1'), ('masterids', None), ('kidids', None),
        ('kidtypes', None), ('kidfreqs', 'GHz')
    ])

    obsinfo_dict = {
        'hdr_vals': hdr_vals, 'hdr_coms': hdr_coms, 'col_vals': col_vals,
        'col_form': col_form, 'col_unit': col_unit
    }
    return obsinfo_dict

#---------------- ANTENNA
def antenna_dict():
    hdr_vals = OrderedDict([('EXTNAME', 'ANTENNA'), ('FILENAME', None)])
    hdr_coms = {
        'EXTNAME': 'name of binary table',
        'FILENAME': 'filename of "Antenna Log"',
        'TTYPE1': 'label for field 1',
        'TFORM1': 'data format of field 1',
        'TTYPE2': 'label for field 2',
        'TFORM2': 'data format of field 2',
        'TTYPE3': 'label for field 3',
        'TFORM3': 'data format of field 3',
        'TUNIT3': 'data unit of field 3',
        'TTYPE4': 'label for field 4',
        'TFORM4': 'data format of field 4',
        'TUNIT4': 'data unit of field 4',
        'TTYPE5': 'label for field 5',
        'TFORM5': 'data format of field 5',
        'TUNIT5': 'data unit of field 5',
        'TTYPE6': 'label for field 6',
        'TFORM6': 'data format of field 6',
        'TUNIT6': 'data unit of field 6',
        'TTYPE7': 'label for field 7',
        'TFORM7': 'data format of field 7',
        'TUNIT7': 'data unit of field 7',
        'TTYPE8': 'label for field 8',
        'TFORM8': 'data format of field 8',
        'TUNIT8': 'data unit of field 8'
    }
    col_vals = OrderedDict([
        ('time', None), ('scantype', None), ('az', None), ('el', None),
        ('ra', None), ('dec', None), ('az_center', None), ('el_center', None)
    ])
    col_form = OrderedDict([
        ('time', '26A'), ('scantype', '4A'), ('az', 'D'), ('el', 'D'),
        ('ra', 'D'), ('dec', 'D'), ('az_center', 'D'), ('el_center', 'D')
    ])
    col_unit = OrderedDict([
        ('time', None), ('scantype', None), ('az', 'deg'), ('el', 'deg'),
        ('ra', 'deg'), ('dec', 'deg'), ('az_center', 'deg'), ('el_center', 'deg')
    ])

    antenna_dict = {
        'hdr_vals': hdr_vals, 'hdr_coms': hdr_coms, 'col_vals': col_vals,
        'col_form': col_form, 'col_unit': col_unit
    }
    return antenna_dict

#---------------- READOUT
def readout_dict():
    hdr_vals = OrderedDict([('EXTNAME', 'READOUT'), ('FILENAME', None)])
    hdr_coms = {
        'EXTNAME': 'name of binary table',
        'FILENAME': 'filename which is readed for READOUT',
        'TTYPE1': 'label for field 1',
        'TFORM1': 'data format of field 1',
        'TTYPE2': 'label for field 2',
        'TFORM2': 'data format of field 2',
        'TTYPE3': 'label for field 3',
        'TFORM3': 'data format of field 3',
        'TUNIT3': 'data unit of field 3',
        'TTYPE4': 'label for field 4',
        'TFORM4': 'data format of field 4',
        'TUNIT4': 'data unit of field 4',
        'TTYPE5': 'label for field 5',
        'TFORM5': 'data format of field 5',
        'TUNIT5': 'data unit of field 5',
        'TTYPE6': 'label for field 6',
        'TFORM6': 'data format of field 6',
        'TUNIT6': 'data unit of field 6',
        'TTYPE7': 'label for field 7',
        'TFORM7': 'data format of field 7',
        'TUNIT7': 'data unit of field 7'
    }
    col_vals = OrderedDict([
        ('starttime', None), ('pixelid', None), ('amplitude', None), ('phase', None),
        ('line_phase', None), ('Tsignal', None), ('Psignal', None)
    ])
    col_form = OrderedDict([
        ('starttime', '26A'), ('pixelid', 'K'), ('Amplitude', '63D'), ('Phase', '63D'),
        ('Line_Phase', '63D'), ('Tsignal', '63D'), ('Psignal', '63D')
    ])
    col_unit = OrderedDict([
        ('starttime', None), ('pixelid', None), ('Amplitude', '??'), ('Phase', '??'),
        ('Line_Phase', '??'), ('Tsignal', 'K'), ('Psignal', 'W')
    ])

    readout_dict = {
        'hdr_vals': hdr_vals, 'hdr_coms': hdr_coms, 'col_vals': col_vals,
        'col_form': col_form, 'col_unit': col_unit
    }
    return readout_dict

#---------------- WEATHER
def weather_dict():
    hdr_vals = OrderedDict([('WEATHER', 'ANTENNA'), ('FILENAME', None)])
    hdr_coms = {
        'EXTNAME': 'name of binary table',
        'FILENAME': 'filename of "Weather Log"',
        'TTYPE1': 'label for field 1',
        'TFORM1': 'data format of field 1',
        'TTYPE2': 'label for field 2',
        'TFORM2': 'data format of field 2',
        'TUNIT2': 'data unit of field 2',
        'TTYPE3': 'label for field 3',
        'TFORM3': 'data format of field 3',
        'TUNIT3': 'data unit of field 3',
        'TTYPE4': 'label for field 4',
        'TFORM4': 'data format of field 4',
        'TUNIT4': 'data unit of field 4',
        'TTYPE5': 'label for field 5',
        'TFORM5': 'data format of field 5',
        'TUNIT5': 'data unit of field 5',
        'TTYPE6': 'label for field 6',
        'TFORM6': 'data format of field 6',
        'TUNIT6': 'data unit of field 6'
    }
    col_vals = OrderedDict([
        ('time', None), ('temperature', None), ('pressure', None),
        ('vapor-pressure', None), ('windspd', None), ('winddir', None)
    ])
    col_form = OrderedDict([
        ('time', '19A'), ('temperature', 'D'), ('pressure', 'D'),
        ('vapor-pressure', 'D'), ('windspd', 'D'), ('winddir', 'D')
    ])
    col_unit = OrderedDict([
        ('time', None), ('temperature', 'deg_C'), ('pressure', 'hPa'),
        ('vapor-pressure', 'hPa'), ('windspd', 'm/s'), ('winddir', 'deg')
    ])

    weather_dict = {
        'hdr_vals': hdr_vals, 'hdr_coms': hdr_coms, 'col_vals': col_vals,
        'col_form': col_form, 'col_unit': col_unit
    }
    return weather_dict

#-------------------------------- MAIN
if __name__ == '__main__':
#---------------- Create the dictionary for DFITS
    dfits_dict = {
        'obsinfo_dict': obsinfo_dict(), 'antenna_dict': antenna_dict(),
        'readout_dict': readout_dict(), 'weather_dict': weather_dict()
    }

#---------------- Write DFITS to the file (yaml)
    with PATH_DFITSDICT.open('w') as f:
        f.write(yaml.dump(dfits_dict, default_flow_style=False))
