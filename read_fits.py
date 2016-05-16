import os
import fitsio

"""
Although astropy.io is a much more convenient way to read and utilize fits data, the nature of the stacked headers
from SWarp result in some data not being contained in one of the fits extensions of the LDAC catalogs from SExtractor.
These functions simply read the text from the fits files themselves (not the compiled information) and returns a list
of all header information.  Three functions can be used to extract DATE-OBS, MJD-OBS and FILTER, all of which are not
accessible when using astropy.io.
"""

# valid filters from data
valid_filters = ['B', 'V', 'R', 'I']


def decode_fitshead(file1, cwd=None):
    """
    Decodes header of fits file.
    :param file1: fits file whose header will be read and decoded
    :param cwd: specify cwd
    :return: list of header information from fits
    """
    if cwd:
        os.chdir(cwd)
    fitsfile = fitsio.read(file1)
    header_list = []
    for i in fitsfile:
        for j in i:
            for k in j:
                header_list.append(bytes.decode(k, 'utf-8'))
    return header_list


def get_dateobs(list1):
    # retrieves DATE-OBS from fits header
    for k in list1:
        if 'DATE-OBS' in k:
            first = k.split('=', 1)[1]
            second = first.split('/', 1)[0]
            third = second.replace('\'', '')
            dateobs = third.strip()
            return dateobs.split('T', 1)[0]
    return False


def get_filter(list1):
    # retrieves FILTER from fits header
    for k in list1:
        if 'FILTER' in k:
            first = k.split('=', 1)[1]
            second = first.split('/', 1)[0]
            third = second.replace('\'', '')
            fitsfilter = third.strip()
            for filt in valid_filters:
                if filt == fitsfilter:
                    return filt
    return False


def get_mjd(list1):
    # retrieves MJD-OBS from fits header
    for k in list1:
        if 'MJD-OBS' in k:
            first = k.split('=', 1)[1]
            second = first.split('/', 1)[0]
            third = second.replace('\'', '')
            mjd = third.strip()
            return mjd
    return False
