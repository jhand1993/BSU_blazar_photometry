import os
import glob
from astropy.io import fits
import numpy as np
import fitsio

import target_data
import config
import read_fits
import ref_cats

td = target_data.TargetData()

# All values from target_data are Byte data type.
target_name = td.target_data[:, 0].tolist()
target_RA = td.target_data[:, 1].tolist()
target_dec = td.target_data[:, 2].tolist()

# Converts Byte data in each list to string data type.
target_name = td.bytes_to_str(target_name)
target_RA = td.bytes_to_str(target_RA)
target_dec = td.bytes_to_str(target_dec)

RA_dict = td.target_dict(target_name, target_RA)
dec_dict = td.target_dict(target_name, target_dec)


def haversine(ra1, dec1, ra2, dec2, radius=1):
    """
    Haversine formula used to calculate great circle distance (geodesic between two sources on celestial sphere).
    :type radius: radius of sphere
    :param ra1: right ascension of first object in degrees
    :param dec1: declination of first object in degrees
    :param ra2: right ascension of second object in degrees
    :param dec2: declination of second object in degrees
    :param radius: Default to unity.  Distance from origin to spherical surface
    :return: distance in degrees
    """
    # Convert all necessary input parameters to radians
    rad_ra1 = np.pi * ra1 / 180
    rad_ra2 = np.pi * ra2 / 180
    rad_dec1 = np.pi * dec1 / 180
    rad_dec2 = np.pi * dec2 / 180

    delta_ra = abs(rad_ra2 - rad_ra1)
    delta_dec = abs(rad_dec2 - rad_dec1)
    # bla = np.sin(delta_dec / 2) ** 2
    # bloo = np.cos(rad_dec1) * np.cos(rad_dec2) * (np.sin(delta_ra / 2) ** 2)
    value = np.sqrt((np.sin(delta_dec / 2) ** 2) + np.cos(rad_dec1) * np.cos(rad_dec2) * (np.sin(delta_ra / 2) ** 2))
    return np.arcsin(value) * 2 * radius


def find_ref(ref_source, cat):
    """
    Returns list of reference stars in catalog using coordinates of reference stars in source_ref.
    :param ref_source: first column must be right ascension, second column must be declination (both in degrees)
    :param cat: input catalog from pipeline. Must be 'B', 'V', 'R' or 'I'.
    :return: [[row number, distance in degrees],...]
    """
    os.chdir(config.blazar_photometry)
    try:
        cat_data = fits.open(cat)
    except:
        raise
    cat_len = len(cat_data[2].data)
    ref_list = []
    for ref_star in range(len(ref_source)):
        closest = 1e10
        uid = float(np.copy(ref_source[ref_star][0]))
        ref_ra = np.copy(ref_source[ref_star][1])
        ref_dec = np.copy(ref_source[ref_star][2])
        for i in range(cat_len):
            cat_ra = np.copy(cat_data[2].data[i][33])
            cat_dec = np.copy(cat_data[2].data[i][34])
            delta = haversine(ref_ra, ref_dec, cat_ra, cat_dec)
            if delta < closest:
                closest = delta
                closest_source = i + 1
        # grabbing FLUX-AUTO and FLUXERR-AUTO from catalog file
        ref_flux = float(np.copy(cat_data[2].data[closest_source - 1][4]))
        ref_fluxerr = float(np.copy(cat_data[2].data[closest_source - 1][5]))
        ref_data = [uid, closest_source, ref_flux, ref_fluxerr]
        ref_list.append(ref_data)
    return ref_list


def find_source(cat):
    """
    Uses methods from target_data class to determine actual source in input catalog.
    :param cat: Input catalog from pipeline
    :return:
    """
    source_ra = float(td.coord_lookup(cat, RA_dict)) * 15.0
    source_dec = float(td.coord_lookup(cat, dec_dict))
    os.chdir(config.blazar_photometry)
    try:
        cat_data = fits.open(cat)
    except:
        raise
    cat_len = len(cat_data[2].data)
    closest = 1e10
    for k in range(cat_len):
        cat_ra = np.copy(cat_data[2].data[k][33])
        cat_dec = np.copy(cat_data[2].data[k][34])
        delta = haversine(source_ra, source_dec, cat_ra, cat_dec)
        if delta < closest:
            closest = delta
            closest_source = k + 1
    # grabbing FLUX-AUTO and FLUXERR-AUTO from catalog file
    source_flux = float(np.copy(cat_data[2].data[closest_source - 1][4]))
    source_fluxerr = float(np.copy(cat_data[2].data[closest_source - 1][5]))
    return [closest_source, source_flux, source_fluxerr]


def mag_fit(ref_source, cat, filt):
    """
    Uses find_ref and find_source to create a best fit line for x-axis being flux and y-axis being magnitude.
    :param ref_source: array containing reference stars
    :param cat: Input catalog from pipeline
    :param filt: Filter of stacked file.  Must be 'B', 'V', 'R' or 'I'.
    :return: [catalog name, source information, reference star information, linear fit coefficients,
              reference aperture flux, reference magnitude, flux error, magnitude error]
    """
    ref_list = find_ref(ref_source=ref_source, cat=cat)
    os.chdir(config.blazar_photometry)
    ref_flux = []
    ref_fluxerr = []
    source = find_source(cat=cat)

    if filt == 'B':
        ref_mag = ref_source[:, 3].tolist()
        ref_magerr = ref_source[:, 4].tolist()
    elif filt == 'V':
        ref_mag = ref_source[:, 5].tolist()
        ref_magerr = ref_source[:, 6].tolist()
    elif filt == 'R':
        ref_mag = ref_source[:, 7].tolist()
        ref_magerr = ref_source[:, 8].tolist()
    elif filt == 'I':
        ref_mag = ref_source[:, 9].tolist()
        ref_magerr = ref_source[:, 10].tolist()
    else:
        print('Invalid filter.  Must be B, V, R or I.')
        return False
    # Remove all values of 999.0 in ref_mag and ref_magerr
    if 999. in ref_mag:
        return False
    while 999. in ref_mag:
        ref_mag.remove(999.)
    while 999. in ref_magerr:
        ref_magerr.remove(999.)
    # Uses uid to populate correct flux for each ref source whose mag or magerr is not 999.0
    for s in range(len(ref_mag)):
        for t in ref_list:
            if s == t[0]:
                ref_flux.append(t[2])
                ref_fluxerr.append((t[3]))
    print(source)
    print(ref_list)
    print(ref_flux)
    print(np.log(ref_flux))
    lin_fit = np.polyfit(x=np.log10(ref_flux), y=ref_mag, deg=1)
    print(lin_fit)
    return [cat, source, ref_list, lin_fit, ref_flux, ref_mag, ref_fluxerr, ref_magerr]


def mag_err(flux, fit):
    top = fit[0] * np.log10(flux + np.sqrt(flux))
    bottom = fit[0] * np.log10(flux - np.sqrt(flux))
    delm = top - bottom
    return delm / 2

# User specifies filename information
os.chdir(config.blazar_photometry)
cat_list = glob.glob('*.cat')
# ref_source_input = input('Enter ref_source: ')
ref_source_input = 'ref_mrk501'
year_input = input('Enter year: ')
ref_source_data = np.copy(ref_cats.ref_mrk501)

# loop to create .txt files containing all data(verbose_file) and data for plotting(data_file)
for catalog in cat_list:
    os.chdir(config.blazar_photometry)
    head_info = read_fits.decode_fitshead(catalog)
    filter_input = read_fits.get_filter(head_info)
    verbose_file = open(ref_source_input + '_' + filter_input + '_' + year_input + '_verbose.txt', 'a')
    data_file = open(ref_source_input + '_' + filter_input + '_' + year_input + '.txt', 'a')
    cat_data_list = mag_fit(ref_source=ref_source_data, cat=catalog, filt=filter_input)

    if not cat_data_list:
        pass
    else:
        for x in cat_data_list:
            verbose_file.write(str(x))
            verbose_file.write(',')
        verbose_file.write('\n')
        verbose_file.close()

        cat_name = cat_data_list[0]
        source_info = cat_data_list[1]
        source_position = source_info[0]
        source_flux_data = source_info[1]
        source_fluxerr_data = source_info[2]
        linear_fit = cat_data_list[3]
        source_mag = linear_fit[0] * np.log10(source_flux_data) + linear_fit[1]
        # source_magerr = linear_fit[0] * source_fluxerr_data
        source_magerr = mag_err(source_flux_data, linear_fit)
        data_file.write(cat_name + ',' + str(filter_input) + ',' + str(source_mag) + ',' + str(source_magerr) + ',')
        data_file.write('\n')
        data_file.close()
