import numpy as np
from astropy.io import fits
import target_data
import glob
import config
import os
import read_fits
# import aplpy
# SExtractor catalog column numbers.  Note that columns in fits file are base 1 while in python they are given as
# base 0.
racol = 33
deccol = 34
fluxcol = 4
fluxerrcol = 5


def haversine(ra1, dec1, ra2, dec2, radius=1, deg=True):
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

    delta_ra = rad_ra2 - rad_ra1
    delta_dec = rad_dec2 - rad_dec1
    # bla = np.sin(delta_dec / 2) ** 2
    # bloo = np.cos(rad_dec1) * np.cos(rad_dec2) * (np.sin(delta_ra / 2) ** 2)
    value = np.sqrt((np.sin(delta_dec / 2) ** 2) + np.cos(rad_dec1) * np.cos(rad_dec2) * (np.sin(delta_ra / 2) ** 2))
    sigma = np.arcsin(value) * 2
    if deg:
        return sigma * radius * 180 / np.pi
    else:
        return sigma * radius


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
        cat_ra = np.copy(cat_data[2].data[k][racol])
        cat_dec = np.copy(cat_data[2].data[k][deccol])
        delta = haversine(source_ra, source_dec, cat_ra, cat_dec)
        if delta < closest:
            closest = delta
            closest_source = k + 1
    # grabbing FLUX-AUTO and FLUXERR-AUTO from catalog file
    source_flux = float(np.copy(cat_data[2].data[closest_source - 1][fluxcol]))
    source_fluxerr = float(np.copy(cat_data[2].data[closest_source - 1][fluxerrcol]))
    return [closest_source, source_flux, source_fluxerr]


def find_obj(chosen, cats):
    """
    searches through cats to find closest source to each source given in chosen.  Returns flags for duplicate finds
    or nonexistent sources.
    """
    bigarr = np.zeros((len(cats), len(chosen), 4), dtype=float)
    # bigarr[0, :, :] = np.copy(chosen)
    dup = np.empty((len(cats), len(chosen)), dtype=bool)
    exists = np.empty((len(cats), len(chosen)), dtype=bool)
    for cat in range(len(cats)):
        # print(cats[cat])
        for ref in range(len(chosen)):
            closest = 2 / 60  # if no target is within two arcminutes, no source found
            closest_source = -1
            ref_ra = np.copy(chosen[ref][2])
            ref_dec = np.copy(chosen[ref][3])
            for source in range(len(cats[cat])):
                cat_ra = np.copy(cats[cat][source][2])
                cat_dec = np.copy(cats[cat][source][3])
                delta = haversine(ref_ra, ref_dec, cat_ra, cat_dec)
                if delta < closest:
                    closest = delta
                    closest_source = source
            # check to see if source in cat is already in bigarr
            if cats[cat][closest_source] in bigarr[cat, :ref, :]:
                # print(cats[cat][closest_source], bigarr[cat, :ref, 0])
                dup[cat, ref] = False
            else:
                dup[cat, ref] = True
            # mark exists as false if no closest source is found for a catalog
            if closest_source == -1:
                exists[cat, ref] = False
                bigarr[cat, ref, :] = np.array([-1., -1., -1., -1.])
            else:
                exists[cat, ref] = True
                bigarr[cat, ref, :] = np.copy(cats[cat][closest_source])
    return bigarr, dup, exists

# instantiate TargetData class to find RA and Dec for obj
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

obj_ra = float(RA_dict[config.obj.lower()]) * 15.
obj_dec = float(dec_dict[config.obj.lower()])

print(obj_ra, obj_dec)

# grab all catalogs with correct object
os.chdir(config.catalogs)
catalogs = np.asarray(glob.glob('*.cat'), dtype=str)
good_cats = np.empty(len(catalogs), dtype=bool)
for i in range(len(catalogs)):
    header = read_fits.decode_fitshead(catalogs[i])
    obj = str(read_fits.get_info('OBJECT', header))
    obj = obj.strip()
    obj = obj.replace('_', '')
    obj = obj.replace(' ', '')
    obj = obj.lower()
    if obj == config.obj.lower():
        good_cats[i] = True
    else:
        good_cats[i] = False
# catalogs = catalogs[good_cats]
mega = []

# create 3D array containing source number, flux, ra and dec for each catalog.
# fig1 = aplpy.FITSFigure(catalogs[0])
for i in catalogs:
    header = read_fits.decode_fitshead(i)
    cat_data = fits.open(i)
    fluxlist = [s[fluxcol] for s in cat_data[2].data]
    sid = np.asarray([s[0] for s in cat_data[2].data])
    sflux = np.asarray(fluxlist / np.sum(fluxlist))
    sra = np.asarray([s[racol] for s in cat_data[2].data])
    sdec = np.asarray([s[deccol] for s in cat_data[2].data])
    # first find and remove object of interest from arrays
    source_ra = float(td.coord_lookup(i, RA_dict)) * 15.0
    source_dec = float(td.coord_lookup(i, dec_dict))
    cat_len = len(cat_data[2].data)
    closest = 1e10
    for k in range(cat_len):
        cat_ra = np.copy(cat_data[2].data[k][racol])
        cat_dec = np.copy(cat_data[2].data[k][deccol])
        delta = haversine(source_ra, source_dec, cat_ra, cat_dec)
        if delta < closest:
            # print(source_ra, cat_ra, source_dec, cat_dec)
            closest = delta
            closest_source = k + 1
    # create mask of all sources that are not object of interest
    notobject = np.asarray([int(x) != int(closest_source-1) for x in sid])
    sid = sid[notobject]
    sflux = sflux[notobject]
    sra = sra[notobject]
    sdec = sdec[notobject]
    srow = np.stack((sid, sflux, sra, sdec), axis=1)
    # create mask of sources within radius r specified in config.py
    havs = [haversine(s[2], s[3], obj_ra, obj_dec) for s in srow]
    print(np.mean(havs), min(havs), max(havs))
    sgood = np.asarray([abs(haversine(s[2], s[3], obj_ra, obj_dec)) < config.r for s in srow])
    # print(len(srow))
    srow = srow[sgood]
    # print(len(srow))
    mega.append(srow)

# clear out sources with flux less than cutoff
cutoff = 1e-3
for i in range(len(mega)):
    cat = mega[i]
    catmask = np.asarray([x > cutoff for x in cat[:, 1]])
    cat = cat[catmask]
    # mega[i] = cat[catmask]
    mega[i] = cat
    # print(np.mean(cat[:, 1]), len(cat[:, 1]))
# find catalog with fewest sources
catcount = len(mega)
shortestcatval = 1e10
shortestcat = -1
for i in range(catcount):
    if len(mega[i]) < shortestcatval:
        shortestcatval = len(mega[i])
        shortestcat = i
shortestcat = mega[shortestcat]
shortlen = len(shortestcat)
# contract all other catalogs to size of shortest catalog:
for i in range(catcount):
    mega[i] = np.copy(mega[i][:shortlen])
megaarr = np.asarray(mega)
# print(np.shape(megaarr))
# sort megaarr and yield flag matrices
megasort, dup, missing = find_obj(shortestcat, megaarr)
# print(megasort)
print(shortlen)
stdcut = 1e-4
refsources = []
for i in range(shortlen):
    objids = megasort[:, i, 0]
    objfluxes = megasort[:, i, 1]
    # print(objfluxes)
    truthmask = np.logical_and(dup[:, i], missing[:, i])
    # print(truthmask)
    objids = objids[truthmask]
    # print(objids)
    objfluxes = objfluxes[truthmask]
    # print(objfluxes)
    objfluxstd = np.std(objfluxes)
    # newline = np.copy(megasort[:, i, :])
    if objfluxstd <= stdcut:
        # newline = np.append(newline, objfluxstd)
        # print(newline)
        refsources.append(objfluxstd)
    print(objfluxstd)
print(len(refsources))
