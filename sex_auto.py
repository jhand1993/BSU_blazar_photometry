import subprocess as sp
import os
import glob

# Set current working directory to /Users/jaredhand/Documents/photometry/data_in/
os.chdir('/Users/jaredhand/Documents/photometry/data_in/')
# str for sextractor script
script = 'sex -CATALOG_NAME %s %s'


def get_files():
    """
    Decompress all FIT.gz files is they exists and returns list of all fit/FIT files in active directory.

    :return: fits+FITS as a list
    """
    try:
        if glob.glob('*.gz*'):
            sp.check_call('gunzip *FIT.gz', shell=True)
    except sp.CalledProcessError:
        raise

    fits = glob.glob('*.fit*')
    FITS = glob.glob('*.FIT*')
    return fits + FITS

files = get_files()


def rm_spaces(list1):
    """
    Removes all whitespace from members of a list (assuming list contains only str objects).

    :param list1: input list where are whitespace is replaced with '_'
    :return: modified list1
    """
    for i in list1:
        if ' ' in i:
            j = i.replace(' ', '_')
            os.rename(i, j)
            i = j
    return list1

files_nospaces = rm_spaces(list1=files)


def sex_call(list1):
    """
    Runs sextractor on all members of input list1.

    :param list1: list of fits/FITS files in cwd
    :return: False if subprocess fails, otherwise returns True
    """
    print(os.getcwd())
    if os.getcwd() != '/Users/jaredhand/Documents/photometry/data_in':
        print('Check current directory.')
        return False
    elif not list1:
        print('No .fit or .FIT files in /data_in folder.')
        return False

    for i in list1:
        cat_name = i.split('.')[0] + '.cat'
        try:
            sp.check_call(script % (cat_name, i), shell=True)
        except:
            raise
    return True

print(os.getcwd())
sex_call(list1=files_nospaces)





