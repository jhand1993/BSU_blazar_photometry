import subprocess as sp
import os
import glob



os.chdir('/Users/jaredhand/Documents/photometry/data_in/')
sp.check_call('')
files = glob.glob('*')

fits = glob.glob('*.fit*')
FITS = glob.glob('*.FIT*')

files = fits + FITS

sp.Popen()


