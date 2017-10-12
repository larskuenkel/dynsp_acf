import numpy as np
from astropy.io import fits
import secondary
import argparse
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description='Commands for program')
parser.add_argument('file', help='The chosen fits file')

args = parser.parse_args()

file = args.file

file_content = fits.open(file)
center_frequency = file_content[0].header['FREQ']
bandwidth =  file_content[0].header['BW']
minFreq = center_frequency - 0.5*bandwidth
maxFreq = center_frequency + 0.5*bandwidth
# name = file_content[0].header['SOURCE']
time = file_content[0].header['TTOT']
mjd_start = file_content[0].header['MJD']
mjd_end = file_content[0].header['MJD'] + time /(60.*60.*24.)
# site = file_content[0].header['ORIGIN']
print repr(file_content[0].header)

data = file_content[0].data
# data = np.rot90(data, 1)
data[data < 0] = 0
# plt.imshow(data, aspect='auto',  vmax=5, vmin=0)
# plt.colorbar()
# plt.show()
np.savetxt('%s_dyn.txt' % args.file, data, fmt='%.5f')

with open('%s_dyn.txt' % args.file, 'a') as myfile:
    myfile.write("\n # - - %f %f %f %f - -" % (minFreq, maxFreq, mjd_start, mjd_end))
