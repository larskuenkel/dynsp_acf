import numpy as np
from astropy.io import fits
import secondary
import argparse
import matplotlib.pyplot as plt
import create_dyn


parser = argparse.ArgumentParser(description='Commands for program')
parser.add_argument('file', help='The chosen fits file')

args = parser.parse_args()

file = args.file

full_arr, name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec = create_dyn.load_from_txt(file)
ini_shape = np.shape(full_arr)
full_arr = np.rot90(full_arr, 1)
full_arr = np.flipud(full_arr)
mid_freq = (maxFreq + minFreq) / 2.0
mid_mjd = (mjd_start + mjd_end) /2.0
bandwidth = maxFreq - minFreq
length = (mjd_end - mjd_start) * 86400.
freq_reso = bandwidth / float(ini_shape[0])
int_time = length /float(ini_shape[1])
ini_header = fits.Header()
ini_header['MJD'] = mid_mjd
ini_header['FREQ'] = mid_freq
ini_header['SOURCE'] = name
ini_header['ORIGIN'] = site
ini_header['BW'] = bandwidth
ini_header['T_INT'] = int_time
ini_header['TTOT'] = length
ini_header['CDELT1'] = freq_reso
ini_header['CDELT2'] = int_time
ini_header['NCHANS'] = ini_shape[0]
ini_header['NCHUNKS'] = ini_shape[1]
hdu = fits.PrimaryHDU(full_arr, header=ini_header)
hdulist = fits.HDUList([hdu])
# print fits.info(hdu)

# fits.setval(hdu, 'MJD', value=mid_mjd)
# print repr(hdulist[0].header)
print repr(hdu.header)

hdulist.writeto('%s_%s_%.2f_%.1f.fits' % (name, site, mid_mjd, mid_freq), overwrite='True')
