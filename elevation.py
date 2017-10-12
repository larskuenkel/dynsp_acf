import numpy as np
from astropy.time import Time
import subprocess
import matplotlib.pyplot as plt
import create_dyn
# from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from scipy.signal import medfilt
from scipy.signal import wiener
import os


def correct_elevation(full_arr, name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec, archive):
    FNULL = open(os.devnull, 'w')
    try:
        full_arr, _, _, _, _, _, _, _, _ = create_dyn.load_from_txt('%s_calibE_dyn.txt' % archive)
        print 'Existing calibrated %s_calibE_dyn.txt is used' % (archive)
    except:
        print 'Correcting the elevation.'
        start_date = mjd2datestring(mjd_start)
        end_date = mjd2datestring(mjd_end)
        # print start_date
        # print end_date
        # print full_arr.shape
        times_mjd = np.linspace(mjd_start, mjd_end, num = full_arr.shape[1])
        # print times_mjd.shape
        # print times_mjd
        times_str = mjd2datestring(times_mjd)
        # print times_str

        try:
            elevations = np.genfromtxt('elev_%s_%s_%s_%s.txt' % (name, site, mjd_start, mjd_end))
            print ('Uses elev_%s_%s_%s_%s.txt as elevations.' % (name, site, mjd_start, mjd_end))
        except:
            corrections = []
            elevations = []
            time_number = len(times_str)
            i = 1
            for time in times_str:
                time_dummy = str(time).split(' ')
                time_mod = time_dummy[0] + 'T' +time_dummy[1]
                output_ini = subprocess.Popen('/home/lkuenkel/pulsar_up/azza.pl -t %s -ra %s -dec %s -site %s' % (time_mod, ra, dec, site), shell=True, stdout=subprocess.PIPE, stderr=FNULL)

                if not i%50:
                    print '%s / %s' % (i, time_number)
                i += 1
                output = output_ini.stdout.read().split()
                n = 0
                for string in output:
                    if string =='EL':
                        m = n
                    n += 1
                el = output[m+2]
                elevation = float(el) * np.pi / 180.0
                # print elevation
                # correction = np.sin(elevation)**-1.39
                # print correction
                elevations.append(elevation)
                # corrections.append(correction)
            np.savetxt('elev_%s_%s_%s_%s.txt' % (name, site, mjd_start, mjd_end), elevations)
        means = []

        full_arr_masked = np.ma.masked_equal(full_arr, 0)
        full_arr_masked = np.ma.fix_invalid(full_arr_masked, fill_value=0)
        for j in range(full_arr.shape[1]):
            means.append(np.ma.mean(full_arr_masked[:,j]))
        # print means
        # print means
        mean_condition = np.asarray(np.nonzero(np.asarray(means) > 0)).squeeze()
        # print mean_condition
        elevs_cleaned = []
        means_cleaned = []
        x_cleaned = []
        x_values = np.linspace(mjd_start, mjd_end, num = full_arr.shape[1])
        for nonzero_mean in mean_condition:
            elevs_cleaned.append(elevations[nonzero_mean])
            means_cleaned.append(means[nonzero_mean])
            x_cleaned.append(x_values[nonzero_mean])
        # medians = np.ma.median(full_arr_masked, axis=0)
        # print np.shape(medians)
        # print np.shape(full_arr)
        # x_values = np.linspace(mjd_start, mjd_end, num = full_arr.shape[1])
        # elevs_cleaned = elevations[~medians.mask]
        # x_cleaned = x_values[~medians.mask]
        # medians_cleaned = medians[~medians.mask]
        # means_cleaned = medians_cleaned
        # mean = np.ma.mean(full_arr_masked)
        # print means
        # means_filt = medfilt(means_cleaned, 15)
        # means_filt = wiener(means_cleaned)
        means_filt = means_cleaned

        def corr_f2(x, a, b):
            return np.sin(x)**-a * b
        inv_mean = np.divide(1, means_filt)
        inv_mean = np.nan_to_num(inv_mean)
        # print inv_mean
        # print means
        # print elevs_cleaned
        # print inv_mean
        opt, cov = curve_fit(corr_f2,elevs_cleaned,inv_mean, p0=(1, 0.3))
        corr2 = []
        fit2 = []
        for elev in elevations:
            corr2.append(np.sin(elev) ** -opt[0] )
            fit2.append(np.sin(elev) ** -opt[0] * opt[1])
        print 'Correction factor for elevation: %s' % -opt[0]
        if -opt[0]>0:
            print 'Value is expected to be negative.'
        x_values = np.linspace(mjd_start, mjd_end, num = full_arr.shape[1])
        full_arr = np.multiply(full_arr, corr2)
        # mean_mod = np.multiply(means, corr2)
        # mean_mod = np.multiply(medians, corr2)
        full_arr_rot = np.fliplr(full_arr)  
        full_arr_rot = np.rot90(full_arr_rot, 1)
        full_arr = np.ma.filled(full_arr, 0)
        np.savetxt('%s_calibE_dyn.txt' % archive, full_arr_rot)



        # plt.plot(x_cleaned, inv_mean, x_values, fit2)
        # plt.show()
        with open('%s_calibE_dyn.txt' % archive, "a") as myfile:
            myfile.write("\n # %s %s %f %f %f %f %s %s" % (name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec))


    return full_arr



def mjd2datestring(mjd):
    time_object = Time(mjd, format='mjd', scale='utc')
    time_object.format = 'iso'
    return time_object

def corr_func(a, mean, elev):
    return np.sin(elev)**-a * mean - b
