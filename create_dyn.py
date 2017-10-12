import numpy as np
import psrchive as psr
import subprocess
import matplotlib.pyplot as plt
import os


def dyn_txt(archive, args):
    arch = psr.Archive_load(archive) #load a certain archive
    minFreq = arch.get_centre_frequency() -  (arch.get_bandwidth()/2)
    maxFreq = arch.get_centre_frequency() +  (arch.get_bandwidth()/2)
    arfreq = arch.get_centre_frequency()
    length = arch.integration_length()
    name = arch.get_source()
    mjd_start = float(arch.start_time().strtempo())
    mjd_end = float(arch.end_time().strtempo())
    diff = maxFreq-minFreq
    coord = arch.get_coordinates()
    ra = coord.ra().getHMS()
    dec = coord.dec().getHMS()
    # print dir(coord)
    ra, dec = coord.getHMSDMS().split(' ')
    # print ra
    # print dec
    # print coord.getHMSDMS()
    # print dir(coord.getRaDec())
    site = arch.get_telescope()
    bpass = ''
    # if args.bandpass:
    #     bpass = '_calbp.ar'
    print 'New %s%s_dyn.txt is being created' % (archive, bpass)
    # if args.bandpass:
    #     print 'Calibrating the bandpass.'
    #     calibrate_bandpass(archive)
    subprocess.call('/usr/share/psrchive-hacked/bin/dynamic_spectra -o %s%s_dyn.txt -f matlab %s%s' % (archive, bpass, archive, bpass), shell=True)
    # if args.remove:
    #     if args.bandpass:
    #         os.remove('%s%s' % (archive, bpass))
    with open("%s%s_dyn.txt" % (archive, bpass), "a") as myfile:
        myfile.write("\n # %s %s %f %f %f %f %s %s" % (name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec))


# def calibrate_bandpass(archive):
#     # Calibrates the bandpass. Based on:
#     # https://github.com/sosl/public_codes/blob/master/pulsar/bandpass/bandpass_correction.py
#     ar = psr.Archive_load(archive)
#     ar.pscrunch()
#     ar.tscrunch()
#     subint = ar.get_Integration(0)
#     (bl_mean, bl_var) = subint.baseline_stats()
#     bl_mean = bl_mean.squeeze()
#     bl_var = bl_var.squeeze()
#     non_zeroes = np.where(bl_mean != 0.0)
#     min_freq = ar.get_Profile(0, 0, 0).get_centre_frequency()
#     max_freq = ar.get_Profile(0, 0, ar.get_nchan()-1).get_centre_frequency()
#     freqs = np.linspace(min_freq, max_freq, ar.get_nchan())
#     # plt.clf()
#     # fig1 = plt.plot(freqs[non_zeroes],bl_mean[non_zeroes])
#     # xlab = plt.xlabel('frequency [MHz]')
#     # ylab = plt.ylabel('power [arbitrary]')
#     # plt.savefig(archive+"_bandpass.png")
#     ar = psr.Archive_load(archive)
#     ar.remove_baseline()
#     bl_mean_avg = np.average(bl_mean[non_zeroes])
#     for isub in range(ar.get_nsubint()):
#         for ipol in range(ar.get_npol()):
#             for ichan in range(ar.get_nchan()):
#                 prof = ar.get_Profile(isub, ipol, ichan) 
#                 if ichan in non_zeroes[0]:
#                     prof.scale(bl_mean_avg / bl_mean[ichan])
#                 else:
#                     prof.set_weight(0.0)
#     ar.unload('%s_calbp.ar' % archive)


def get_dyn(arg, args, combined):# bandpass=False, remove=False):
    if arg.endswith('.txt'):
        full_arr, name,site,  minFreq, maxFreq, mjd_start, mjd_end, ra, dec = load_from_txt(arg) 
        print '%s is loaded.' % arg
    else:
        try:
            full_arr, name,site,  minFreq, maxFreq, mjd_start, mjd_end, ra, dec = load_from_txt('%s_dyn.txt' % (arg))
            print 'Existing %s_dyn.txt is used' % (arg)
        except (IOError, ValueError):
            dyn_txt(arg, args)
            full_arr, name,site,  minFreq, maxFreq, mjd_start, mjd_end, ra, dec = load_from_txt('%s_dyn.txt' % (arg))
    if not combined:
        if args.edges:
            # full_arr[0,:] = 0
            # full_arr[-1,:] = 0
            # full_arr[-2,:] = 0
            full_arr = remove_edges(full_arr)
    return  full_arr, name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec


def load_from_txt(txt):
    with open(txt, 'r') as myfile:
        lines = myfile.readlines()
    meta_data = lines[-1]
    _, _, name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec = meta_data.split(' ')
    mjd_start = float(mjd_start)
    mjd_end = float(mjd_end)
    minFreq = float(minFreq)
    maxFreq = float(maxFreq)
    full_arr = np.genfromtxt('%s' % txt)
    full_arr = np.rot90(full_arr, 3)
    full_arr = np.fliplr(full_arr)
    return full_arr, name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec


def remove_edges(full_arr):
    mean_diff = []
    for j in range(full_arr.shape[0] - 1 ):
        mean_diff.append(np.abs(np.ma.mean(full_arr[j,:]) - np.ma.mean(full_arr[j+1,:])))
    mean_diff_std = np.std(mean_diff) * 2
    for channel in (0, 1, -1, -2):
        if mean_diff[channel] > mean_diff_std:
            # print channel
            full_arr[channel, :] = 0
    return full_arr


def set_range(full_arr, args, minFreq, maxFreq, mjd_start, mjd_end):
    hour_length = (mjd_end - mjd_start) * 24.0
    ini_shape = np.shape(full_arr)
    freqs = np.linspace(minFreq, maxFreq, ini_shape[0])
    times = np.linspace(0, hour_length, ini_shape[1])
    range_parameter = np.asarray(args.range)
    def find_nearest(array,value):
        near_index = None
        if value:
            value = float(value)
            near_index = (np.abs(array-value)).argmin()
        return near_index
    near_minfreq = find_nearest(freqs, range_parameter[0])
    near_maxfreq = find_nearest(freqs, range_parameter[1])
    near_hourstart = find_nearest(times, range_parameter[2])
    near_hourend = find_nearest(times, range_parameter[3])

    full_arr = full_arr[near_minfreq: near_maxfreq, near_hourstart: near_hourend]


    if not near_minfreq == None:
        minFreq = freqs[near_minfreq]
    if not near_maxfreq == None:
        maxFreq = freqs[near_maxfreq]
    if not near_hourstart == None:
        mjd_start = mjd_start + times[near_hourstart] / 24.0
    if not near_hourend == None:
        mjd_end = mjd_start + (times[near_hourend] - times[near_hourstart]) / 24.0

    return full_arr, minFreq, maxFreq, mjd_start, mjd_end


def write_dyn(full_arr, full_name, name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec):
    full_arr = np.fliplr(full_arr)
    full_arr = np.rot90(full_arr, 1)
    full_arr = np.ma.filled(full_arr, 0)
    if not full_name.endswith('.txt'):
        full_name += '.txt'
    np.savetxt('%s' % full_name, full_arr,  fmt='%.5f')
    with open('%s' % full_name, "a") as myfile:
        myfile.write("\n # %s %s %f %f %f %f %s %s" % (name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec))