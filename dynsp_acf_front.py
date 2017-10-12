import argparse
import create_dyn
from functions import *
import output
import combine_dyn
import acf
import secondary
import elevation
import warnings
import wavelength
import slice_dyn
import numpy as np



def standard_processing(archive, args, combined):
    full_arr, name, site, minFreq, maxFreq, mjd_start, mjd_end, archive, ra, dec = preparing_dyn(archive, args, combined)
    if not args.create:
        if args.slice!=[1,1]:
            full_arr, minFreq, maxFreq, mjd_start, mjd_end = slice_dyn.slice_dyn(full_arr, args,  minFreq, maxFreq, mjd_start, mjd_end)
            for i in range(0, np.shape(full_arr)[0]):
                pseudo_archive_name = archive + '_part_%s' % (i+1)
                standard_analysis(full_arr[i], args, name, site, minFreq[i], maxFreq[i], mjd_start[i], mjd_end[i], pseudo_archive_name, ra, dec) 
        else:        
            standard_analysis(full_arr, args, name, site, minFreq, maxFreq, mjd_start, mjd_end, archive, ra, dec)



def preparing_dyn(archive, args, combined):
    full_arr, name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec = create_dyn.get_dyn(archive, args, combined)

    if args.elevation:
        full_arr = elevation.correct_elevation(full_arr, name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec, archive)

    if args.range:
        full_arr, minFreq, maxFreq, mjd_start, mjd_end = create_dyn.set_range(full_arr, args, minFreq, maxFreq, mjd_start, mjd_end)

    if args.wavelength:
        full_arr = wavelength.wavelength_to_frequency(full_arr, minFreq, maxFreq)

    full_arr = np.nan_to_num(full_arr)
    if args.normal_t:
        for i in range(full_arr.shape[1]):
            maxi = np.max(full_arr[:,i])
            for j in range(full_arr.shape[0]):
                full_arr[j,i] = full_arr[j,i]/maxi
                if np.isnan(full_arr[j,i]):
                    full_arr[j,i] = 0

    #full_arr_masked = np.ma.masked_equal(full_arr, 0)
    if args.range:
        archive = archive + '_' + ','.join(args.range) 


    return full_arr, name, site, minFreq, maxFreq, mjd_start, mjd_end, archive, ra, dec



def standard_analysis(full_arr, args, name, site, minFreq, maxFreq, mjd_start, mjd_end, archive, ra, dec):

    if args.savedyn:
        # if args.range or args.slice:
        create_dyn.write_dyn(full_arr, archive, name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec)

    diff, length, mjd, mhzperbin, secperbin, ar_freq = acf.help_variables(full_arr, minFreq, maxFreq, mjd_start, mjd_end)
    full_arr, action = acf.remove_mean(full_arr, args, minFreq, diff)
    hwhm_freq, hwhm_freq_mod, freq_error, e_time, e_time_mod, time_error,curv, curv_error = (np.nan,)*8
    full_arr_masked = np.ma.masked_equal(full_arr, 0)
    modulation = np.ma.std(full_arr_masked) / np.ma.mean(full_arr_masked)
    if not args.acf:
        acf_norm, length, middle_f, middle_t, freqlag_ticks, midACF_freq, extrapo_x, extrapo_y, timelag_ticks, midACF_time, extrapo_x2, extrapo_y2, \
            hwhm_freq, freq_error, e_time, time_error= acf.analyse_acf(action, diff, length, args, funcs, fill_factor, minFreq)
        hwhm_freq_mod = np.sqrt(hwhm_freq**2 - mhzperbin**2)
        e_time_mod = np.sqrt(e_time**2 - (secperbin/3600)**2)

        print ("Calculated Values: \n freq_diss = %f +- %f MHz, time_diss = %f +- %f min" %
            (hwhm_freq_mod, freq_error, e_time_mod*60, time_error*60))
    
    if args.secondary:

        sec_spec, sec_axes, curv, curv_error, arc_strength, power_q1, power_q2  = secondary.compute_secondary(action, args, minFreq, maxFreq, diff, archive, mhzperbin, secperbin, name, mjd, length, ar_freq)
    # Plotting functions in matplotlib
    if args.write:
        with open("%s" % args.write, "a") as myfile:
            myfile.write("\n %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f"
                % (name, mjd, ar_freq, hwhm_freq, hwhm_freq_mod, freq_error, e_time*60, e_time_mod*60, time_error*60,curv, curv_error, modulation, arc_strength, power_q1, power_q2))
    if args.save or not args.view:
        output.plot_dyn(full_arr, args, length, diff, minFreq, maxFreq, ar_freq, name, mjd, archive)
    if not args.acf:
        output.plot_acf(full_arr, args, acf_norm, length, minFreq, maxFreq, middle_f, middle_t, 
            diff, freqlag_ticks, midACF_freq, extrapo_x, extrapo_y, timelag_ticks, midACF_time, extrapo_x2, extrapo_y2, archive, 
            hwhm_freq_mod, freq_error, e_time_mod, time_error, mhzperbin, secperbin, ar_freq, name, mjd, site)
    if args.secondary:
        output.plot_secondary(sec_spec, args, mhzperbin, secperbin, archive, ar_freq, name, diff, length, mjd, site, sec_axes, curv)

    output.show_plots(args)


parser = argparse.ArgumentParser(description='Select the archive.')
parser.add_argument('archive', nargs='+', help='The chosen archives or lists of archives (-l)')
parser.add_argument('-n','--normal_t', action='store_true', help='Normalize the dynamic spectrum in time.')
parser.add_argument('-m','--normal_f', action='store_true', help='Normalize the dynamic spectrum in frequency.')
parser.add_argument('-s', '--save', action='store_true', help='Save the images.')
parser.add_argument('-v', '--view', action='store_true', help='Do not view the plots.')
parser.add_argument('-t', '--trange', type=float, default=0, help='Sets the time range on acf plots (in minutes).')
parser.add_argument('-f', '--frange', type=float, default=0, help='Sets the frequency range on acf plots (in MHz).')
parser.add_argument('-w', '--write', type=str, help='Writes the parameters in a file with given path.')
parser.add_argument('-c', '--combine', action='store_true', help='Combines multiple archive to one dynamic spectrum.')
parser.add_argument('-l', '--list', action='store_true', help='Interprets the ar argument as .txt`s containing the archives.\
                                                                In combination with -c each line are the archives that get combined')
parser.add_argument('-p', '--secondary', action='store_true', help='Calculates the secondary spectrum')
parser.add_argument('-e', '--elevation', action='store_true', help='Corrects for the elevation')
parser.add_argument('-g', '--edges', action='store_true', help='Removes egdes from high frequency data.')
parser.add_argument('-q', '--wavelength', action='store_true', help='Converts the dynamic spectrum from frequencies to wavelengths')
parser.add_argument('-u', '--curvature', action='store_true', help='Computes the curvature of the scintillation arc.')
parser.add_argument('-r', '--range', type=lambda s: [item for item in s.split(',')], default=None,  
     help="Truncate the secondary spectrum to a certain range. Format: '-r min_freq, max_freq, start_hour, end_hour' \
     Frequencies in MHz, time in hours from the start of the observation")
parser.add_argument('-a', '--acf', action='store_true', help='Suppresses the acf calculation.')
parser.add_argument('-i', '--slice', nargs=2, type=int, default=[1,1], help="Slices the dynamic spectrum into several chunks.\
    This happens after using --range. Format: '-i number_of_chunks_in_freq, number_of_chunks_in_time'")
parser.add_argument('-z', '--gaussian_filter', nargs = 2, type = float, help='Applies a filter to the seconday spectrum, \
    Arguments define the size of the filter in percent in the axes' )
parser.add_argument('-x', '--curvature_new', action='store_true', help='Computes the curvature of the scintillation arc.')
parser.add_argument('-k', '--create', action='store_true', help='Only create the dyn spectrum files.')
parser.add_argument('-F', '--filling', type=float, default=0.4, help='Sets the filling factor.')
parser.add_argument('-H', '--hough', action='store_true', help='Computes the curvature of the scintillation arc via the Hough transform.')
parser.add_argument('-O', '--options', type=str, help='Stores the options for the Hough transform.')
parser.add_argument('-S', '--savedyn', action='store_true', help='Saves the dynamic spectrum after applying -r and -i.')
parser.add_argument('-E', '--extent', action='store_true', help='Saves the extent of the parabolic arc features in the secondary spectrum.')




warnings.simplefilter('ignore', RuntimeWarning)
warnings.simplefilter('ignore', UserWarning)

args = parser.parse_args()
# ar_path = str(os.path.abspath(args.ar)) 

funcs = {'func_f' : func_e_norm,
            'func_t' : func_g_norm}

# np.seterr(invalid='ignore')
# t_scale = args.tscale
# f_scale = args.fscale
fill_factor = args.filling
archives = args.archive
if args.list:
    archives = []
    for file in args.archive:
        print file
        with open(file) as l:
            archives.append(l.read().splitlines())
if args.combine:
    if not args.list:
        dyn_name = combine_dyn.stitch_dyn(archives, args)
        standard_processing(dyn_name, args, 1)
    else:
        for lis in archives:
            for row in lis:
                archives = row.split()
                dyn_name = combine_dyn.stitch_dyn(archives, args)
                standard_processing(dyn_name, args, 1)

else:
    for archive in archives:
        standard_processing(archive, args, 0)
