import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
import subprocess
import os
import functions
from matplotlib.patches import Ellipse
from scipy.signal import medfilt
#import gaussian_filter

def new_func(y,a):
    return (y / a) ** 0.5



def compute_secondary(action, args, minFreq, maxFreq, diff, archive, mhzperbin, secperbin, name, mjd, length, ar_freq):
    #Computes the secondary spectrum
    #action = acf.remove_mean(action, args, minFreq, diff)
    hour_length = length / 3600
    action = np.ma.filled(action, np.ma.median(action))
    # plt.imshow(action, aspect='auto')
    # plt.colorbar()
    # plt.show()
    # action = action - np.mean(action)
    axis_1 = next_greater_power_of_2(action.shape[0])
    axis_2 = next_greater_power_of_2(action.shape[1])

    #padded = np.zeros((axis_1,axis_2))
    padded = np.zeros(np.shape(action))
    padded[:action.shape[0],:action.shape[1]] = action[:,:]
    sec_initial = np.fft.fft2(padded)
    # print np.shape(sec_initial)
    sec_i_abs = (np.abs(sec_initial))**2
    sec_shifted = np.fft.fftshift(sec_i_abs)
    sec_median = np.median(sec_shifted)
    if args.gaussian_filter:
        #if not args.curvature_new:
        shape1, shape2 = np.shape(sec_shifted)
        filter1, filter2 = args.gaussian_filter
        filter1 = int(filter1 * shape1 *.005)
        filter2 = int(filter2 * shape2 *.01)
        sec_shifted = gaussian_filter(sec_shifted, sigma=(filter1,filter2))


    shape_f = int(np.shape(sec_shifted)[0]/2.0) 
    shape_f_full = np.shape(sec_shifted)[0]
    shape_t = np.shape(sec_shifted)[1]
    t_end = int(1/7.0*shape_t)
    x_vals = np.linspace(-10, 10, shape_f_full)
    outer_medians = medfilt(np.median(sec_shifted[:,0:t_end], axis = 1), kernel_size=41)
    opt, cov = curve_fit(functions.linear_sym, x_vals, outer_medians)
    med_correction = - functions.linear_sym(x_vals, *opt) + opt[1] 
    for j in range(np.shape(sec_shifted)[0]):
        sec_shifted[j,:] = sec_shifted[j,:] + med_correction[j]
    sec_log = 10. * np.log10(sec_shifted / np.max(sec_shifted))      
    sec_cut = sec_log[sec_log.shape[0]/2:sec_log.shape[0],:]
    # plt.plot(outer_medians)
    # plt.plot(functions.linear_sym(x_vals, *opt))
    # plt.plot(np.median(sec_shifted[:,0:t_end], axis = 1))

    global_median = np.median(sec_shifted)
    power_q1 = np.mean(sec_shifted[shape_f+30:shape_f+shape_f*3/4,shape_t/4:shape_t/2-30]) - global_median
    power_q2 = np.mean(sec_shifted[shape_f+30:shape_f+shape_f*3/4,shape_t/2+30:shape_t*3/4]) - global_median
    print power_q1
    print power_q2
    # plt.show()

    sec_final = sec_cut
    # np.savetxt('%ss.txt'%minFreq, sec_final)
    sec_axes = secondary_axes(sec_cut, mhzperbin, secperbin)
    a_value, a_error, arc_strength = (np.nan,)*3

    if args.extent:
        sec_means = np.mean(sec_final[30:np.shape(sec_final)[1]/2,:], axis=0)
        sec_means = medfilt(sec_means, kernel_size=13)
        counter_1 = 0
        f_start = int(2/7.0*sec_final.shape[0])
        f_end = int(5/7.0*sec_final.shape[0])
        med_reduced = np.median(sec_final[f_start:f_end,0:t_end])
        mid_shape = len(sec_means)/2
        breaking_1 = (np.max(sec_means[mid_shape/3:mid_shape-10])-np.mean(sec_means[mid_shape/3:mid_shape-10])) *0.3 + np.mean(sec_means[mid_shape/3:mid_shape-10])
        breaking_2 = (np.max(sec_means[mid_shape+10:-mid_shape/3])-np.mean(sec_means[mid_shape+10:-mid_shape/3])) *0.3 + np.mean(sec_means[mid_shape+10:-mid_shape/3])
        print breaking_1, breaking_2
        while 1:
            if counter_1 == len(sec_means)/2:
                counter_1 = np.nan
                break
            counter_1 += 1
            if sec_means[counter_1] >= breaking_1:
                break
        counter_2 = 0
        while 1:
            if counter_2 == len(sec_means)/2:
                counter_2 = np.nan
                break
            counter_2 += 1
            if sec_means[-counter_2] >= breaking_2:
                break
        print counter_1, counter_2
        complete_axis = np.shape(sec_final)[1] /2
        conv_factor = 50 / float(complete_axis)
        q1_extent = (complete_axis - counter_1 ) * conv_factor
        q2_extent = (complete_axis - counter_2 ) * conv_factor
        print q1_extent, q2_extent
        with open('extents.txt', "a") as myfile:
            myfile.write("\n %s %f %f" % (mjd, q1_extent, q2_extent))
        plt.plot(sec_means)
        plt.show()

    if args.curvature_new:

        a_value, error = new_curvature(sec_log, args, name, mjd, hour_length, ar_freq, diff)

    if args.hough:
        # med_freq = np.median(sec_final[:,0:100], axis=1)
        # med_time = np.median(sec_final, axis=0)
        # shape_f = np.shape(sec_final)[0]
        # shape_t = np.shape(sec_final)[1]
        # f_start = int(3/7.0*shape_f)
        # f_end = int(4/7.0*shape_f)
        # t_end = int(1/7.0*shape_f)
        # x_vals = np.arange(0, shape_f)
        # opt, cov = curve_fit(functions.linear, x_vals, np.median(sec_final[:,0:t_end], axis = 1))
        # med_correction = - functions.linear(x_vals, *opt) + opt[1] 
        # for j in range(shape_f):
        #     sec_final[j,:] = sec_final[j,:] + med_correction[j]
        # sec_final = sec_final[:,:]
        f_start = int(2/7.0*sec_final.shape[0])
        f_end = int(5/7.0*sec_final.shape[0])
        std_reduced = np.std(sec_final[f_start:f_end,0:t_end])
        # print std_reduced
        # med_reduced = np.median(sec_final[f_start:f_end,0:t_end])
        med_reduced = global_median
        curvature_command = find_command(sec_axes, sec_shifted, archive, args)
        print curvature_command
        if not os.path.isfile('./' + archive + '_sec.txt'):
            print 'Writing secondary spectrum to disk'
            rewrite_data(sec_shifted, archive)
        output_ini = subprocess.Popen('%s' % (curvature_command), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = output_ini.stdout.read()
        output_err = output_ini.stderr.read()
        print output
        print output_err
        hough = np.genfromtxt('./%s_hough.txt' % archive)
        hough_x = hough[:,0]
        hough_y = np.divide(hough[:,1], hough[:,2])
        # # y_values = y_values / np.max(y_values)
        # # print y_values, np.shape(y_values)
        max_pos = np.argmax(hough_y)
        max_value = hough_y[max_pos]
        arc_strength = max_value - med_reduced
        a_value_ini = float(output.split()[3])
        break_point = max_value - (max_value - med_reduced) /20.0
        # break_point = max_value - 0.1
        counter_upper = 0
        # print break_point
        # print max_value
        while 1:
            # print counter_upper
            # print hough_y[max_pos + counter_upper]
            counter_upper += 1
            if max_pos + counter_upper == len(hough_x):
                counter_upper += -1
                break
            # print counter_upper
            # print hough_y[max_value + counter_upper]
            if hough_y[max_pos + counter_upper] <= break_point:
                break
        diff_upper = hough_x[max_pos + counter_upper] - hough_x[max_pos]
        counter_lower = 0
        while 1:
            counter_lower += 1
            if max_pos - counter_lower == -1:
                counter_lower += -1
                break
            if hough_y[max_pos - counter_lower] <= break_point:
                break
        diff_lower = hough_x[max_pos] - hough_x[max_pos - counter_upper]
        max_error = max(diff_upper, diff_lower)

        a_value = a_value_ini
        a_error = max_error
        print a_value, a_error
        # med_freq2 = np.median(sec_final[:,0:100], axis=1)
        # med_std = np.std(sec_final[:,0:100], axis=1)
        # plt.plot(med_freq2)
        # plt.show()
        # plt.plot(med_std)
        # plt.show()
        sec_med = np.median(sec_final)
        title_data = '%s \n %s MJD:%.2f  Duration:%.2fh Frequency:%.2fMHz  Bandwith:%.2fMHz' % (args.options, name, mjd, length/ 3600.0, ar_freq, diff)
        plt.close('all')
        hough_plot = plt.figure(5, figsize=(12,6))
        # ax_big = hough_plot.add_subplot(111,frameon=False)
        hough_plot.suptitle('%s' % title_data, fontsize=10)
        # ax_big.set_yticklabels('')
        # ax_big.set_xticklabels('')
        ax1 = hough_plot.add_subplot(121)
        ax1.plot(hough_x, hough_y)
        ax1.axhline(med_reduced, color='r')
        ax1.axhline(break_point, color='r')
        if med_reduced - 0.5 <= max_value +0.5:
            ax1.set_ylim(med_reduced - 100,max_value +100)
        ax1.set_title('%s %s'%(a_value, a_error))
        ax2 = hough_plot.add_subplot(122)
        sec_im = ax2.imshow(sec_final, cmap='binary', aspect='auto', extent=[sec_axes[2], sec_axes[3], 0, sec_axes[1]], vmax=sec_med+6, vmin=sec_med+1, origin='lower')
        ax2.autoscale(False)
        x_values = np.linspace(sec_axes[2], sec_axes[3], np.shape(sec_final)[1])
        y_values = functions.parabola(x_values, a_value)
        file = open('dummy.gpi', 'r')
        file2 = open('dummy2.gpi', 'r')
        dummy_lines = file.readlines()
        dummy2_lines = file2.readlines()
        x_dist = float(dummy2_lines[4].split()[-4].strip('('))
        y_dist = float(dummy2_lines[4].split()[-2])
        x_values2 = (x_values + x_dist)
        y_values2 = (y_values - y_dist)
        x_values3 = (x_values - x_dist)
        y_values3 = (y_values + y_dist)
        x_values2_filtered = x_values2[x_values2>=0]
        y_values2_filtered = y_values2[x_values2>=0]
        x_values3_filtered = x_values3[x_values3>=0]
        y_values3_filtered = y_values3[x_values3>=0]
        vertical = float(dummy_lines[27].split()[0])
        horizontal = float(dummy_lines[36].split()[1])
        ell_height = 2*float(dummy_lines[20].split('*')[0])
        ell_width = 2*float(dummy_lines[20].split('/')[1].split(')')[0])
        ax2.axvline(vertical, color='r')
        ax2.axvline(-vertical, color='r')
        ax2.axhline(horizontal, color='r')
        ell = Ellipse(xy=(0,0), width=ell_width, height=ell_height, angle=0)
        ax2.add_artist(ell)
        ell.set_alpha(1)
        ell.set_fill(False)
        ell.set_edgecolor('red')
        ax2.plot(x_values, y_values)
        ax2.plot(x_values2_filtered, y_values2_filtered, color='green')
        ax2.plot(x_values3_filtered, y_values3_filtered, color='green')
        ax2.plot(-x_values2_filtered, y_values2_filtered, color='green')
        ax2.plot(-x_values3_filtered, y_values3_filtered, color='green')
        # ax2.set_title('%s' % title_data)
        ax1.set_xlabel(r'Curvature [$s^{3}$]')
        ax1.set_ylabel(r'Arc Strength')
        ax2.set_xlabel(r'Fringe Frequency [$10^{-3}$Hz]')
        ax2.set_ylabel(r'Delay [$\mu$s]')


        # hough_plot.subplots_adjust(top=0.75)
        hough_plot.colorbar(sec_im, ax=ax2, use_gridspec=True)
        hough_plot.tight_layout()
        hough_plot.subplots_adjust(top=0.87)
        # plt.show()
        hough_plot.savefig('./%s_hough.png' % archive, bbox_inches='tight',dpi=300)
        os.remove('dummy.gpi')
        os.remove('dummy2.gpi')
        os.remove('dummy.txt')
        os.remove('%s_hough.txt' % archive)

    return sec_final, sec_axes, a_value, a_error, arc_strength, power_q1, power_q2



def find_command(sec_axes, sec_shifted, archive, args):
    reso_0 = (sec_axes[1]-sec_axes[0])/( float(np.shape(sec_shifted)[0] -1))
    reso_1 = (sec_axes[3]-sec_axes[2])/( float(np.shape(sec_shifted)[1] -1))
    middle_0, middle_1 = np.shape(sec_shifted)
    middle_0 = middle_0 / 2.0
    middle_1 = middle_1 / 2.0
    x_range_0 = -25
    x_range_1 = 25
    y_range_0 = 0
    y_range_1 = 100
    new_name = archive + '_sec.txt'
    path = os.path.dirname(os.path.abspath(__file__))
    if args.options:
        part_string = args.options
    else:
        part_string = '--xrange=%i,%i --yrange=%i,%i --mask=0.5,3 --omask=25.0,15.0 --curves=0.2:7:800 --pdist=2,1 --q1' % (x_range_0, x_range_1, y_range_0, y_range_1)
    command_string = ('%s/parabfit --orig=%i,%i --res=%f,%f --out=%s_hough.txt \
--ssgpi=dummy.gpi --hggpi=dummy2.gpi --crop=dummy.txt %s %s'\
        %(path, middle_1, middle_0, reso_1, reso_0, archive, part_string, new_name))
    return command_string



def secondary_axes(secondary, mhzperbin, secperbin):
    conj_freq = np.fft.fftfreq(np.shape(secondary)[0], mhzperbin)
    max_conj_freq = np.max(conj_freq)
    min_conj_freq = np.min(conj_freq)
    # print mhzperbin
    # print secperbin
    # print max_conj_freq
    conj_time = np.fft.fftfreq(np.shape(secondary)[1], secperbin) * 1000
    max_conj_time = np.max(conj_time)
    min_conj_time = np.min(conj_time)
    # print max_conj_time  
    return min_conj_freq, max_conj_freq, min_conj_time, max_conj_time


def next_greater_power_of_2(x):  
    return 2**(x-1).bit_length()


def rebin(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    # print ''.join(evList)
    return eval(''.join(evList))


def rewrite_data(sec_spec, archive):
    # sec_spec = sec_spec / np.max(sec_spec)
    sec_rot = np.rot90(sec_spec, 1)
    shape = np.shape(sec_rot)
    # reso_0 = conj_time_range/( float(shape[0] -1))
    # reso_1 = conj_freq_range/( float(shape[1] -1))
    # x_max, y_max = np.unravel_index(sec_rot.argmax(), sec_rot.shape)
    # plt.imshow(sec_resampled, aspect='auto')
    # plt.show()
    new_fname = archive + '_sec.txt'
    open(new_fname, 'w').close()
    for index, x in np.ndenumerate(sec_rot):
        with open("%s" % new_fname, "a") as myfile:
            myfile.write("%s %s %.3f\n" % (index[0]+1, index[1]+1, x))


def new_curvature(sec_log, args, name, mjd, hour_length, ar_freq, diff):
        #get size and define full secondary spec
    sec_log_shape = np.shape(sec_log)
    sec_cut_plot = np.fliplr(sec_log)

    #define full array lengths
    x_axis_length = sec_log_shape[1]
    y_axis_length = sec_log_shape[0]
    med = np.median(sec_cut_plot)

    #choose start & stop values, get index of x-axis as difference
    start = 9.8/18.0*x_axis_length #left border
    stop = 14/18.0*x_axis_length #right border
    bottom = y_axis_length / 10 #exclue the bottom

    #mask the array
    sec_cut_mask = np.ma.array(sec_cut_plot, mask = True)
    sec_cut_mask.mask[bottom : ,start : stop ] = False

    #array of equally spaced y-values
    y_vals = np.linspace(0, 100, y_axis_length)

    #finds maximum x-value point on the y-axis
    x_maximums = (np.argmax(sec_cut_mask, axis=1)- x_axis_length / 2.0) * 100.0 / x_axis_length
    #gets rid of 10 bottom x-value points on y-axis
    mean_condition = np.asarray(np.nonzero(np.asarray(x_maximums) > 0)).squeeze()
    y_vals_cleaned = []
    x_max_cleaned = []
    for nonzero_max in mean_condition:
        y_vals_cleaned.append(y_vals[nonzero_max])
        x_max_cleaned.append(x_maximums[nonzero_max])


    opt = []
    opt, cov = curve_fit(new_func, y_vals_cleaned, x_max_cleaned)
    error = np.sqrt(np.diag(cov)[0])
    x_parabola_points = new_func(y_vals, opt[0])

    a_value = opt[0] #* 277.77

    curve_plot = plt.figure(7)
    plt.imshow(sec_cut_plot, extent=[-50,50,0,100], aspect='auto', origin='lower', vmin=med, vmax=med+10, cmap='binary')
    plt.autoscale(False)
    title_data = '\n %s MJD:%.2f  Duration:%.2fh\nFreqency:%.2fMHz  Bandwith:%.2fMHz  k=%.4f' % (name, mjd, hour_length, ar_freq, diff, a_value)        
    plt.title('Curve Fit %s' % title_data)
    plt.xlabel(r'fringe freq ($10^{-3}$Hz)')
    plt.ylabel(r'delay ($\mu$s)')
    plt.colorbar()
    plt.scatter(x_max_cleaned, y_vals_cleaned)
    plt.plot(x_parabola_points, y_vals, color='red')
    plt.plot(-x_parabola_points, y_vals, color='red')
    if not args.view:
        plt.show()
    if args.save:
        if args.gaussian_filter:
            curve_plot.savefig('Curve_Fit_freq=%.1f_mjd=%.2f_a=%s_filtered.png' % (ar_freq, mjd, a_value) , bbox_inches='tight')
        else:
            curve_plot.savefig('Curve_Fit_freq=%.1f_mjd=%.2f_a=%s_.png' % (ar_freq, mjd, a_value) , bbox_inches='tight')
    return a_value, error