import numpy as np
import matplotlib.pyplot as plt
import functions
from matplotlib.patches import Ellipse
# import secondary



def plot_dyn(full_arr, args, length, diff, minFreq, maxFreq, ar_freq, name, mjd, ar):
    hour_length = length / 3600.0
    title_data = '\n %s MJD:%.2f  Duration:%.2fh\nFrequency:%.2fMHz  Bandwith:%.2fMHz' % (name, mjd, hour_length, ar_freq, diff)
    plt.close('all')

    dyn_plot = plt.figure(1)
    plt.imshow(full_arr,origin='lower',aspect='auto',extent=[0,length/60/60,minFreq,maxFreq],cmap='jet', interpolation='None')
    plt.title('%s' % title_data)
    # plt.title('Dyamic Spectrum J0814+7429 MJD:57166')
    plt.xlabel('Time (hours)')
    plt.ylabel('Frequency (MHz)')
    plt.colorbar(use_gridspec=True)
    dyn_plot.tight_layout()
    if args.save:
        dyn_plot.savefig('%s_dyn.png' % ar, bbox_inches='tight', dpi=300)



def plot_acf(full_arr, args, acf_norm, length, minFreq, maxFreq, middle, middle2, 
        diff, freqlag_ticks, midACF_freq, extrapo_x, extrapo_y, timelag_ticks, midACF_time, extrapo_x2, extrapo_y2, ar, 
        hwhm_freq, freq_error, e_time, time_error, mhzperbin, secperbin, ar_freq, name, mjd, site):
    t_range = args.trange
    f_range = args.frange
    if t_range != 0:
        t_scale = length*60 / (t_range*4)
    else:
        t_scale = 1
    if f_range != 0:
        f_scale = diff / (f_range*4)
    else: f_scale = 1
    hour_length = length / 3600.0
    title_data = '\n %s MJD:%.2f  Duration:%.2fh\nFreqency:%.2fMHz  Bandwith:%.2fMHz' % (name, mjd, hour_length, ar_freq, diff)


    min_y_value = np.min(acf_norm[int(middle):int(middle+middle*0.5*f_scale),int(middle2-middle2*0.5*t_scale):int(middle2+middle2*0.5*t_scale)])

    acf_2dplot = plt.figure(2)
    plt.imshow(acf_norm, origin='lower', aspect='auto', extent=[-length/60,length/60,-diff,+diff], interpolation='None', vmin=-0.1, vmax=1)
    plt.title('%s' % title_data)
    # plt.title('Autocorrelation Function J0814+7429 MJD:57166')
    if t_range != 0:
        plt.xlim(-t_range, t_range)
    else:
        plt.xlim(-length/120, length/120)
    plt.ylim(0, diff/2)
    if f_range != 0:
        plt.ylim(ymax=f_range)
    plt.xlabel('Time Lag (min)')
    plt.ylabel('Frequency Lag (MHz)')
    plt.colorbar(use_gridspec=True)
    acf_2dplot.tight_layout()

    slice_plot = plt.figure(3)
    ax_big = slice_plot.add_subplot(111,frameon=False)
    ax_big.set_title('%s' % title_data)
    ax_big.set_yticklabels('')
    ax_big.set_xticklabels('')

    ax1 = slice_plot.add_subplot(121)
    try:
        ax1.plot(freqlag_ticks,midACF_freq,extrapo_x,extrapo_y)
    except NameError:
        ax1.plot(freqlag_ticks,midACF_freq)
    if f_range != 0:
        ax1.set_xlim(-f_range,f_range)
    else:
        ax1.set_xlim(-diff/2,diff/2)
    ax1.set_ylim(-.2,1.1)
    ax1.set_xlabel('Frequency Lag (MHz)')
    ax1.set_ylabel('Autocorrelation')
    ax1.set_title(r'$\nu_{d}$  = %.2f $\pm$ %.2f MHz' % (hwhm_freq, freq_error))
    ax1.title.set_fontsize(11)

    ax2 = slice_plot.add_subplot(122)
    try:
        ax2.plot(timelag_ticks*60,midACF_time,extrapo_x2*60,extrapo_y2)
    except NameError:
        ax2.plot(timelag_ticks*60,midACF_time)
    ax2.set_xlabel('Time Lag (min)')
    ax2.set_title(r'$\tau_{d}$  = %.2f $\pm$ %.2f min' % (e_time*60, time_error*60))
    ax2.set_ylim(-.2,1.1)
    if f_range != 0:
        ax2.set_xlim(-t_range,t_range)
    else:
        ax2.set_xlim(-length/120,length/120)
    ax2.title.set_fontsize(11)
    slice_plot.tight_layout()
    if args.save:
        acf_2dplot.savefig('%s_2dacf.png' % ar, bbox_inches='tight')
        slice_plot.savefig('%s_acf_slice.png' % ar, bbox_inches='tight')


    
def plot_secondary(sec_spec, args, mhzperbin, secperbin, ar, ar_freq, name, diff, length, mjd, site, sec_axes, curv):
    hour_length = length / 3600.0
    title_data = '\n %s MJD:%.2f  Duration:%.2fh\nFrequency:%.2fMHz  Bandwith:%.2fMHz' % (name, mjd, hour_length, ar_freq, diff)

    sec_plot = plt.figure(4)
    ax = sec_plot.add_subplot(111, aspect='auto')
    sec_mean = np.mean(sec_spec)
    # sec_2 = np.copy(sec_spec)
    # sec_2[0:40,:] = -100
    # mid_shape = np.shape(sec_2)[1]
    # sec_2[:,mid_shape -40: mid_shape+40] = -100
    # max_cbar = np.amax(sec_2)
    # min_conj_freq, max_conj_freq, _, max_conj_time = secondary.secondary_axes(sec_spec, mhzperbin, secperbin)
    # max_conj_time = 1000* max_conj_time
    sec_im = ax.imshow(sec_spec, cmap='binary', aspect='auto', extent=[sec_axes[2], sec_axes[3], 0, sec_axes[1]], vmax=sec_mean+8, vmin=sec_mean+2, origin='lower')
    # ax.set_xlim(-5,5)
    # ax.set_ylim(0,40)
    ax.autoscale(False)
    # if args.hough:
    #     x_values = np.linspace(sec_axes[2], sec_axes[3], np.shape(sec_spec)[1])
    #     y_values = functions.parabola(x_values, curv)
    #     file = open('dummy.gpi', 'r')
    #     file2 = open('dummy2.gpi', 'r')
    #     dummy_lines = file.readlines()
    #     dummy2_lines = file2.readlines()
    #     x_dist = float(dummy2_lines[4].split()[-4].strip('('))
    #     y_dist = float(dummy2_lines[4].split()[-2])
    #     x_values2 = (x_values + x_dist)
    #     y_values2 = (y_values - y_dist)
    #     x_values3 = (x_values - x_dist)
    #     y_values3 = (y_values + y_dist)
    #     vertical = float(dummy_lines[27].split()[0])
    #     horizontal = float(dummy_lines[36].split()[1])
    #     ell_height = 2*float(dummy_lines[20].split('*')[0])
    #     ell_width = 2*float(dummy_lines[20].split('/')[1].split(')')[0])
    #     ax.axvline(vertical, color='r')
    #     ax.axvline(-vertical, color='r')
    #     ax.axhline(horizontal, color='r')
    #     ell = Ellipse(xy=(0,0), width=ell_width, height=ell_height, angle=0)
    #     ax.add_artist(ell)
    #     ell.set_alpha(1)
    #     ell.set_fill(False)
    #     ell.set_edgecolor('red')
    #     ax.plot(x_values, y_values)
    #     ax.plot(x_values2, y_values2, color='green')
    #     ax.plot(x_values3, y_values3, color='green')
    ax.set_title('%s' % title_data)
    ax.set_xlabel(r'Fringe Frequency ($10^{-3}$Hz)')
    ax.set_ylabel(r'Delay ($\mu$s)')
    sec_plot.colorbar(sec_im, use_gridspec=True)
    sec_plot.tight_layout()
    if args.save:
        sec_plot.savefig('%s_sec.png' % ar, bbox_inches='tight', dpi=300)

def show_plots(args):

    if not args.view:
        plt.show()
    # if not args.acf:
    #     dyn_plot.clear()
    #     acf_2dplot.clear()
    # if sec_spec != []:
    #     sec_plot.clear()