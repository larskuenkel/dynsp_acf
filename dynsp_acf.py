#!/bin/usr/python

import psrchive as psr
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import subprocess
import argparse
import os
import scipy.ndimage


def func(x, a, b, c):
    # Exponential function
    # 3 Parameters
    return a * np.exp(-b * x) + c

def func2(x, b):
    # Normalised Gaussian function
    # 1 Parameter
    return np.exp(-b * x**2)

def func3(x, a, b):
    # Exponential function
    # 2 Parameters
    return a * np.exp(-b * x)


def fit_func(function, xperbin, value_array, mid, end_point):
    # Fits a function to the central part of the acf. 
    # Replaces center value with extrapolated value.
    num_para = function.func_code.co_argcount - 1  # Number of parameters of the fitting function
    fitx = np.linspace(xperbin,end_point*xperbin,end_point)  # X-values starting at index of 1
    vals = np.copy(value_array)[mid+1:mid+end_point+1]  # Array of values we want to fit to
    extrapol_x = np.linspace(0,end_point*xperbin,end_point + 1) #Array of x-values including zero point
    opt, cov = curve_fit(function,fitx,vals,p0=num_para*[0])  # Fits a curve and gives values for a b and c as indices of opt
    extrapol_y = function(extrapol_x,*opt)  # Extrapolated values
    adj_array = np.copy(value_array)  # The following lines replace the 0-value and normalize
    adj_array[mid] = extrapol_y[0]
    adj_array_norm = adj_array / max(adj_array)
    extrapol_y_norm = extrapol_y / max(adj_array)
    initial_peak = value_array[mid] /max(adj_array)  # Shows where the initial 0-value was
    return extrapol_x, extrapol_y_norm, adj_array_norm, initial_peak, opt, cov

def get_hwhm_new(function, *parameters):
    # Finds half maximum for normalized function
    find_value = lambda x_value: function(x_value, *parameters) - 0.5
    x_hwhm, status = leastsq(find_value, [1.0])
    return x_hwhm


parser = argparse.ArgumentParser(description='Select the archive.')
parser.add_argument('ar', help='The chosen archive')
parser.add_argument('-n','--normal', action='store_true', help='Normalize the dynamic spectrum.')
parser.add_argument('-s', '--save', action='store_true', help='Save the images.')
parser.add_argument('-v', '--view', action='store_true', help='Do not view the plots.')
parser.add_argument('-t', '--tscale', type=float, default=1, help='Scales down the plotted time scale in acf plots.')
parser.add_argument('-f', '--fscale', type=float, default=1, help='Scales down the plotted frequency scale in acf plots.')
args = parser.parse_args()
ar_path = str(os.path.abspath(args.ar)) 

arch = psr.Archive_load(args.ar) #load a certain archive
minFreq = arch.get_centre_frequency() -  (arch.get_bandwidth()/2)
maxFreq = arch.get_centre_frequency() +  (arch.get_bandwidth()/2)
length = arch.integration_length()
diff = maxFreq-minFreq
splitname = args.ar.split('.')
t_scale = args.tscale
f_scale = args.fscale
# pulsar = splitname[0]
# arfreq = '%s.%s' % (splitname[1], splitname[2])
# mjd = '%s.%s' % (splitname[3], splitname[4])
try:
    full_arr = np.genfromtxt('./%s.txt' % args.ar)
    print 'Existing %s.txt is used' % args.ar 
except IOError:
    print 'New %s.txt is being created' % ar_path
    subprocess.call('/usr/share/psrchive-hacked/bin/dynamic_spectra -o %s.txt -f matlab %s' % (args.ar, ar_path), shell=True)
    print 'Done'
    full_arr = np.genfromtxt('./%s.txt' % args.ar)
full_arr = np.rot90(full_arr, 3)
full_arr = np.fliplr(full_arr)

mhzperbin = diff / (full_arr.shape[0] - 1)
secperbin = length / (full_arr.shape[1] - 1)

full_arr_med = scipy.ndimage.filters.median_filter(full_arr, size=(5,5), mode='nearest')


np.putmask(full_arr, full_arr==0, full_arr_med)

# This normalizes each sub-int
if args.normal:
    for i in range(full_arr.shape[1]):
        maxi = np.max(full_arr[:,i])
        for j in range(full_arr.shape[0]):
            full_arr[j,i] = full_arr[j,i]/maxi
            if np.isnan(full_arr[j,i]):
                full_arr[j,i] = 0



#Creating ACF 

    
#Get acf for above along with plotting ACF time and frequency for 6 hour obs and hwhm
action = full_arr
action = action - np.mean(action)
padded = np.zeros((full_arr.shape[0]*2 - 1,full_arr.shape[1]*2 - 1)) #Creates an array of zeros of 2x dimensions
padded[:full_arr.shape[0],:full_arr.shape[1]] = action[:,:] #Fills padded with previous array
middle = padded.shape[0]/2
middle2 = padded.shape[1]/2

#Setting up array
padf = np.fft.fft2(padded) #2d fft
padf_sq = (np.abs(padf))**2
acfpad = np.fft.ifft2(padf_sq) #inverse fft
acfpad = np.real(np.fft.fftshift(acfpad)) #real part of shifted
acfpad_norm = acfpad/np.amax(acfpad)

midACF_freq = acfpad[:,middle2] / acfpad[middle,middle2] #normalizes the ACF in frequency domain
midACF_time = acfpad[middle,:]


midACF_time[middle2] = (midACF_time[middle2-1] + midACF_time[middle2+1]) / 2 
#Replaces zero lag freq spike with average around it
midACF_time = np.asarray(midACF_time) / midACF_time[len(midACF_time)/2]  #Normalizes array

freqlag_ticks = np.linspace(-diff,diff,len(midACF_freq)) #Sets up tick labels
timelag_ticks = np.linspace(-length/60/60,length/60/60,len(midACF_time)) #Sets up tick labels



#This next part is to fit the exponential to the points after the zero lag spike and extrapolate the zero-point #from that


# for i in range(middle-2):
#     #This finds the point in the ACF freq array where the function starts to increase
#     if midACF_freq[middle+i+2] > midACF_freq[middle+i]:
#             end = i
#             # print end
#             break

for i in range(middle-2):
    if midACF_freq[middle+i+2] < 0.05:
            end = i
            break

for i in range(middle2-2):
    if midACF_time[middle2+i+2] < 0.4:
            end2 = i
            break

extrapo_x, extrapo_y_norm, adj_midACF_freq_norm, initial_freq_peak, parameters_1, errors_1 = \
    fit_func(func3, mhzperbin, midACF_freq, middle, end)
extrapo_x2, extrapo_y2_norm, adj_midACF_time_norm, initial_time_peak, parameters_2, errors_2 = \
    fit_func(func2, secperbin/60/60, midACF_time, middle2, end2)
# hwhm if it can't be diectly derived from the fitting parameters
# hwhm_freq = get_hwhm_new(func3, *parameters_1)
hwhm_freq = np.log(2)/float(parameters_1[1])
cov_diag = np.sqrt(np.diag(errors_1))
freq_error = np.log(2) * parameters_1[1]**(-2) * cov_diag[0]
# Error if a * np.exp(-b * x) + c is used
# freq_error = np.sqrt((cov_diag[0]*parameters_1[0]*np.exp(-parameters_1[1]*hwhm_freq))**2+(cov_diag[1]*parameters_1[0]*\
#     hwhm_freq*np.exp(-parameters_1[1]*hwhm_freq))**2+(cov_diag[2])**2) 
hwhm_time = np.sqrt(1/float(parameters_2))
cov_diag2 = np.sqrt(np.diag(errors_2))
time_error = 0.5*float(parameters_2[0]**(-1.5))*cov_diag2[0]
print ("Calculated Values: \n hwhm_freq = %s +- %sMHz, 1/e_time = %s +- %s min" %
 (hwhm_freq, freq_error, hwhm_time*60, time_error*60))


#Plotting functions in matplotlib, pretty self explanatory
dyn_plot = plt.figure(1)
plt.imshow(full_arr,origin='lower',aspect='auto',extent=[0,length/60/60,minFreq,maxFreq],cmap='jet')
plt.title('Dyamic Spectrum %s' % args.ar)
plt.xlabel('Time (hours)')
plt.xticks(np.round(np.linspace(0,round(length/60/60,1),7,endpoint=True),1))
plt.ylabel('Freq (MHz)')
plt.colorbar()
plt.yticks(np.round(np.linspace(minFreq,maxFreq,9),1))
# plt.savefig('%s_dyn.png' % args.ar, bbox_inches='tight')

acf_2dplot = plt.figure(2)
plt.imshow(acfpad_norm, origin='lower', aspect='auto', extent=[-length/60,length/60,-diff,+diff], interpolation='None')
plt.title('Autocorrelation Function %s' % args.ar)
plt.xlim(-length/60/2*t_scale, length/60/2*t_scale)
plt.ylim(0,diff*f_scale*0.5)
plt.xlabel('Time Lag (min)')
plt.ylabel('Frequency Lag (MHz)')
plt.colorbar()

slice_plot = plt.figure(3)
ax_big = slice_plot.add_subplot(111,frameon=False)
ax_big.set_title('%s Autocorrelation Functions' % args.ar)
ax_big.set_yticklabels('')
ax_big.set_xticklabels('')
ax_big.title.set_y(1.05)

ax1 = slice_plot.add_subplot(121)
ax1.plot(freqlag_ticks,adj_midACF_freq_norm,extrapo_x,extrapo_y_norm)
ax1.plot(0,initial_freq_peak,'g^')
ax1.set_xlim(-diff*f_scale*0.5,diff*f_scale*0.5)
ax1.set_xlabel('Frequency Lag (MHz)')
ax1.set_ylabel('Autocorrelation')
ax1.set_title('Autocorr Freq: hwhm  = %s +- %s MHz' % (str(round(hwhm_freq,2)), str(round(freq_error,2))))
ax1.title.set_fontsize(11)

ax2 = slice_plot.add_subplot(122)
ax2.plot(timelag_ticks*60,adj_midACF_time_norm,extrapo_x2*60,extrapo_y2_norm)
# ax2.plot(0,initial_time_peak,'g^')
ax2.set_xlabel('Time Lag (min)')
ax2.set_title('Autocorr Time: 1/e  = %s +- %s min' % (str(round(hwhm_time*60,2)), str(round(time_error*60,2))))
ax2.set_ylim(-.2,1)
ax2.set_xlim(-length/60/2*t_scale,length/60/2*t_scale)
ax2.title.set_fontsize(11)

if args.save:
    dyn_plot.savefig('%s_dyn.png' % args.ar, bbox_inches='tight')
    acf_2dplot.savefig('%s_2dacf.png' % args.ar, bbox_inches='tight')
    slice_plot.savefig('%s_acf_slice.png' % args.ar, bbox_inches='tight')
if not args.view:
    plt.show()


# with open("/scratch/lkuenkel/J0814/parameters_2.txt", "a") as myfile:
#     myfile.write("\n %s %s %s %s"
#         % (mjd, arfreq, hwhm_freq, hwhm_time*60))

