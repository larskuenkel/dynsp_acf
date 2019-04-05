import numpy as np
from functions import *
import scipy.ndimage
import scipy.signal
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
# from scipy.optimize import least_squares


def help_variables(full_arr, minFreq, maxFreq, mjd_start, mjd_end):
    diff = maxFreq- minFreq
    length = (mjd_end - mjd_start) * 86400
    mjd = (mjd_end + mjd_start)/2.0
    mhzperbin = diff / (full_arr.shape[0] - 1)
    secperbin = length / (full_arr.shape[1] - 1)
    ar_freq = (maxFreq + minFreq) / 2
    return diff, length, mjd, mhzperbin, secperbin, ar_freq


def analyse_acf(full_arr, diff, length, args, funcs, fill_factor, minFreq):
    mhzperbin = diff / (full_arr.shape[0] - 1)
    secperbin = length / (full_arr.shape[1] - 1)

    rfi_frac = np.count_nonzero(full_arr) / float(np.size(full_arr))
    # plt.clf()
    # # full_arr_flat = np.ndarray.flatten(full_arr[full_arr != 0])
    # counts, edges = np.histogram(full_arr[full_arr != 0], bins=90)
    # # plt.hist(full_arr_flat, bins=60)
    # bar_width = edges[1] - edges[0]
    # plt.bar(edges[:-1], counts, bar_width)
    # bar_width = edges[1] - edges[0]
    # # plt.xlim(min(edges), max(edges))
    # print np.ma.mean(full_arr_masked)
    # print np.ma.median(full_arr_masked)
    # max_dist = (edges[np.argmax(counts)+1] + edges[np.argmax(counts)]) /2
    # print (edges[np.argmax(counts)+1] + edges[np.argmax(counts)]) /2
    # print edges
    # plt.show()

    #action = remove_mean(full_arr, args, minFreq, diff)
    action = np.ma.filled(full_arr, 0)

    # plt.imshow(action, aspect='auto')
    # plt.colorbar()
    # plt.show()

    # Correcting the sampling
    sampling_test = np.ones(np.shape(full_arr))
    np.putmask(sampling_test, full_arr==0, 0)

    sampling = scipy.signal.fftconvolve(sampling_test, np.flipud(np.fliplr(sampling_test)), mode = 'full')
    sampling[sampling == 0] = 1 #Avoids dividing by zero

    #Setting up array
    acf = scipy.signal.fftconvolve(action, np.flipud(np.fliplr(action)), mode = 'full')

    middle_f = acf.shape[0]/2
    middle_t = acf.shape[1]/2


    acf = np.divide(acf, sampling)

    acf[:,middle_t] = (acf[:,middle_t-1] + acf[:,middle_t+1]) / 2
    midACF_freq = acf[:,middle_t] / acf[middle_f,middle_t] #normalizes the ACF in frequency domain
    midACF_time = acf[middle_f,:] / acf[middle_f,middle_t]
    acf_norm = acf/np.amax(acf[int(middle_f):int(middle_f+middle_f*0.5),int(middle_t-middle_t*0.5):int(middle_t+middle_t*0.5)])

    # midACF_time = np.asarray(midACF_time) / midACF_time[len(midACF_time)/2]  #Normalizes array

    freqlag_ticks = np.linspace(-diff,diff,len(midACF_freq)) #Sets up tick labels
    timelag_ticks = np.linspace(-length/60/60,length/60/60,len(midACF_time)) #Sets up tick labels


    # Estimating the noise in the acf for the weighting function
    try:
        filtered = scipy.signal.wiener(acf_norm[int(middle_f):int(middle_f+50),int(middle_t-50):int(middle_t+50)])
        residual = filtered - acf_norm[int(middle_f):int(middle_f+50),int(middle_t-50):int(middle_t+50)]
    except:
        filtered = scipy.signal.wiener(acf_norm)
        residual = filtered - acf_norm
    # plt.imshow(residual, aspect='auto',vmin=-0.01, vmax=0.01)
    # plt.colorbar()
    # plt.show()
    noise = np.std(residual[int(0.25*residual.shape[0]):int(0.75*residual.shape[0]),int(0.25*residual.shape[1]):int(0.75*residual.shape[1])])
    # print noise
    start = 1
    for i in range(1, middle_f-2):
        start += 1
        if midACF_freq[middle_f+i+2] < 0.5:
                break

    end = start
    for i in range(start, middle_f-2):
        end += 1
        if midACF_freq[middle_f+i+2] < 0.0:
                break

    start2 = 1
    for i in range(1, middle_t-2):
        start2 += 1
        if midACF_time[middle_t+i+2] < 0.5:
                break

    end2 = start2
    for i in range(start2,middle_t-2):
        end2 += 1
        if midACF_time[middle_t+i+2] < 0.0:
                break

    try:
        best_end = find_fit2(funcs['func_f'], mhzperbin, midACF_freq, sampling[middle_f:,middle_t], noise, middle_f, end, start_point=start, index=0)
        extrapo_x, extrapo_y, midACF_freq, parameters_1, cov_1 = \
            fit_func(funcs['func_f'], mhzperbin, midACF_freq, sampling[middle_f:,middle_t], noise, middle_f, best_end, index=0)
        hwhm_freq = globals()[funcs['func_f'].__name__ + '_hwhm'](parameters_1)
    except:
        hwhm_freq = np.nan
        extrapo_x = 0
        extrapo_y = 0
    try:
        best_end2 = find_fit2(funcs['func_t'], secperbin/60/60, midACF_time, sampling[middle_f,middle_t:], noise, middle_t, end2, start_point=start2, index=1)
        extrapo_x2, extrapo_y2, midACF_time, parameters_2, cov_2 = \
            fit_func(funcs['func_t'], secperbin/60/60, midACF_time, sampling[middle_f,middle_t:], noise, middle_t, best_end2, index=1)
        e_time = globals()[funcs['func_t'].__name__ + '_e'](parameters_2)
    except:
        e_time = np.nan
        extrapo_x2 = 0
        extrapo_y2 = 0
    try:
        print 'Fitting regions:'
        print ('%s - %s: %s, %s - %s: %s' % (start, end, best_end ,start2 , end2, best_end2))
    except:
        print 'Error during fit.'
    if np.isnan(hwhm_freq + e_time):
        freq_error, time_error = np.nan, np.nan
    else:
        stat_error = (fill_factor * ( diff * length * rfi_frac / (hwhm_freq * (e_time*60*60))))**(-0.5)
        # print stat_error
        error_1 = np.sqrt(np.diag(cov_1))
        error_2 = np.sqrt(np.diag(cov_2))
        freq_error = globals()[funcs['func_f'].__name__ + '_h_err'](parameters_1, error_1)
        freq_error = np.sqrt(freq_error ** 2 + (hwhm_freq * stat_error) ** 2)
        time_error = globals()[funcs['func_t'].__name__ + '_e_err'](parameters_2, error_1)
        time_error = np.sqrt(time_error ** 2 + (e_time * stat_error) ** 2)
    return acf_norm, length, middle_f, middle_t, freqlag_ticks, midACF_freq, extrapo_x, extrapo_y, timelag_ticks, midACF_time, extrapo_x2, extrapo_y2, \
        hwhm_freq, freq_error, e_time, time_error


def fit_func(function, xperbin, value_array, sampling, noise, mid, end_point, index):
    # Fits a function to the central part of the acf. 
    fitx = np.linspace(xperbin,end_point*xperbin,end_point)  # X-values starting at index of 1
    sampling = sampling[0:fitx.shape[0]]
    vals = np.copy(value_array)[int(mid+1):int(mid+end_point+1)]  # Array of values we want to fit to
    err_acf = np.divide(np.max(sampling),sampling)
    err_acf = (err_acf**0.5) * noise
    extrapol_x = np.linspace(0,end_point*xperbin,end_point + 1) #Array of x-values including zero point
    # plt.errorbar(fitx, vals, yerr=err_acf)
    # plt.show()
    opt, cov = curve_fit(function,fitx,vals,sigma=err_acf)  # Fits a curve and gives values for a b and c as indices of opt
    extrapol_y = function(extrapol_x,*opt)  # Extrapolated values
    squares = np.sum((function(fitx, *opt) - vals)**2)
    return extrapol_x, extrapol_y, value_array, opt, cov


def get_hwhm_new(function, *parameters):
    # Finds half maximum for normalized function
    find_value = lambda x_value: function(x_value, *parameters) - 0.5
    x_hwhm, status = leastsq(find_value, [1.0])
    return x_hwhm


def find_fit2(function, xperbin, value_array, sampling, noise, mid, end_point, index, start_point=10):
    parameter = []
    errs = []
    if index == 0:
        function_err = function.__name__ + '_h_err'
        function_hwhm = function.__name__ + '_hwhm'
    else:
        function_err = function.__name__ + '_e_err'
        function_hwhm = function.__name__ + '_e'  
    start = start_point
    for end in range(start, end_point + 1):
        _, _, _, para, cov = fit_func(function, xperbin, value_array, sampling, noise, mid, end, index)
        try:
            err = globals()[function_err](para, cov)
        except ValueError:
            err = np.nan
            pass
        errs.append((end,float(err)))
    errs = np.asarray(errs)
    x_values = np.arange(start, end_point) * xperbin
    # plt.clf()
    # graph_error = plt.figure(1)
    # plt.plot(x_values, errs[:,1])
    # plt.show()
    # print errs
    fit_end = errs[np.nanargmin(errs[:,1]),0]
    return fit_end


def remove_mean(full_arr, args, minFreq, diff):
    #Remove an interpolated channel-dependent mean.
    try:
        full_arr_masked = np.ma.masked_equal(full_arr, 0)
        chans = np.linspace(minFreq, minFreq + diff, num=full_arr.shape[0])
        means = []
        for j in range(full_arr.shape[0]):
            means.append(np.ma.mean(full_arr_masked[j,:]))

        mean_condition = np.asarray(np.nonzero(np.asarray(means) > 0)).squeeze()
        chans_cleaned = []
        means_cleaned = []

        for nonzero_mean in mean_condition:
            chans_cleaned.append(chans[nonzero_mean])
            means_cleaned.append(means[nonzero_mean])

        def err(x, chans, means):
            return x[0] * np.power(chans, x[1]) - means
        def power_law(x, chans):
            return x[0] * np.power(chans, x[1])

        params, err_mean = leastsq(err, [10.,-3.], args=(chans_cleaned, means_cleaned))#, loss='soft_l1')

        mean_extra2 = []
        for channel in chans:
            mean_extra2.append(power_law(params, channel))
        # plt.plot(chans, means, chans, mean_extra2)
        # plt.show()
        if args.normal_f:
            for j in range(full_arr.shape[0]):
                full_arr[j,:] = full_arr[j,:] / mean_extra2[j]
                full_arr_masked = np.ma.masked_equal(full_arr, 0)
            action = full_arr_masked - np.ma.mean(full_arr_masked)
        else:
            full_arr_masked = np.ma.masked_equal(full_arr, 0)
            action = full_arr_masked
            for j in range(full_arr.shape[0]):
                action[j,:] = action[j,:] - mean_extra2[j]
            action = action - np.ma.mean(action)
    except:
        print 'Fitting the frequency dependent mean failed.'
        full_arr_masked = np.ma.masked_equal(full_arr, 0)
        action = full_arr_masked - np.ma.mean(full_arr_masked)

    return full_arr, action
