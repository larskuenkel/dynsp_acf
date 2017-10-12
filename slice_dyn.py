import numpy as np 



def slice_dyn(full_arr, args,  minFreq, maxFreq, mjd_start, mjd_end):
    freq_div = args.slice[0]
    time_div = args.slice[1]
    new_arr = []
    new_minFreq = []
    new_maxFreq = []
    new_mjd_start = []
    new_mjd_end = []
    ini_shape = np.shape(full_arr)
    new_shape = np.divide(ini_shape, args.slice)
    freqs = np.linspace(minFreq, maxFreq, ini_shape[0])
    times = np.linspace(mjd_start, mjd_end, ini_shape[1])

    for freq_slice in range(0, freq_div):
        for time_slice in range(0, time_div):
            try:
                new_arr.append(full_arr[new_shape[0]*freq_slice : new_shape[0]*(freq_slice+1)-1, new_shape[1]*time_slice : new_shape[1]*(time_slice+1)-1])
                new_minFreq.append(freqs[new_shape[0]*freq_slice])
                new_maxFreq.append(freqs[new_shape[0]*(freq_slice+1)-1])
                new_mjd_start.append(times[new_shape[1]*time_slice])
                new_mjd_end.append(times[new_shape[1]*(time_slice+1)-1])
            except:
                pass

    return new_arr, new_minFreq, new_maxFreq, new_mjd_start, new_mjd_end

