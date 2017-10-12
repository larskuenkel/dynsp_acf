import numpy as np
import create_dyn


def stitch_dyn(ar, args):
    dummy = []
    full_name = ''
    for archive in ar:
        dummy.append(create_dyn.get_dyn(archive, args, 0))
    dummy = np.asarray(dummy)
    name = dummy[0,1]
    site = dummy[0,2]
    minFreq = np.min(dummy[:,3])
    maxFreq = np.max(dummy[:,4])
    mjd_start = np.min(dummy[:,5])
    mjd_end = np.max(dummy[:,6])
    ra = dummy[0,7]
    dec = dummy[0,8]
    total_length = (mjd_end - mjd_start) * 86400
    total_bandwith = maxFreq - minFreq
    secperbin = (dummy[0,6] - dummy[0,5]) * 86400 / (dummy[0,0].shape[1] - 1)
    mhzperbin = (dummy[0,4] - dummy[0,3]) / (dummy[0,0].shape[0] - 1)
    totaltimebins = np.ceil(total_length / secperbin) + 1
    totalfreqbins = np.ceil(total_bandwith / mhzperbin) + 1
    full_arr = np.zeros((int(totalfreqbins), int(totaltimebins)))
    dummy[:,5] = np.subtract(dummy[:,5], dummy[0,5]) * 86400 / secperbin
    dummy[:,5] = np.around(dummy[:,5].astype(np.double))
    dummy[:,3] = np.subtract(dummy[:,3], dummy[0,3]) / mhzperbin
    dummy[:,3] = np.around(dummy[:,3].astype(np.double))
    edge = ''
    if args.edges:
   		edge = 'edge'
    full_name = '%s_%.2f_%.2f_%s_dyn.txt' % (name, (mjd_end + mjd_start)/2, (maxFreq + minFreq)/2, edge)
    # print dummy
    for dyn in dummy:
        full_arr[int(dyn[3]):int(dyn[0].shape[0]+dyn[3]),int(dyn[5]):int(dyn[0].shape[1]+dyn[5])] = dyn[0]
    print 'Creating combined dynamic spectrum %s' % full_name
    create_dyn.write_dyn(full_arr, full_name, name, site, minFreq, maxFreq, mjd_start, mjd_end, ra, dec)
    return ('%s' % full_name)
