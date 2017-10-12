import numpy as np
import scipy as sp
# from scipy import ndimage, signal, fft
import matplotlib.pyplot as plt



def gaussian_window_convolve(sec_spec, args):

	shape1, shape2 = np.shape(sec_spec)

	input1, input2 = args.gaussian_filter

	input1 = int(input1 * shape1 *.005)

	input2 = int(input2 * shape2 *.01)

	#call shape of array and use it divided by 10 in place of numbers 5 and 10. 

	filtered = sp.ndimage.filters.gaussian_filter(sec_spec, {input1, input2})


	#filtered = signal.convolve(sec_spec, filtered, mode= 'same') / sum(filtered)



	return filtered