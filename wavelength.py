import numpy as np

from scipy import interpolate

import matplotlib.pyplot as plt

import create_dyn

'''


#This function takes the dynamic spectrum with frequency domain and switches it to a wavelength domain
 
def frequency_to_wavelength(full_arr, minFreq, maxFreq):
 
	frequencies = np.linspace(minFreq,maxFreq,np.shape(full_arr)[0])
	wavelengths = 3.0*10**8 / frequencies
	minWav = np.min(wavelengths)
	maxWav = np.max(wavelengths)
	evens = np.linspace(np.floor(wavelengths[0]),np.floor(wavelengths[-1]),2^8) #Evenly spaced wavelengths
	full_arr_flat = np.ndarray.flatten(full_arr)
	print np.shape(full_arr_flat)
	times = np.linspace(0, 100, np.shape(full_arr)[1])
	print np.shape(times)
	print np.shape(full_arr)
	coord = []

	print(wavelengths)
	
	for wave in wavelengths:
		for time in times:
			coord.append((wave,time))

	print np.shape(coord)

	x = np.shape(times)[0]*1j
	y = np.shape(frequencies)[0]*3j
	
	grid_x, grid_y = np.mgrid[minWav:maxWav:y, 0:100:x]

	print grid_x


	full_arr = interpolate.griddata(coord, full_arr_flat, (grid_x, grid_y), method='linear')

	full_arr = np.flipud(full_arr)

	return full_arr

'''
"""

 #Credit for this code goes to Bill Coles

def frequency_to_wavelength(full_arr, minFreq, maxFreq):

  #Function to resample an array from equal spacing in frequency to equal spacing in wavelength.

    full_arr_in = full_arr


    frequencies = np.linspace(minFreq,maxFreq,np.shape(full_arr_in)[0])

    #print(frequencies)

    nf,nt = full_arr_in.shape

    #changes the size of the grid

    nf = nf
    nt = nt

    full_arr = np.zeros((nf,0))

    wavelengths = 300000000 /frequencies

    #print (wavelengths)

    lambdaeq = np.arange(wavelengths[0],wavelengths[-1],-(wavelengths[0]-wavelengths[-1])/nf)
    for it in range(nt):
        #full_arr[:,it] = interpolate.interp1d(wavelengths[::-1],full_arr_in[::-1,it])(lambdaeq[:nf])
        full_arr[:,it] = interpolate.interp1d(wavelengths[::-1],full_arr_in[::-1,it])(lambdaeq[:nf])



    return full_arr

"""

#Credit for this code goes to Bill Coles

def frequency_to_wavelength(full_arr, minFreq, maxFreq):

  #Function to resample an array from equal spacing in frequency to equal spacing in wavelength.

    full_arr_in = full_arr


    frequencies = np.linspace(minFreq,maxFreq,np.shape(full_arr_in)[0])

    #print(frequencies)

    nf,nt = full_arr_in.shape

    #changes the size of the grid

 

    full_arr = np.zeros((nf,nt))

    wavelengths = 300000000 /frequencies

    #print (wavelengths)

    lambdaeq = np.arange(wavelengths[0],wavelengths[-1],-1*(wavelengths[0]-wavelengths[-1])/nf)
    for it in range(nt):
        #full_arr[:,it] = interpolate.interp1d(wavelengths[::-1],full_arr_in[::-1,it])(lambdaeq[:nf])
        full_arr[:,it] = interpolate.interp1d(wavelengths[::-1],full_arr_in[::-1,it])(lambdaeq[:nf])


    return full_arr















'''









 #Credit for this code goes to Jason Rosenblum


def find_nearest(my_array, target):
    
	#returns nearest above and nearest below
    
	diff = my_array - target #Relates all values in array to the target value
	mask_neg = np.ma.less_equal(diff, 0) #Masks values below 0
	mask_pos = np.ma.greater_equal(diff, 0) #Masks values above zero
    
	# We need to mask the negative differences and zero
	# since we are looking for values above
    
	masked_diff_pos = np.ma.masked_array(diff, mask_neg)
	masked_diff_neg = np.ma.masked_array(diff, mask_pos)
 
	#Returns the greatest negative value and the smallest positive value (ie the values on either side)
	return masked_diff_pos.argmin(),masked_diff_neg.argmax() 
 
def linearize(x,x1,x2,y1,y2):

	#Returns the y value that corresponds to the nearest slope
	#Arguments:
	#	x1,y1,x2,y2 are the two points
	#	x is the array of values
	return (((y2-y1)/(x2-x1))*(x-x1)+y1)
    
def interpol(x1,x2,y1,y2):

	#Interpolates the values in between the wavelengths in order to get the point that matches up
	#I can explain better later

	return linearize(np.linspace(x1,x2,101),x1,x2,y1,y2)
 

#Below maps the dynamic spectrum onto the wavelength vs time axes. The way it works is to take all of the frequencies, turn them into wavelengths then turn the wavelengths array into an evenly spaced array.
 #then other stuff

 
 
frequencies = np.linspace(minFreq,maxFreq,488)
wavelengths = 3.0*10**8 / frequencies
evens = np.linspace(np.floor(wavelengths[0]),np.floor(wavelengths[-1]),488) #Evenly spaced wavelengths
pos = [] #Array of values greater than target
neg = [] #Array of values less than target
 
wave_adj = np.zeros(full_arr.shape)
 
for i in range(len(evens)):
    
	#Find closest values to evens[i], above and below, in the original wavelength array
	pos,neg = find_nearest(wavelengths,evens[i]) 
	
	#Find the equivalent percentage between pos and neg that evens[i] is 
	percent = int((evens[i]-wavelengths[pos])/np.abs(wavelengths[pos]-wavelengths[neg]))
    	
	
	for j in range(full_arr.shape[1]):
		#This finds the equivalent y-value by interpolating a straight line between pos and neg
		#then finding what the new value for intensity should be to put into the new ds (wave_adj)
        		value = interpol(wavelengths[pos],wavelengths[neg],full_arr[pos,j],full_arr[neg,j])[percent]
 
    		wave_adj[i,j] = value
 
#This plots the new dynamic spectrum so you can see the change
plt.clf()
plt.imshow(wave_adj,aspect='auto',origin='lower')


'''