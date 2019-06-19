#Alan Long 6/16/16
#Last edited: Alan Long 5/17/18

#This code takes data and filters it against a noise file through
#a Fourier transform, a method called Weiner filtering, then returning
#the filtered signal. It accepts two arrays of the same length noise and unfiltered
# and two ints exponent and fixedfreq, and returns an array of the same length

#This code is based on Aya's Matlab code

#Note of caution: There is a known issue with the first and last few points of the filtered results. 
#If you plot them over the original data, you may notice some curled up or curled down edges
#Usually, we include some extra points out side our region of interest 
#and just trim of the first and last 200 points after filtering. 

#import necessary modules
import numpy as np
import numpy.fft as fft                         
import math
import scipy.io as io


def weinerfilter(noise, unfiltered,fixedfreq,exponent): #fixedfreq is a value where the amplitude is taken. Recomended 3000. exponent is recomended 2

  unfiltered=np.array(unfiltered)
  noise=np.array(noise)
  

    #We take the fourier transform of both the signal and the noise
    unfil_fourier=fft.fft(unfiltered)
    noise_fourier=fft.fft(noise)

    #We now make a theoretical signal based on the given frequency with exponent scaling
    ampl=np.absolute(unfil_fourier[fixedfreq]**2)
    fake_sig=ampl*((fixedfreq+1)**exponent)*(np.arange(1.0,len(unfil_fourier)+1)**(-exponent))
    #The weiner coefficents are essentially the inverse of the noise contribution
    #scaled by our theoretical signal
    weiner_coefficents=1/(1+np.absolute(noise_fourier)**2/fake_sig)
    weiner_coefficents[0]=1
    #We now scale the signal by these values, decreasing the higher frequencies
    #based on their inclusion in the noise
    fil_fourier=weiner_coefficents*unfil_fourier
    #Finally we inverse fourier transorm to get a filtered signal. The ifft may
    #sometimes leave small amounts of imaginary components due to rounding error,
    #this is taken out by taking the real portion.
    filtered=(np.real(fft.ifft(fil_fourier))).tolist()

    return filtered
