"""
This file contains an example script for analyzing output from the MQDS package
"""
import numpy as np
import analysis as anls
import mqds_class as mqds
import matplotlib.pyplot as plt

method = input( "Please insert method used to calculate the linear \
response function formatted as t(fs), Re( R(t) ), Im( R(t) ) :  " )

data = mqds.LinearResponse(method)

# retrieve the lower bound on frequency domain over which to ft
wmin = input( "Please insert the lower bound on frequency domain (cm-1):  " ) 
wmin = float(wmin)
# retrieve the upper bound on frequency domain over which to ft
wmax = input( "Please insert the upper bound on frequency domain (cm-1):  " ) 
wmax = float(wmax)
# retrieve the number of discrete points over which to fourier transform
wpts = input( "Please insert the number of discrete frequency points :  " ) 
wpts = int(wpts)
tpts = int(len(data.time))

real_response = anls.apply_cos_window(data.real) 
imag_response = anls.apply_cos_window(data.imag) 

w = data.wrange(wmin, wmax, wpts)

test = anls.one_d_ft(data.real, data.imag, data.time, w)

plt.plot(w, -test.imag )
plt.show()
