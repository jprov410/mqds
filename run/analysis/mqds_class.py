import os
import glob
#import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt

__all__ = ['LinearResponse','ThirdOrderResponse','DensityMatrix']

class DensityMatrix(object):
    """
    class for mqds density matrix data
    """
    def __init__(self, method):
        """
        takes reduced density matrix method as input
        e.g.) pldm, sqc, ipldm
        """
        prefix = method + '.'
        self.time = None
        size = int( np.sqrt( len( glob.glob( prefix+'*' ) ) ) )
        for filename in os.listdir('.'):
            if filename.startswith(prefix):
                if self.time is None:
                    with open(filename, 'r') as f:
                        self.time = []
                        for line in f:
                            self.time.append( float( line.split()[0] ) )
                    break    
        self.time = np.array( self.time )
        self.rho = np.empty([size,size,len(self.time)], dtype=complex)
        for filename in os.listdir('.'):
            if filename.startswith(prefix):
                suffix = filename.split('.')[1]
                index1, index2 = suffix.split('-')
                index1, index2 = int(index1), int(index2)
                self.rho[index1-1,index2-1] = self.matrix_element( filename )
        self.shape = self.rho.shape

    def matrix_element(self, filename):
        """
        Takes filename as an argument to retrieve data for each density matrix element
        """
        real, imag = [], []
        with open(filename,'r') as f:
            for line in f:
                real.append( float( line.split()[1] ) )
                imag.append( float( line.split()[2] ) )
        real = np.array( real )
        imag = np.array( imag )
           
        return real + 1j * imag

    def __getitem__(self, i):
        return self.rho[i]

class Response(object):
    """
    Superclass for response functions calculated with mqds
    """
    def ___init___(self):
        self.type = 'Response Function'

    def wrange(self, wmin=0.0, wmax=500.0, wpts=250):
        """
        defines frequency range over which to Fourier transform the
        linear response function
        """
        w = []
        for i in range(0, wpts+1):
            w.append( float( i * (wmax - wmin) / wpts + wmin ) )
        w = np.array( w )
        return w

class LinearResponse(Response):
    """
    class for mqds linear response function data
    """
    def __init__(self, method):
        """
        takes method used to compute the response function
        as an argument
        """
        infile = method + '_linrespfunc.out'
        self.time, self.real, self.imag = [], [], []
        with open(infile, 'r') as f:
            for line in f:
                self.time.append( float( line.split()[0] ) )
                self.real.append( float( line.split()[1] ) )
                self.imag.append( float( line.split()[2] ) )
                
        self.time = np.array( self.time )
        self.real = np.array( self.real )
        self.imag = np.array( self.imag )


class ThirdOrderResponse(Response):
    """
    class for mqds nonlinear response function data
    """
    def __init__(self, method):
        """
        takes filename that contains the response function 
        as an agument
        """
        prefix = method + '_nonlin'
        self.time1, self.time2, self.time3 = [], [], [] 
        self.real, self.imag = [], []
        with open(infile, 'r') as f:
            for line in f:
                self.time1.append( float( line.split()[0] ) ) 
                self.time2.append( float( line.split()[1] ) ) 
                self.time3.append( float( line.split()[2] ) ) 
                self.real.append( float( line.split()[3] ) ) 
                self.imag.append( float( line.split()[4] ) ) 

        self.time1 = np.array( self.time1 )
        self.time2 = np.array( self.time2 )
        self.time3 = np.array( self.time3 )
        self.real = np.array( self.real )
        self.imag = np.array( self.imag )
