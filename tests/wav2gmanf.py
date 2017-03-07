#! /usr/bin/env python 
import scipy as __sci
import scipy as sci
import numpy as __np
import numpy as np

import pylab as py

import numpy.fft.fftpack as __ft
import numpy.fft.fftpack as ft

import Charm as chr

import cPickle as cpk
import time as tm
#========================================
# Compare phase in synthesis and analysis
# 
#========================================
N=1000
alpha=-1.90
# Up to +1.9 gmanf is fine, fix for a>2 is  DONE !
phi=  -(alpha)/2
tau=2*N/5.
scl=17


#                                                              __  
# Positive Waveform is used with minus sign (@ zero has max \_/  \_/ )
# Same as trace generator function
spf=chr.fspline(N,alpha=alpha,phi=phi,loc=tau,tnorm=0,domain=None)
spt=np.fft.ifft(spf)/ chr.norm( np.fft.ifft(spf) )
src= chr.wavelet('gaus',wtype='real',param={'p':2})
srcf= src.wavf(N,s=scl,domain=None,tnomr=0)[0]
wav=chr.outtype(spf*srcf,domain='time',tnorm=1)[0]

gobj=chr.gmanf(N,tau,scl *1.4 ,alpha-2.00000,phi)
gwav=gobj.data(domain='time',mode='complex',tnorm=1)[0]



py.clf()
py.plot(wav,'b',lw=1)
py.plot(gwav,'--r',lw=1)

py.show()

