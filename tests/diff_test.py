#! /usr/bin/env python 
import scipy as __sci
import scipy as sci
import numpy as __np
import numpy as np

import pylab as py

import numpy.fft as __ft
import numpy.fft as ft

import Charm as chr
from Charm import *
import cPickle as cpk
import time as tm

#============================
py.close('all')
# py.figure()
# for ph in np.arange(0,1,0.1):
#     f=chr.gmanf(1000,512,10,-1,ph)
#     f.plot()
# py.show()

dNRM=0
const=-1j

f=chr.gmanf(1126,579,126.066,-128.07,1.4764)
# f=chr.gmanf(1000,200,10,-1.1,.505)
print f.__param__()
py.figure()
for k in range(5):
    py.subplot(510+k+1)
    if k==0:
        py.plot(f.data(tnorm=dNRM)[0],lw=1)
        py.plot(np.imag(f.data(tnorm=dNRM)[0]),'r',lw=1)
    elif k==4:
         py.plot( f.eval(att='d'+str(k-1))(dphi=const,tnorm=dNRM)[0],lw=1)
         py.plot(np.imag(f.eval(att='d'+str(k-1))(dphi=const,tnorm=dNRM)[0]),'r',lw=1)

    else:
        py.plot( f.eval(att='d'+str(k-1))(tnorm=dNRM)[0],lw=1)
        py.plot(np.imag(f.eval(att='d'+str(k-1))(tnorm=dNRM)[0]),'r',lw=1)


py.show()



py.figure()
G=f.grad(tnorm=dNRM)[0]
for k in range(5):
    py.subplot(510+k+1)
    if k==0:
        py.plot(f.data(tnorm=dNRM)[0],lw=1)
        py.plot(np.imag(f.data(tnorm=dNRM)[0]),'r',lw=1)
    else :
        py.plot(G[:,k-1],lw=1)
        py.plot(np.imag(G[:,k-1]),'r',lw=1)
py.show()












# f=chr.gmanf(1000,500,10,-3.4,0)
# py.figure()
# rng=[0,1,2,3,4]
# # Imaginary parts are all zero
# for k in rng:
#     if k==0:
# #         py.plot(f.data(tnorm=1)[0],lw=1)
#         py.plot(np.imag(f.data()[0]),'r',lw=1)
#     else:
# #         py.plot( f.eval(att='d'+str(k-1))(tnorm=1)[0],lw=1)
#         py.plot(np.imag(f.eval(att='d'+str(k-1))()[0]),'r',lw=1)

# py.legend()
# py.show()
