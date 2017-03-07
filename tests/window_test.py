
import scipy as sci
import numpy as np
import pylab as py
import numpy.fft as ft
import scipy.integrate as integ
__sci=sci
__np=np
__py=py
import Charm as chr
import cPickle as cpk
import time as tm
import pdb
#===================================================
# Windowing
#===================================================

N=500
tau=loc=250
phi=0.83

sigma=4.2
alpha=-1.74

f=chr.gmanf(N,tau,sigma,alpha,phi);data=f.data()[0]; data/=max(abs(data))

if sigma<3:
    width=60
    ordl=ordr=4
    left=right=0.1
else:
    width=50
    ordl=ordr=6
    left=right=0.1


width*=sigma
t=__sci.arange(0,N,dtype=int)
loc=int(loc)

"# Left & Right half of the window function"
# ordl*=sigma
# ordr*=sigma
locl=left*width
locr=right*width
b0= (1 +(t/locl)**(2*ordl)) **-1 
bn=(1 +(t/locr)**(2*ordr))**-1


mask=__sci.zeros(N)
# Check Left & Right part of window 
mask[loc:0:-1]=b0[0:loc]
mask[loc:N]=bn[0:N-loc]
lenl=loc
lenr=N-loc

wdata=data*mask

py.plot(data,'y');py.plot(wdata,'--b');py.plot(mask,'--g')
py.show()




W=chr.wbutterworth(data,tau,sigma,ordl=6,ordr=6,width=65)
# py.plot(data,'y');py.plot(W[0],'--g');py.plot(W[1],'--k')

"""
# Windowing test

# W=chr.wboxcar(data,tau,sigma,lhalf=5,rhalf=5)
# W=chr.wgaussian(data,tau,sigma,scl=4)

# W=chr.wblackman(data,tau,sigma,left=3,right=3,taper=3)
"""

