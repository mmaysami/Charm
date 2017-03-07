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
from os.path import join
#============================
py.close('all')
prefix=chr.__path__[3]

mymodel=open(join(prefix,'synt_model.pyd'),'r+')
model=cpk.load(mymodel)
mymodel.close()
trace=model['trace']
events=model['events']
attrib_org=model['attrib_org']
events_det=sci.load(join(prefix,'events_det.pyd'))


`ind=0
N=model['N']
print 'Attributes:',attrib_org[ind,:]
I=events_det[ind,:]
In=I/chr.norm(I)

# Create a range of values for alpha and sigma around actua values
al=np.linspace(-2.5,2.5,100)
al+=attrib_org[ind,2]
scl=np.linspace(-15,15,1000)
scl+=attrib_org[ind,1]
scl=scl[scl>0]

x=8.13147785426

misfit=np.zeros(x.shape)
f=chr.gmanf(N)
f.update(attrib_org[ind,:])
# I=f.data(tnorm=1)[0]
f.alpha+=1.5
# f.sigma+=7
f.__param__()
for t in range(len(x)):
    f.sigma=x[t]
    misfit[t]=chr.mse(f.data(tshift=0,tnorm=1)[0] , In)


py.plot(x,misfit);py.show()
