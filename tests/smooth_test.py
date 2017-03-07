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
scl=5
fi=0


f=chr.gmanf(1000,307,11.6,-3.54,1.9)
df=f.data(tnorm=1)[0]
noise=0.005*sci.randn(1000)* max(abs(df))
# noise-=np.mean(noise)
df+=noise*0

g=f.copy(s=scl,fi=fi)
dg=g.data(tnorm=1)[0]
dfs=chr.smooth(df,s=scl,fi=fi,tnorm=1)[0]


f.fi=fi
f.sigma=np.sqrt(f.sigma**2+scl**2)
df2=f.data(tnorm=1)[0]

py.plot(df2,'g',lw=1)
py.plot(dg,'--b',lw=1)
py.plot(dfs,'--r',lw=1)

py.show()
