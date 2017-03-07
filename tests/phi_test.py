#! /usr/bin/env python 
import scipy as __sci
import scipy as sci
import numpy as __np
import numpy as np

import pylab as py

import numpy.fftpack as __ft
import numpy.fftpack as ft

import Charm as chr
from Charm import *
import cPickle as cpk
import time as tm

#============================


f=chr.gmanf(1000,307,11.6,-3.54,1.9)
df=f.data(tnorm=1)[0]

df2=chr.phase_shift(df,phi=-0.9,tnorm=1)

g=f.copy()
g.phi=1
dg=g.data(tnorm=1)[0]

py.plot(df2,'b',lw=1)
py.plot(dg,'r',lw=1)

py.show()
