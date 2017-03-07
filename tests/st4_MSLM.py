#! /usr/bin/env python 
"""
Multi-Dimensional Nonlinear parametric inversion for Estimation of seismic singularities


  Author:  Mohammad Maysami, SLIM EOS UBC
           Seismic Laboratory for Imaging and Modeling
           Department of Earch & Ocean Sciences
           The University of British Columbia

  Copyright (C) 2006 The University of British Columbia at Vancouver
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

"""
from sys import stderr,stdout
import cPickle as cpk
import copy  as cp
import Charm as chr
import Charm.Core.Misc as misc
import scipy as sci
import numpy as np
import numpy.fft as ft
import time as tm
from os import path
try:
    import pylab as py
except ImportError:
    from Charm import Pylab_Error
    py=None
    

prefix=chr.__path__[4]
mymodel=open(path.join(prefix,'synt_model.pyd'),'r+')
model=cpk.load(mymodel)
mymodel.close()
attrib_org=model['attrib_org']
events_det=sci.load(path.join(prefix,'events_det.pyd'))
attrib_det=sci.load(path.join(prefix,'attrib_det.pyd'))        
py.close('all')

#=====================================
#            Estimation
#=====================================
fi0=0; NRM,dNRM=1,1; dphi=-1j

niter=500   # Number of iters
ind=0
numb=[1,2]

est_dict={'ind':ind,'prm_ind':numb,'dphi':dphi,'flatphi':1,'niter':niter,'intg':fi0,'show':1}

#---------------------------------------
# Initial Guess & Events
param=attrib_org[ind,:]
guess=attrib_det[ind,:]
guess[2]=-1.2
In=events_det[ind,:]


#---------------------------------------
smth=[10,5,0]
tol=[5e-6,1e-6,5e-7,1e-7,5e-8]


print "\n Smoothness Level:"+str(smth)+" & Inverting for ",numb
print "\n Actual Parameters ",param
tm_est=tm.time()
attrib_est,err,f,ind=chr.MSLMestimate(events_det,guess,ind,smth,tol,verb=1,show=1,LM_dict=est_dict)
tm_est=tm.time()-tm_est
print "\n Estimation time in second(s): ",tm_est


