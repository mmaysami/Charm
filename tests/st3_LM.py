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
import os.path as path
import cPickle as cpk
import copy  as cp
import Charm as chr
import Charm.Core.Misc as misc
import scipy as sci
import numpy as np
import numpy.fft as ft
import time as tm
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
N=1000
fi0=0; NRM,dNRM=1,1; dphi=-1j
niter=500

con=0* np.array([0.75,0.75,0.75,0.75,0])
test='OW'  # 'gmanf' for gmanf2gmanf match or OW from windowed event


ind=4
param=[238,18,-3.78,0.34]   # actual attributes
numb=[1,2]

est_dict={'prm_ind':numb,'dphi':dphi,'flatphi':1,'niter':niter,'intg':fi0,'verb':0,'show':0,'smth':0}

#---------------------------------------
s0=[0]
tol=[5e-8]

print "\n Smoothness Level:"+str(s0)+" & Inverting for ",numb

guess=np.array(param) 
guess.put(numb,[10.7,-1.3])
I  = chr.gmanf(n=N); I.update(param)
In = I.data(tshift=0,tnorm=NRM)[0]

if test != 'gmanf':
    param=attrib_org[ind,:]
    guess=attrib_det[ind,:]
    guess[2]=-1
    In=events_det[ind,:]

        
runtime=0
tm_est=tm.time()
print "\n Actual Parameters ",param
for i in range(len(s0)):

    if s0[i]==s0[-1]:
        est_dict['show']=1
        est_dict['verb']=1
    est_dict.update({'smth':s0[i],'tol':tol[i]})
    (attrib_est,err,f,ind)=chr.LMestimate(events_det,guess,ind=ind,**est_dict)
    guess=0.5* (guess+attrib_est)

tm_est=tm.time()-tm_est
print "\n Estimation time in second(s): ",tm_est

