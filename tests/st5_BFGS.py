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
niter=500
fi0=0; NRM,dNRM=1,1; dphi=-1j
s0=0;tol=5e-8
test='gmanf'  # 'gmanf' for gmanf2gmanf match or OW from windowed event


ind0=3
param=[238.0,18.0,-3.78,0.34]   # actual attributes
numb=[0,1,2,3]


est_dict={'prm_ind':numb,'dphi':dphi,'flatphi':1,'niter':niter,'intg':fi0,'verb':1,'show':1,'smth':s0,'tol':tol}

#---------------------------------------
print "\n Smoothness Level:"+str(s0)+" & Inverting for ",numb

guess=np.array(param) 
#numb=[1,2]
guess.put(numb,[225,10.7,-1.3,0.25])
I  = chr.gmanf(n=N); I.update(param)
In = I.data(tshift=0,tnorm=NRM)[0]
events=chr.vector(In,'row')
ind=0

if test != 'gmanf':
    ind=ind0
    param=attrib_org[ind,:]
    guess=attrib_det[ind,:]
    guess[2]=-1
    events=events_det
    In=events[ind,:]/chr.norm(events[ind,:])

runtime=0
tm_est=tm.time()
print "\n Actual Parameters ",param
attrib_est,err,f,ind=chr.LBFGSestimate(events,guess,ind=ind,**est_dict)
tm_est=tm.time()-tm_est
print "\n Estimation time in second(s): ",tm_est



#========================
# GRAPHS
#========================
from   sys   import stderr
import os 
py.close('all')
chr.Core.Steps.__est_show("LBFGS",ind,numb,5,guess,attrib_est,In,f,
           err,NRM,verb=1,show=1,sepfig=1,log=stderr)

prefix2=chr.__path__[5]+'/Example2/'
py.figure(1);py.savefig(prefix2+'EstI.eps')
py.figure(2);py.savefig(prefix2+'EstF.eps')
