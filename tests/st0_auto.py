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
import Charm as chr
import scipy as sci
import numpy as np
import pylab as py
import numpy.fft as ft
import cPickle as cpk
import copy as cp
import time as tm
from sys import stderr

#============ Parameters =============
niter=120    # Number of iters
dNRM=1      # Makes H=I if 1 
NRM=1       # Makes MANIFOLD.data normalized if 1
tol=1e-6
VERB=0
N=1000

df_fi=0
smth=5
dphi=-1j

param=[238.0,15.0,-2.78,0.25]   # actual attributes
act_param=cp.copy(param)

# Initial Guess
numb=[1,2]
guess=np.array(act_param) 
guess.put(numb,[11.0,-1.5])
numb=[1,2]
#---------------------------------------
print "Inverting for ",numb
I  = chr.gmanf(n=N)
I.update(param)
In = I.data(tshift=0,tnorm=NRM)[0]


VERB=1;SHOW=0

tm_LSest=tm.time()
(attrib_est,err,f,ind)=chr.LSestimate(In,guess,ind=0,prm_ind=numb,smth=smth,intg=df_fi,dphi=dphi,flatphi=1,niter=niter,tol=tol,verb=VERB,show=SHOW,log=stderr)
tm_LSest=tm.time()-tm_LSest

tm_LMest=tm.time()
(attrib_est,err,f,ind)=chr.LMestimate(In,guess,ind=0,prm_ind=numb,dphi=dphi,flatphi=1,niter=2*niter,tol=tol,verb=VERB,show=SHOW,log=stderr)
tm_LMest=tm.time()-tm_LMest

tm_BFGSest=tm.time()
(attrib_est2,err2,f2,ind2)=chr.BFGSestimate(In,guess,ind=0,prm_ind=numb,dphi=dphi,flatphi=1,niter=10*niter,tol=tol,verb=VERB,show=SHOW,log=stderr)
tm_BFGSest=tm.time()-tm_BFGSest

print "\n BFGS Estimation time iin second(s): ",tm_BFGSest
print "\n\n LS   Estimation time iin second(s): ",tm_LSest
print "\n\n LM   Estimation time iin second(s): ",tm_LMest
