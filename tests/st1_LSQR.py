#! /usr/bin/env python 
"""
2D Nonlinear parametric inversion for Estimation of seismic singularities


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
import SLIMsoft as ss
import numpy.fft as ft
import cPickle as cpk
import time as tm

py.close('all')
#============ Parameters =============
VERB=0
FIG=0

N=1024
niter=50    # Number of iters
dNRM=1      # Makes H=I if 1 
NRM=1       # Makes MANIFOLD.data normalized if 1
thr=1e-6 
dphi=1j
param=[tau,sigma,alpha,phi]=[238,15,-2.78,0.25]  # Actual Values
# Initial Guess
guess=np.array(param) #[258,6,-1.8,1.5]
numb=[1,2]
guess.put(numb,[16,-0.1])

USEMIN=False  # True to STOP when reaches Min. Error 
dphi=1j

# ti,si : very good
# ti,ai : ti very close OR  wierd wrt ai (248,-1.3..5 works but >-1.8 NO) !
# si,ai : ~good

# ti,pi : Need switch between [+j]mostly, -j and Wierd
# si,pi : [-j]mostly, and si <sigma
# ai,pi : [+j]mostly, Wierd 

#---------------------------------------
#py.close('all')
err=list()
Theta=list()
Theta.append(guess.take(numb))


I  = chr.gmanf(n=N,tau=tau,sigma=sigma,alpha=alpha,phi=phi)
In = I.data(tshift=0,tnorm=NRM)[0]
[ti,si,ai,pi]=guess
f=chr.gmanf(n=N,tau=ti,sigma=si,alpha=ai,phi=pi)
err.append(chr.mse(f.data(tshift=0,tnorm=NRM)[0] , In ))

# Show Initial Statistics
print '\n  --< Initial Guess >--'
print '(tau,sigma,alpha,phi) = ',(tau,sigma,alpha,phi)
print '( ti , si , ai,pi  ) = ',guess
print '   Error   : ',err[-1],'\n'
nfig=py.figure()
nfig=nfig.number
py.plot(In,'--b',f.data(tshift=0,tnorm=NRM)[0],'r',lw=1)
py.legend(['Image','Gaussian Estimator'])
py.show()
tm.sleep(0.1)


# Computing Jacobian & Hessian (Column Vect.)
E= np.matrix(In - f.data(tshift=0,tnorm=NRM)[0]).T
G=f.grad(tshift=0,tnorm=dNRM,dphi=dphi)[0]
G=np.matrix(G.take(numb,axis=1))
# b=2*(G.T * E)
# A=2*np.matrix (sci.diag(sci.diag( G.T * G  )))

b=E
A=G  
# Determining new parameter values
iters=0
dTheta,res2,rank,s=np.linalg.lstsq(A,b)
# dTheta,res2,iters=ss.solvers.solvelsqr(A, b, x0=None, maxiter=niter, tol=1e-15*thr)

Theta.append( Theta[-1]+np.array(dTheta).ravel()  )
guess.put(numb,Theta[-1])
guess[3]%=2
f.update(guess)
err.append(chr.mse(f.data(tshift=0,tnorm=NRM)[0],In))


# Show Final Results
Theta=sci.array(Theta)
err=sci.array(err)
print "--------------------------- \n"
print ' --< Final Values >--' 
print 'Iteration # ',iters
print '(tau,sigma,alpha,phi) = ',(tau,sigma,alpha,phi)
print '(  ti , si,  ai , pi ) = ',guess
ind=np.argmin(err) 
minguess=guess
minguess.put(numb,Theta[ind])
minguess[3]%=2
print 'Min  error Parameters  = ',minguess
print 'Final Error Value =',err[-1]
print ' Min. Error       =',min(err)

g=f.copy()
g.update(minguess)
py.figure(nfig)
py.clf()
py.plot(In,'--b',f.data(tshift=0,tnorm=NRM)[0],'r',g.data(tshift=0,tnorm=NRM)[0],'g',lw=1)
py.title('Image & Gausssian Manifold Estimation',size=16)
py.legend(['Image','Gaussian Estimator','Min. Error Estimate'])
py.show()
