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
from sys import stderr as __stderr
from scipy import optimize as __opt
import Charm as chr
import scipy as sci
import numpy as np
import pylab as py
import numpy.fft as ft
import cPickle as cpk
import copy as cp
import time as tm

# py.close('all')
prefix=chr.__path__[1]
# wtraces=sci.load(prefix+'/events_det.pyd')
#============ Parameters =============
niter=120    # Number of iters
dNRM=1      # Makes H=I if 1 
NRM=1       # Makes MANIFOLD.data normalized if 1
thr=1e-6
VERB=0
N=1000


df_fi=0
integ=0
if integ ==0:
    df_fi=0


param=[238,15,-2.78,0.25]   # actual attributes
act_param=cp.copy(param)
if integ:
    act_param[2]+=1;act_param[3]=(act_param[3]-0.5)%2

# Initial Guess
numb=[1,2]
guess=np.array(act_param) 
guess.put(numb,[16,-0.1])
numb=[1,2]



dphi=1j       # Useful only when inverting for phase 
STEP=1        # Find step size  for firection "d"
USEMIN=False  # True to STOP when reaches Min. Error 
"""
# ti,si : very good
# ti,ai : ti very close OR  wierd wrt ai (248,-1.3..5 works but >-1.8 NO) !
# si,ai : ~good

# ti,pi : Need switch between [+j]mostly, -j and Wierd
# si,pi : [-j]mostly, and si <sigma
# ai,pi : [+j]mostly, Wierd 
"""

#---------------------------------------
iter=0
err=list()
Theta=list()
Theta.append(guess.take(numb))


I  = chr.gmanf(n=N)
I.update(param)
In = I.data(tshift=0,tnorm=NRM)[0]
Is=chr.smooth(In,s=0,fi=df_fi,tnorm=NRM)[0]
if integ:
    Is =  np.cumsum(In)
    Is -= np.mean(Is)
    Is /= chr.norm(Is)


"""
#Pick segmented event to be characterized
# ind=4
# In= wtraces[ind,:]     
# if NRM==1:
#     In/= chr.norm(In)
"""

[ti,si,ai,pi]=guess.copy()
f=chr.gmanf(n=N,tau=ti,sigma=si,alpha=ai,phi=pi,fi=df_fi)
err.append(chr.mse(f.data(tshift=0,tnorm=NRM)[0] , Is ))

cond=(err[iter]<=err[max(0,iter-2)] or iter<niter/5.)
while (err[iter]>thr and (cond or not(USEMIN)) and iter<niter ):

    iter+=1
    # Computing Jacobian & Hessian (Column Vect.)
    E= np.matrix(f.data(tshift=0,tnorm=NRM)[0] - Is).T
    G=f.grad(tnorm=dNRM,tshift=0,dphi=dphi)[0]
    G=np.matrix(G.take(numb,axis=1))
    J=(G.T * E)
    H=np.matrix (sci.diag(sci.diag( G.T * G  )))
    #     H=2*(G.T * G)

    # Determining new parameter values (Direction:d & Step:a0)
    mu=0.1
    a0=1
    d=np.matrix(H.I * J)
    guess.put(numb,Theta[iter-1]-a0*d)
    guess[3]%=2
    f.update(guess)
    E_new=f.data(tshift=0,tnorm=NRM)[0]-Is
    # Find step size along direction d
    while chr.norm(E_new)  > chr.norm(E+mu*a0*E.T*G*d) and STEP and a0>5e-4:
        a0=0.5*a0
        guess.put(numb,Theta[iter-1]-a0*d)
        guess[3]%=2
        f.update(guess)
        E_new= f.data(tshift=0,tnorm=NRM)[0]-Is

        while __misc.norm(E_new) > __misc.norm(E+mu*a0*E.T*G*d) and a0>1e-3:
            a0 *= 0.50
            param.put(prm_ind,Theta[iter-1]-a0*d)
            if flatphi:
                param[3]%=2
            # Print & Show results            f.update(param)
            E_new= f.data(tshift=0,tnorm=NRM)[0]-In


    Theta.append( Theta[iter-1] - a0*d )
    err.append(chr.mse(f.data(tshift=0,tnorm=NRM)[0],Is))

    # Condition for optimization (Usedwhen USEMIN=1)
    cond=(err[iter]<=3 * err[max(0,iter-2)] or iter<niter/5)
"""
#     # Print and show results of current step
#     if VERB:
#         print 'iter #',iter
#         print 'err',err[iter]
#         print '( ti , si , ai , pi ) = ',guess, '\n'
"""











#===================================
#       Statistics & Results
#===================================
print " Actual Value  ",act_param
print " Initial Guess ",guess
print " Estimation    ",Theta[-1]
print " Error ", err[-1]


ind=np.argmin(err) 
minguess=guess.copy()
minguess.put(numb,Theta[ind]);minguess[3]%=2
f.update(guess)
g=f.copy()
g.update(minguess)
py.figure()
py.clf()
py.plot(Is,'-b',g.data(tshift=0,tnorm=NRM)[0],'-r')#,f.data(tshift=0,tnorm=NRM)[0],'--g',lw=1.5)
py.title('Image & Gausssian Manifold Estimation',size=16)
py.legend(['Image / Cumsum:'+str(integ),'Min Estimator'])#,'Estimation'])
py.show()
