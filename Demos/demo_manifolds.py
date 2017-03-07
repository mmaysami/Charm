#!/usr/bin/env python
import os,sys
import numpy.fft as ft
import numpy as np
import scipy as sci
import pylab as py
import Charm as chr
from Charm import Misc as misc

tsize=py.rcParams['axes.titlesize']
lsize=py.rcParams['axes.labelsize']

dir = chr.__path__[5]
prefix=os.path.join(dir,'Intro')
if not os.access(prefix,os.F_OK):
    os.makedirs(prefix)

py.close('all')
N=1000
tau=500
sigma=15
alpha=0
phi=0
#---------------
tau1=230
sigma1=60
alpha1=-1.5
phi1=0.25

f=chr.gmanf(N,tau,sigma,alpha,phi)

py.figure()
f.plot(label=False,tamp=0,tnorm=1,lw=2,c='b')
py.title("Reference object",size=tsize)
py.xlabel('[tau,sigma,alpha,phi] ='+str([f.tau,f.sigma,f.alpha,f.phi]),size=lsize)
py.savefig(os.path.join(prefix,"gman_ref.eps"))

f.tau=tau1
py.figure()
f.plot(label=False,tamp=0,tnorm=1,lw=2,c='k')
py.title("Shift",size=tsize)
py.xlabel('[tau,sigma,alpha,phi] ='+str([f.tau,f.sigma,f.alpha,f.phi]),size=lsize)
f.tau=tau
py.savefig(os.path.join(prefix,"gman0_tau.eps"))

f.sigma=sigma1
py.figure()
f.plot(label=False,tamp=0,tnorm=1,lw=2,c='r')
py.title("Shrinkage/Expansion",size=tsize)
py.xlabel('[tau,sigma,alpha,phi] ='+str([f.tau,f.sigma,f.alpha,f.phi]),size=lsize)
f.sigma=sigma
py.savefig(os.path.join(prefix,"gman1_sigma.eps"))

f.alpha=alpha1
f.phi=-alpha1/2.
py.figure()
f.plot(label=False,tamp=0,tnorm=1,lw=2,c='m')
py.title("Fractional differentiation",size=tsize)
py.xlabel('[tau,sigma,alpha,phi] ='+str([f.tau,f.sigma,f.alpha,f.phi]),size=lsize)
f.phi=phi
f.alpha=alpha
py.savefig(os.path.join(prefix,"gman2_alpha.eps"))


f.phi=phi1
py.figure()
f.plot(label=False,tamp=0,tnorm=1,lw=2,c='g')
py.title("Instantaneous phase",size=tsize)
py.xlabel('[tau,sigma,alpha,phi] ='+str([f.tau,f.sigma,f.alpha,f.phi]),size=lsize)
f.phi=phi
py.show()
py.savefig(os.path.join(prefix,"gman3_phi.eps"))
