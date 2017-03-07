#!/usr/bin/env python
import os,sys
import numpy as np
import numpy.fft as ft
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
loc=200
scl=16
pad=10
src=chr.wavelet()

#===========================
# Show Transition-Reflectivity & Trace
#      for different alphas
#===========================
step=0.1
kt=np.arange(0,1+step,2*step)
lw=sci.ones(kt.shape,dtype='int');lw[0]=lw[-1]=2
cl=['b','g','m','y','c','r']
# for k in range(len(kt)-1,-1,-1):
#     py.figure(1)
#     py.plot(chr.spline(N,kt[k],loc),lw=lw[k],c=cl[k])
#     ref=chr.fspline(N,kt[k],phi=None,loc=N/2,domain='freq')
#     w=chr.gauswavf(N,scl,domain='freq')[0]
#     wave=ft.ifft(w*ref)
#     py.figure(2)
#     py.plot(ft.ifft(ref),lw=lw[k],c=cl[k])
#     py.figure(3)
#     py.plot(wave/chr.norm(wave),lw=lw[k],c=cl[k])


# py.figure(1)
# py.axis([0,N,-0.05,1.05])
# py.title('Different order of singularities (0 to 1, step=0.1)',size=tsize)
# py.xlabel('Singularity orders = '+str(kt),size=lsize)
# py.ylabel('Fluctuation ',size=lsize)
# py.savefig(os.path.join(prefix,"fractional_splines.eps"))

# py.figure(2)
# py.axis([0,N,-1.2,1.2])
# py.title('Model for reflectivity',size=tsize)
# py.xlabel('Singularity orders = '+str(kt),size=lsize)
# py.ylabel('Amplitude',size=lsize)
# py.savefig(os.path.join(prefix,"Transitions.eps"))

# py.figure(3)
# py.title('Seismic signal with Ricker',size=tsize)
# py.xlabel('Transition orders = '+str(kt),size=lsize)
# py.ylabel('Amplitude',size=lsize)
# py.savefig(os.path.join(prefix,"seismic_events.eps"))
#==================================
# Transitin and events in same fig
#==================================
scl=32
pad=5
for k in range(len(kt)-1,-1,-1):
    ref=chr.spline(N,kt[k],N/2)
    Fref=chr.fspline(N,kt[k],phi=None,loc=N/2,domain='freq')
    w=chr.gauswavf(N,scl,domain='freq')[0]
    wave=ft.ifft(w*Fref)
    py.figure(5)
    py.subplot(121)
    py.plot(ref,lw=lw[k],c=cl[k])
    py.subplot(122)
    py.plot(wave/max(abs(wave)),lw=lw[k],c=cl[k])
    
py.figure(5)
py.subplot(121)
py.axis([0,N,-0.05,1.05])
py.title('Model for Reflections',size=tsize)
py.ylabel('Amplitude',size=lsize)
py.xlabel('Transition orders = '+str(kt))
py.subplot(122)
py.title('Seismic signals with Ricker',size=tsize)
# py.ylabel('Amplitude',size=lsize)
py.savefig(os.path.join(prefix,"spline2seismic.eps"))

#=======================================
# Spiky Deconv. for Non Spiky reflectors
#=======================================
N=1000;K=8;Delta=60
s0=10                      # ==> Scale of source signature    5
amp=(4,7)                  # Amp Ranges 
rnga=(-2,1)              # Alpha ranges
rngp=(0,2)
rngs=(8,12)                # Scale range for gaussian trace generato
src=chr.wavelet('gaus',param={'p':2})

(trace,events,attrib_org,attrib_vect,k,delta,regul,model_dict)=chr.rand_trace(N,K,Delta,scl=s0,rnga=rnga,rngp=rngp,src=src,amp=amp,regul=0,uniform=0)
spk=attrib_vect[0,:]# np.sign(2*np.random.rand(N)-1)
spk[spk==chr._Null]=0
py.figure()
t=chr.src_conv(spk,scl=s0)
t2=ft.ifft(ft.fft(spk)*src.wavf(N,s0,domain='freq')[0])
py.plot(trace/max(trace),'--g',lw=2)
py.plot(t2/max(t2),'b',lw=2)
py.plot(spk/max(spk),'r',lw=2)
#### Includes phase and alpha both 
py.legend(['Seismic trace','Spiky reflections','Spiky Model'])
py.title('Spiky deconv. Modeling',size=tsize)
py.xlabel("Alphas ="+str(np.int32(100*attrib_org[:,2])/100.),size=lsize)
py.ylabel('Amplitude',size=lsize)
py.savefig(os.path.join(prefix,"CWT_trace.eps"))
py.show()
##========================================
# #CWT  MML (Use after generation of data)
##========================================
# wr=chr.wavelet('gaus',wtype='complex',param={'p':2})
# scl0=2**np.linspace(1.5,5)
# (mml,cw,lnpos,ext)=wr.mmle(trace,s=scl0)
# acw=abs(cw)
# ind=lnpos[1]
# d=acw[ind[:,0],ind[:,1]]


# py.jet()
# py.matshow(py.rot90(acw-mml*(mml>0)*acw));#py.colorbar()
# py.title('Modulus of CWT & MMLs',size=tsize)
# py.ylabel('Scale index',size=lsize)
# py.xlabel('Location',size=lsize)
# py.savefig(os.path.join(prefix,"CWT_MML.eps"))


# py.figure()
# py.loglog(d[::-4],lw=2);py.grid()
# py.title('Modulus of CWT along MML',size=tsize)
# py.ylabel('log |W|',size=lsize)
# py.xlabel('log(Scale index)',size=lsize)
# # py.xlabel('Scale index (Coarse -> fine)',size=lsize)
# py.show()
# py.savefig(os.path.join(prefix,"CWT_ACWonMML.eps"))

# py.figure()
# py.plot(trace,'b',lw=2)
# py.title('Seismic trace',size=tsize)
# py.ylabel('Amplitude',size=lsize)
# py.savefig(os.path.join(prefix,"CWT_trace.eps"))



# #======================================
# #  Random layers of subsurface
# #======================================
# x=np.arange(0,2000)/2000.
# y=np.zeros((4,len(x)))
# y[0,:]=2*x**5-x**3-x**2-0.5*x+1
# y[1,:]=x**3-4*x**2+x+0
# y[2,:]=-x**4-3*x**3+2*x**2-1.5
# y[3,:]=x**4+x**3-4*x**2+x-3
# py.figure()
# for k in range(4):
#     py.plot(x,y[k,:])


# x=np.arange(-1000,1000)/1000.
# y=np.zeros((4,len(x)))
# y[0,:]=2*x**5-x**3-x**2-0.5*x+1
# y[1,:]=-x**3*(x<0)+x**3-2*x**2+0.2*x
# y[2,:]=-5*x**4-3*x**3+2*x**2+2*x+1.5*(x-0.5)*(x>0.5)-2
# y[3,:]=x**4-x**3-4*x**2+x-3
# py.figure()
# for k in range(4):
#     py.plot(x,y[k,:])
