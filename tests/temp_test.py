#! /usr/bin/env python 

import scipy as sci
import numpy as np
import pylab as py
import numpy.fft as ft
import scipy.integrate as integ
__sci=sci
__np=np
__py=py
__ft=ft
import Charm as chr
import cPickle as cpk
import time as tm
from Charm.Core import Misc as __misc
from Charm import _eps
#===================================================
# Smoothing real traces
#===================================================
D=sci.load('data/08finmig.pyd')
N=D.shape[0]
N-=N%2

rx=4
t=D[0:N,rx]
# t=t[2/0.004:4.5/0.004]

Ft=ft.fft(t) 
Ft[300:N+1-300]=0
tl=ft.ifft(Ft).real

print 'Trace #'+str(rx)+' of Dataset !'
I=t
s=3
tnorm=0
Nrm=__misc.norm(I)
N=__sci.size(I)
FI=__ft.fft(I)[0:N/2+1]        # fft=[0,+f,-f]

# w = 2*__sci.pi/N * __sci.arange(eps,N/2+1)
# OP=__sci.exp(-(w*s)**2/ 2)*FI
# F=__sci.concatenate(( OP[0:N/2],[OP[N/2]] ,__sci.conj(OP[N/2-1:0:-1] ) ))
# Is=ft.ifft( F)
# if tnorm != 1:
#     Is*=Nrm/__misc.norm(Is)
# ts=Is

ts=chr.smooth(t,s=s,fi=0,tnorm=0)[0]
tmav=chr.smooth(t,3)[0]

# py.figure()
# py.plot( tl,'y', ts,'r', t,'--b')
# py.legend(['Freq Fixed Data','Smooth Data','Data'])
# py.show()

py.figure()
py.subplot(211);py.plot( t,'b')
py.subplot(212);py.plot( tmav,'r')
py.show()


# D6=chr.char(type='trace',input=tl,user=0)
# D20=chr.char(type='trace',input=ts,user=1)

# for (name, value) in D20.items():
#      exec('%s = value' %name) 
# #===================================================
# # Integration
# #===================================================
# py.close('all')
# prefix='./data/'
# wtraces=sci.load(prefix+'wtraces.pyd')

# mymodel=open(prefix+'model.pyd','r+')
# model=cpk.load(mymodel)
# mymodel.close()
# trace=model['trace']
# events=model['events']
# attrib_org=model['attrib_org']


# N=1000
# ind=1
# fi=1
# scl=10

# f=chr.gmanf(N)
# f.update(attrib_org[ind,:])
# d=wtraces[ind,:]
# # d=f.data(tnorm=1)[0]+0.005*sci.randn(N)
# # d-=np.mean(d)
# d/=chr.norm(d)


# py.figure();
# f.plot(tnorm=1,tamp=0,c='b')
# py.plot(d,'r')

# g=f.copy(s=scl,fi=fi)
# f.alpha+=fi;f.s=scl
# ds=chr.smooth(d,s=scl,fi=fi,tnorm=1)[0]

# # dsum=np.cumsum(d)
# # dsum-=np.mean(dsum)
# # dsum/=chr.norm(dsum)
# # dsum=chr.phase_shift(dsum, phi=0.5,tnorm=1)

# py.figure();
# f.plot(tnorm=1,tamp=0,c='g')
# g.plot(tnorm=1,tamp=0,c='b',ls='--')
# py.plot(ds,'r')
# # py.plot(dsum,'--k')
# f.alpha-=fi

# py.show()




#===================================================
#  Manifolds data for alpha>0
#===================================================
# py.close()
# eps=1e-15

# N=1000
# tau=500
# sigma=18
# alpha=-5.8
# phi=1.30#-0.5*alpha

# ###################
# # Put DC Values of integrated waveform and manifold to zero
# ###################
# DC=0

# beta=1+alpha
# DCfit= -0.0223*(10*beta+1)**3 + 0.6651*(10*beta+1)**2 - 1.7708*(10*beta+1) + 2.5537

# df_NRM=1
# s=0
# fi=1
# #----------------
# sigma+=s
# alpha+=fi

# w=2*__sci.pi/N * __sci.arange(eps,N/2+1)
# F= (w)**(-alpha) *__sci.exp(1j*__sci.pi*phi) *__sci.exp(-(w*sigma)**2 /2)* __sci.exp(-1j*w* tau) 
# # if alpha > 0:
# F[0] = DC

# data0=chr.outtype(F,'time','complex',0,tnorm=df_NRM,stype='real')[0]

# data=np.cumsum(data0)
# data-=np.mean(data)
# data/=chr.norm(data)

# #=================================================
# alpha+=1;phi-=0.5
# w=2*__sci.pi/N * __sci.arange(eps,N/2+1)
# G= (w)**(-alpha) *__sci.exp(1j*__sci.pi*phi) *__sci.exp(-(w*sigma)**2 /2)* __sci.exp(-1j*w* tau) 
# # if alpha > 0:
# G[0] = DC

# gdata=chr.outtype(G,'time','complex',0,tnorm=df_NRM,stype='real')[0]
# print chr.norm(gdata)

# py.figure()
# py.plot(data-gdata)
# # py.plot(gdata,'r')


#=====================================================
#  Check Integration over Manifolds
#=====================================================
# a=-.8
# N=1000
# s=10
# df=np.zeros((10,N))
# dg=np.zeros((10,N))
# f=chr.gmanf(N,N/2.,s,a,0)
# g=chr.gmanf(N,N/2.,s,a,-a/2.)
# src=chr.gauswavf(N,s,p=2,domain='freq')[0]

# df[0,:]=ft.ifft(f.data(tnorm=1,domain='freq')[0]*src)
# df[0]/=chr.norm(df[0,:])
# dg[0,:]=ft.ifft(g.data(tnorm=1,domain='freq')[0]*src)
# dg[0]/=chr.norm(dg[0,:])

# for i in range(8):
#     df[i+1,:]=np.cumsum(df[i,:])
#     df[i+1]/=chr.norm(df[i+1,:])

# # py.figure()
# # py.plot(df2)
# # py.title('Zero-Phase')


# for i in range(8):
#     dg[i+1,:]=np.cumsum(dg[i,:])
#     dg[i+1]/=chr.norm(dg[i+1,:])
# # py.figure()
# # py.plot(dg4)
# # py.title('Phase = -alpha/2')

