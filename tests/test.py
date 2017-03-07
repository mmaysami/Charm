#! /usr/bin/env python 
import scipy as __sci
import scipy as sci
import numpy as __np
import numpy as np

import pylab as py

import numpy.fft as __ft
import numpy.fft as ft

import Charm as chr
# import StbChr as stbchr
import cPickle as cpk
import time as tm

#============================
import os as __os
import numpy as __np
import cPickle as __cpk
import scipy.interpolate as __intp
from   Charm.Core import Misc as __misc,API as __API
from   Charm    import __path__, __Disp_Err, _eps
from   numpy    import ma  as __ma,fft as __ft
from   scipy.io import read_array as __read_array
from   os.path  import join as __join
from Logs import *
#=========================
# Well log interpolate
#=========================
keys=['tvdss','dtco','rhoz']
log_dict=read_pyd()
# Read log values from dictionary
null=log_dict['null']
for name in keys:
    exec('%s = log_dict[name]' %name)  

rng1=__np.nonzero(dtco - null)[0]
rng2=__np.nonzero(rhoz - null)[0]
(lb,ub)=__np.maximum(rng1[0],rng2[0]),__np.minimum(rng1[-1],rng2[-1])

# Mask data when any of them is Null
mask=__np.logical_or(dtco==null,rhoz==null)
log_intp_dict={'null':log_dict['null']}
for name in keys:
    exec("tmp= __ma.masked_where(mask,%s)" %name)
    exec("%s_cmp= tmp.compressed()" %name)
    exec("%s = %s[lb:ub+1]" %(name,name))  
    exec("intp_ob = __intp.interp1d(tvdss_cmp,%s_cmp,bounds_error=True,fill_value=null)" %name)  
    exec("%s_i=intp_ob(tvdss)"%name)
    exec("log_intp_dict.update({name:%s_i})" %name)


#========================================
#  G. Manifolds 
#========================================
# py.figure()
# for ph in np.arange(0,1,0.1):
#     f=chr.gmanf(1000,512,10,-1,ph)
#     f.plot()

# py.show()

# f3=stbchr.gmanf(1000,500,10,-3.4)
# print f3.__param__()
# py.figure()
# rng=[0,3]
# c=0
# for k in rng:
#     c+=1
#     py.subplot(len(rng)*100+10+c)

#     if k==0:
#         py.plot(f3.data()[0],lw=1)
#         py.plot(np.imag(f3.data()[0]),'r',lw=1)
#     else:
#         py.plot( f3.eval(att='d'+str(k-1))()[0],lw=1)
#         py.plot(np.imag(f3.eval(att='d'+str(k-1))()[0]),'r',lw=1)

# py.show()


# f4=chr.gmanf(1000,500,10,-3.4)
# print f4.__param__()
# py.figure()
# rng=[0,3]
# c=0
# for k in rng:
#     c+=1
#     py.subplot(len(rng)*100+10+c)
#     if k==0:
#         py.plot(f4.data()[0],lw=1)
#         py.plot(np.imag(f4.data()[0]),'r',lw=1)
#     else:
#         py.plot( f4.eval(att='d'+str(k-1))()[0],lw=1)
#         py.plot(np.imag(f4.eval(att='d'+str(k-1))()[0]),'r',lw=1)


# py.show()

#========================================
#  Testing different type of Splines
#========================================
# N=1000
# py.figure(1)
# # r_spk,k,d,reg=chr.spiketrain(N,k=1,delta=100,amp=(2,10),regul=0,uniform=0)
# r_spk=np.zeros(N);r_spk[N/2.]=1
# fr=np.fft.fft(r_spk)

# alpha=-0.2
# phi=0#-0.5*(alpha+1)
# spf=chr.fspline(N,alpha=alpha,phi=phi,loc=0,tnorm=1,domain='freq',caus=1)

# d=fr*spf
# td=np.fft.ifft(d)

# alpha=alpha
# phi=-0.5+phi#-0.5*(alpha-1)
# spf=chr.fspline(N,alpha=alpha,phi=phi,loc=0,tnorm=1,domain='freq',caus=1)
# d=fr*spf
# td2=np.fft.ifft(d)
# py.clf()
# py.plot(td,'b',lw=1)
# py.plot(td2,'--r',lw=1)
# py.plot(r_spk/max(abs(r_spk)) *max(abs(td)),'g')
# py.show()

# #========================================
# #  Testing different type of Splines
# #========================================
# N=1200
# alpha=-0.34
# nbp=12
# stp=2./nbp
# # Phase is between any real numbers [a,a+2] + phase bias set to phaseadd
# phiadd= -(alpha+1)/2
# rng= np.arange(-1+3*stp,1+3*stp,stp) +phiadd
# # rng=np.concatenate(( np.arange(-1,1,stp),np.arange(-1,1,stp) ))+phiadd

# py.figure(2)
# py.clf()
# py.figure(3)
# # py.clf()

# scl=64

# for k in range(len(rng)):
#     phi=rng[k]
#     sp=chr.fspline(N,alpha=alpha,phi=3+phi,loc=2*N/5.,tnorm=1,domain='time')

#     spf=chr.fspline(N,alpha=alpha,phi=0+phi,loc=2*N/5.,tnorm=1,domain='freq')
#     src=chr.wavelet('gaus',wtype='real',param={'p':2})
#     srcf=-src.wavf(N,s=scl,domain='freq')[0]
#     evn=np.fft.ifft(spf*srcf)

#     py.figure(2)
#     py.plot(sp)
#     py.figure(3)
#     if k<nbp:
#         py.plot(evn)
#     else:
#         py.plot(evn,'--',lw=2)
# py.legend()
# py.show()


#==================================
# Test Hilbert Transform
#==================================
# py.cla()
# t=np.linspace(-11,24,1000)
# d= np.random.randn(100)
# # H=np.sin(t)


# Fd=__ft.fft(d)
# w = __np.arange(len(d)/2+1)
# w[-1]=0
# w = __np.concatenate( ([w,-w[-2:0:-1]]))
# Hd = - 1.0j * __np.sign(w)*Fd
# hd=__ft.ifft(-Hd).real

# hfftd=sci.fftpack.hilbert(d)

# # py.plot(H,'g',lw=1)
# py.plot(hd,'b',lw=1)
# py.plot(hfftd,'r',lw=1)
