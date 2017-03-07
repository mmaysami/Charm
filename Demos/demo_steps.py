#!/usr/bin/env python
"""
AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""
import pylab as py
import Charm as chr
tsize=py.rcParams['axes.titlesize']
lsize=py.rcParams['axes.labelsize']



#====================================
#     Generate a synthetic trace  
#====================================
N=1000
k=11
delta=60
trace,events,attrib_org,attrib_vect,k,delta,regul,model_dict=chr.rand_trace(N,k,delta)

py.figure()
py.plot(trace)
py.title('Synthetic trace')
#====================================
#   Detection & Segmentation steps  
#====================================
attrib_det,mml,cw,lnpos,ext,max_points=chr.detect(trace)
events_det,masks=chr.segment(trace,attrib_det)


chr.mmlimage(max_points,mml,abs(cw))
#====================================
#   Estimate single events with BFGS 
#====================================
ind=0         # index of event to be characterized
alpha= -1     # guess for alpha
guess=attrib_det[0,:]
guess[2]=alpha

est_dict={'dphi':-1j,'flatphi':1,'tol':1e-7,'verb':1,'show':1}
est_dict.update({'prm_ind':[1,2],'niter':500,'smth':0,'intg':0})
chr.BFGSestimate(events_det, guess, ind=0,**est_dict)



