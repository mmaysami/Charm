#!/usr/bin/env python

"""
      Synthetic seismic modeling
    =================================
=>  read_elog(elog=_df_elog,SI=False):
        Read elastic well Logs from files and save it 
        as a python dictionary in a PYD file (with metric/SI units)


=>  read_epyd(pyd=_df_pyd):
        Read metric well Logs from PYD(cPickles data) file to a dictionary


=>  signature(N=1000):
        Extract source signature from sea bottom (Approximation)


=>  endmember(pyd=_df_pyd,PC=0.3116,BETA=0.41):
        Generates dictionary of end-member properties for 
        Opal A-Opal CT binary mixture from input elastic logs (Approximation)


    Note:    ____  Available Logs ____=

    Depth (mbsf) : meter below sea floor
    Density(g/cc)
    Vp(m/s)	
    Vs(m/s)


    Note:     Opal Definition
         Microcrystalline Opal
         Opal-CT has been interpreted as consisting of clusters of stacking 
         of cristobalite and tridymite over very short length scales. 

         Non-crystalline Opal
         Two broad categories of non-crystalline opals, sometimes just referred
         to as opal-A (amorphous). Non-crystalline silica in siliceous sediments
         is reported to gradually transform to opal-CT as a result of diagenesis
         , due to the increasing overburden pressure in sedimentary rocks,

AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""

import rsf       as __rsf
import os        as __os
import Charm     as __chr
import numpy     as __np
import numpy.fft as __ft
import pylab     as __py
import cPickle   as __cpk
from   scipy.io   import read_array as __read_array
from   os.path    import join as __join
# Dictionary of logs, their index in matrix, and conversion const. to SI unit
_metric = {
    'TVDSS':[0,1],
    'RHO'  :[1,1e3],
    'VP'   :[2,1],
    'VS'   :[3,1]}
_data_path=__chr.__path__[2]     # Abs. path to data
_df_elog="logs_904A_Nohdr.txt"
_df_pyd="elogs_904A.pyd"
_df_file=__chr.abs_file(_df_elog, _data_path)



def read_elog(elog=_df_elog,SI=False):
    """
    Read elastic well Logs from files and save it 
    as a python dictionary in a PYD file (with metric/SI units)

    Input:
         elog     : Path/Name of log file 
         SI       : SI units flag of source file logs 
                    (if False converts to metric)
    Output:
         log_dict : Dictionary of logs        
    """
    data=__read_array(_df_file)
    elog_dict={}
    for (name,value) in _metric.items():
        ind,const  = value
        mask= __np.ones(data[:,ind].shape)
        mask +=  (const-1)* (SI==False)
        elog_dict.update({name.lower():(mask*data[:,ind])[::-1] })
    # Save Logs to PYD File
    mylog=open(__join(_data_path,_df_pyd),'w')
    __cpk.dump(elog_dict,mylog)
    mylog.close()
    return elog_dict



def read_epyd(pyd=_df_pyd):
    """
    Read metric well Logs from PYD(cPickles data) file to a dictionary

    Input:
         pyd    : Path/Name of PYD file of logs
    Output:
         log_dict : Dictionary of logs & null value
    """
    if not __os.path.isabs(pyd):
        pyd=__chr.abs_file(pyd, _data_path)
    mylog=open(pyd,'r')
    log_dict=__cpk.load(mylog)
    mylog.close()
    return log_dict


def endmember(pyd=_df_pyd,PC=0.3116,BETA=0.41):
    """
    Generates dictionary of end-member properties for 
    Opal A-Opal CT binary mixture from input elastic logs (Approximation)

    Input:
         pyd  : File that has elastic logs for Opal transition
         PC   : Critical volume fraction in HB percolation model
         BETA : Fractional order of switch in HB percolation model
    Output:
         dic : Dictionary of end-members of binary mixture
    """
    elog=read_epyd(pyd)
    for (name,value) in elog.items():
        exec('%s = value' %name)

    ind1=__np.where(__np.logical_and(tvdss>=515,tvdss<=523))[0]
    ind2=__np.where(__np.logical_and(tvdss>=527,tvdss<=535))[0]

    # Elastic properties of end members of binary mixture
    # Opal A  & Opal CT
    VP1=__np.average(vp[ind1])
    VP2=__np.average(vp[ind2])
    VS1=__np.average(vs[ind1])
    VS2=__np.average(vs[ind2])
    DEN1=__np.average(rho[ind1])
    DEN2=__np.average(rho[ind2])

    G1=DEN1*(VS1**2)
    G2=DEN2*(VS2**2)

    K1=DEN1*(VP1**2) - 4*G1/3.
    K2=DEN2*(VP2**2) - 4*G2/3.

    DENC=(1-PC)*DEN1+PC*DEN2

    dic=dict(
        A=1,
        CT=2,
        DEN1 =DEN1,
        DEN2 =DEN2,
        K1   =K1,
        K2   =K2,
        G1   =G1,
        G2   =G2,
        DENC =DENC,
        PC   =PC,
        BETA =BETA,
        LP   =1,
        HP   =2,
        VP1  =VP1,
        VP2  =VP2,
        VS1  =VS1,
        VS2  =VS2)

    """ # Show from which interval the averages are taken
     __py.plot(tvdss,vp)
     __py.plot(tvdss[ind1],vp[ind1])
     __py.plot(tvdss[ind2],vp[ind2])
     __py.show()"""
    return dic


def signature(N=1000,thr=1e-3,dt=4e-3):
    """
    Extract source signature from sea bottom (Approximation)

    Input:
         N   : Length of signature
         thr : Percentage of max freq. element to pick the bandwidth
         dt  : Time sampling in second
    Output:
         W     : Source signature as a Gaussian manifold
         amp : Max amplitude of signature
                   L2norm(events_det)*max(abs(gManif(tnorm=1) ))
         bw   : Bandwidth of signature in Hertz
    """
    infile='xSeaBtm01.rsf'
    data,hdr0=__chr.sfread3d(infile)

    # Well= 830 - 840
    d=data[:,830:840,40:45]
    amp = d[2,:,:].transpose()
    sigma = d[3,:,:].transpose()
    alpha = d[4,:,:].transpose()
    phi   = d[5,:,:].transpose()
    amp_ma=__np.ma.masked_equal(amp,-100)
    sigma_ma=__np.ma.masked_equal(sigma,-100)
    alpha_ma=__np.ma.masked_equal(alpha,-100)
    phi_ma=__np.ma.masked_equal(phi,-100)
    
    amp=__np.array(__np.ma.average(amp_ma))
    scl=__np.array(__np.ma.average(sigma_ma))
    alp=__np.array(__np.ma.average(alpha_ma))
    phi=__np.array(__np.ma.average(phi_ma))

    W=__chr.gmanf(N,0,scl,alp,phi)
    FW=abs( W.data(domain='freq',mode='complex')[0][0:0.5*N])
    bw=len(__np.nonzero(FW>thr*max(FW))[0]) /(N*dt)
    
    return W,amp,bw


def line(p0,p1,x):
    """
    XXXX
    """
    m=  (p1[1]-p0[1]) * 1.0 / (p1[0]-p0[0]) 
    b=p0[1]-m*p0[0]

    y=m*x+b
    return y




# ************************************************************
#                Main Body (Running Script)
# ************************************************************
if __name__ == "__main__":

    
    PC   = 0.43 
    BETA = 0.81

    well=range(810,860)
    dic= endmember(_df_pyd,PC,BETA)
    for (name,value) in dic.items():
        exec(name+"="+str(value))

    data_file=['xBFGSW01.rsf', 'xLMW01.rsf','xLMW21.rsf'][2]
    data,hdr=__chr.sfread3d(data_file)
    opal_st=(2.87-hdr['o1'])/hdr['d1'] 
    section=data[0,:,:].transpose()

    """
    End points of 'Diaenetic Event'
      Lambda=V/f=2000/50=40m
      Z_int=0.5*100ms *2000m/s=100m
      dz=0.1524       # 0.5ft
      z1,z2=2158 # ??? ,2190
    """
    dt=hdr['d1']   # 0.004
    t1,t2,t3 = 2.79166, 2.88333, 2.925
    p1,p2,p3 = 0.1, 0.76, 0.3

    Nt=section.shape[0]
    #-------------------------
    #    Using Layer Code
    #-------------------------

    zb,ze=0,1
    zo1,zo2=0,1   #300,420
    dz= (ze-zb)/1000.    #0.1524 

    z=__np.arange(zb,ze,dz)
    zz=dz*__np.ones(z.shape)
    ind1=__np.argmin(abs(z-zo1))
    ind2=__np.argmin(abs(z-zo2))
    zl=z[0:ind1]
    z0=z[ind1:ind2]
    zr=z[ind2:]
    

    # define p,Vp,Den of binary mixture
    ps0=p0=__np.linspace(p1,p2,len(z0))
    vp0 =__chr.vp_model(ps0,'HB',dic)
    den0=__chr.vp_model(ps0,'DEN',dic)
    (VP_L,VP_U,VS_L,VS_U,RHO)=__chr.bounds(ps0,mode='RV',dic=dic)

    # Extend from both sides
    den=line([z0[0],den0[0]],[z0[-1],den0[-1]],z)
    vpl=line([z0[0],vp0[0]],[z0[2],vp0[2]],zl)
    vpr=line([z0[-2],vp0[-2]],[z0[-1],vp0[-1]],zr)
    vp=__np.concatenate([vpl,vp0,vpr])


    # Compute Reflection Coefficients
    Rd,Ru,T=__chr.layers(z,vp,den,Nt=Nt,dt=1e0*dt,const=1) 
    RC,TC= __chr.ref_coeff(z,vp,den,const=ze-zb)
    R=RC 

    # Compute Seismic Signature nad Signal 
    W , amp , bw = signature(N=len(R))
    event_norm=amp/ max(abs(W.data(tnorm=1)[0]))
    FW=W.data(domain='freq',mode='complex',tnorm=1)[0]
    signal0=event_norm*__np.real(__ft.ifft(FW*__ft.fft(R)))
    W.n=Nt

    # Crop edge effects and apply required shift and boosting
    boost=15 
    cent=__np.argmax(abs(signal0))
    signal,mask=__chr.wboxcar(signal0,cent,15)
    if len(signal) > Nt:
        signal=signal[cent-Nt/2.:cent+Nt/2.]

    trace=- __chr.vector(signal*boost,'col')
    seismic=trace*__np.ones((len(trace),len(well)))

    
    print " "
    print "BW of Signature:      ",bw
    print " "
    print "SUM(dz/Vp):           ",__np.sum(dz/vp)
    print " "
    print "SUM(dz/Vp[0:switch]): ",__np.sum(dz/vp[0:ind1])
    print " "
    print "END    DIFFS:         ",(z[-1]-z[0])/(vp[-1]-vp[0])
    print "SWITCH DIFFS:         ",(z[ind2]-z[ind1])/(vp[ind2]-vp[ind1])
    print " "

    CHR_dict=__chr.char('trace',trace,user=0,log=open('mlog.txt','w+'))
    trace_est=CHR_dict['trace_est']
    attrib_vect_est=CHR_dict['attrib_vect_est']
    
       
    #============================
    #     Plotting Results
    #============================
    _im_dict ={'interpolation':'nearest','aspect':'auto'}
    fn1='elastic.pdf'
    fn2='FigMode2.pdf'
    fn3='synthetic.pdf'
    fn4='wellseismic.pdf'
    fn5='welltie.pdf'

    __py.figure();
    __py.subplot(411);__py.plot(z0,ps0,'g')
    __py.ylabel("Volume Fraction");__py.title ("Transition model")
    __py.subplot(412);__py.plot(p0,den0,'--r');__py.ylabel("Density")
    __py.xlabel("Volume fraction (p)")
    __py.subplot(413);__py.plot(p0,vp0);__py.ylabel("Vp")
    __py.plot(p0,VP_L,'--g'),__py.plot(p0,VP_U,'--r')
    __py.subplot(414);__py.plot(p0[0:-1],__np.diff(vp0));__py.ylabel("Diff Vp")
    __py.savefig(fn1)


    __py.figure()
    __py.subplot(311);__py.plot(z,vp);__py.ylabel('Vp')
    __py.subplot(312);__py.plot(z[0:-1],__np.diff(vp));__py.ylabel("Diff Vp")
    __py.subplot(313);__py.plot(__np.arange(len(R))*dt,R);__py.ylabel('Rd')  
    __py.savefig(fn2)


    __py.figure();
    __py.subplot(131);W.plot(tshift=1,label=0,tnorm=0,tamp=amp,c='r')
    __py.title ("Source Signature")
    __py.subplot(132);__py.plot(trace,'k')
    __py.title ("Seismic Trace")
    __py.subplot(133);__py.imshow(seismic,cmap=__py.cm.gray,**_im_dict)
    __py.title ("Seismic slice");__py.colorbar()
    __py.xlabel('Boost: '+str(boost))
    __py.savefig(fn3)


    data[:,well,:]=0
    data[0,well,:]=trace.transpose()
    data[1,well,:]=trace_est.transpose()
    data[2,well,:]=attrib_vect_est[0,:]
    data[3,well,:]=attrib_vect_est[1,:]
    data[4,well,:]=attrib_vect_est[2,:]
    data[5,well,:]=attrib_vect_est[3,:]
    __py.figure()
    __chr.showgray(data[0,:,:].transpose(),label=True,tags=hdr)
    __py.title ("Seismic Section");__py.colorbar()
    __py.savefig(fn4)


    __chr.sfwrite(data.transpose(),'welltie.rsf',hdr)
    __chr.show(data,im_dict=_im_dict,showflg=0,tags=hdr)
    __py.savefig(fn5)


    att=data[-2,620:690,opal_st:opal_st+13].transpose()
    MSK_att=__chr.showmasked(att,limit=(-4,0),out=1)
    print  "Ref. Coeff. ",(den[-1]*vp[-1]-den[0]*vp[0])/(
        den[-1]*vp[-1]+den[0]*vp[0])
    print  "Source Order:       ",W.alpha
    print  "Real Seismic Order: ",__np.ma.average(MSK_att)
    print  "Synthetic Order:    ",CHR_dict['attrib_est'][:,2]


# import scipy as __sci
# __sci.io.savemat('mdata.mat',dict(vp=vp,z=z,den=den))
