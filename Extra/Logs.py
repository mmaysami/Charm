#! /usr/bin/env python 
__doc__="""

          Well Logs
    =====================
    Read and handle well log data
    
    read_las    -- Load well log data from LAS Ascii file
    read_pyd    -- Load well log data from python data file

    depth2time  -- Calculate logs versus 2-way travel time
    impedance   -- Calculate impedance vs. time or depth
    layers      -- Find 1-D acoustice profile 


    Note:    ____  Available Logs ____

                          dz=0.5 ft=0.1524 m
    TVDSS.(FT)            True Vertical  Depth (Sub-Sea) 
    GR.(GAPI)             Gravity
    HCAL.(IN)             Calliper (Well Radious)
    HLLD.(OHM)            Resistivity     
    DTCO.(USec/Ft)        Sonic (DTCO) 
    DTSM.(USec/Ft)        Sonic 
    DT1.(USec/Ft)         Sonic 
    RHOZ.(Gr/Cm3)	  Density (DRHO) 
    NPHI. 		  Neutron Porosity


AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""
# __all__=['read_las','read_pyd','depth2time','impedance','layers']


import os as __os
import numpy as __np
import cPickle as __cpk
import scipy.interpolate as __intp
import Charm.Misc as __misc,Charm.API as __API
from   Charm  import __path__, __Disp_Err, _eps
from   numpy    import ma  as __ma,fft as __ft
from   scipy.io import read_array as __read_array
from   os.path  import join as __join

try:
    import pylab as __py
except ImportError:
    __py=None


# Time-Depth Matching Info
_event_marks={
    'Seabed'                     :0,
    'Diaenetic Event'            :1,
    'Strachan Fan (top Breydon)' :2,
    'Intra Fan DHI (GWC)'        :3,
    'Base Strachan Fan'          :4,
    'Top High Amplitude Event'   :5,
    'Balder Tuff'                :6,
    'Top Lava'                   :7,
    'Base Lava'                  :8,
    'Total Depth'                :9}
_tvd_marks   = [
    1620.6216,
    2189.9880,
    2612.7456,
    2647.1880,
    2790.4440,
    2983.6872,
    3364.0776,
    3791.1024,
    3830.1168,
    4110.5328]  # [m]
_twt_marks   = [
    2.230,
    2.861,
    3.396,
    3.437,
    3.565,
    3.741,
    4.060,
    4.317,
    4.363,
    4.472]  # [sec]
# Dictionary of logs, their index in matrix, and conversion const. to SI unit
_metric = {
    'TVDSS':[0,0.3048],
    'DTCO' :[4,1e-6/0.3048],
    'RHOZ' :[7,1e3],
    'GR'   :[1,1],
    'HCAL' :[2,0.0254],
    'HLLD' :[3,1],
    'DTSM' :[5,1e-6/0.3048],
    'DT1'  :[6,1e-6/0.3048],
    'NPHI' :[8,1] }


# Finding path to data files
_log_null=-999.25
_data_path=__path__[2]     # Abs. path to data
_df_las="logs21441_Nohdr.las"
_df_file=__misc.abs_file(_df_las, _data_path)
_df_pyd="logs21441.pyd"
_df_pfile=__misc.abs_file(_df_pyd, _data_path)
_df_seismic='veritas08.rsf'
_df_sfile=__misc.abs_file(_df_seismic, _data_path)



#==========================================
#                Functions
#==========================================
def read_las(las=_df_las,SI=False,null=_log_null):
    """
    Read well Logs from Log Ascii Standard files and save it 
    as a python dictionary in a PYD file (with metric/SI units)

    Input:
         las  : Path/Name of LAS file of logs
         SI    : SI units flag of source LAS file logs 
                    (if False converts to metric)
         null : Null Value of measured data       
    Output:
         log_dict : Dictionary of logs & null value         
    """
    if not __os.path.isabs(las):
        las=__misc.abs_file(las, _data_path)
    data=__read_array(las)
    log_dict={'null':null}
    for (name,value) in _metric.items():
        ind,const  = value
        mask= __np.ones(data[:,ind].shape)
        mask +=  (const-1)* (SI==False)* (data[:,ind]!=null)
        log_dict.update({name.lower():mask*data[:,ind]})
    # Save Logs to PYD File
    mylog=open(__join(_data_path,_df_pyd),'w')
    __cpk.dump(log_dict,mylog)
    mylog.close()
    return log_dict



def read_pyd(pyd=_df_pyd):
    """
    Read metric well Logs from PYD(cPickles data) file to a dictionary

    Input:
         pyd    : Path/Name of PYD file of logs
    Output:
         log_dict : Dictionary of logs & null value

    """
    if not __os.path.isabs(pyd):
        pyd=__misc.abs_file(pyd, _data_path)

    mylog=open(pyd,'r')
    log_dict=__cpk.load(mylog)
    mylog.close()
    return log_dict


def interpolate(log=_df_pyd,keys=['tvdss','dtco','rhoz']):
    """
    Read metric well Logs from PYD(cPickles data) file to a dictionary

    Input:
         log    : Path/Name of PYD file of logs
         keys   : The logs type to be read from file and written in 
                  output dictionary(if None means all keys)
    Output:
         log_intp_dict : Dictionary of interpolated 
                         logs(based on input argument) & null value

    """
    log_dict=read_pyd(log)
    if keys==None:
        keys=log_dict.keys()
        keys.remove('null')

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
    return log_intp_dict


def depth2time(log=_df_pyd,mark_ind=2):
    """
    Convert logging depths to two-way travel time

    =====  Events mark  =====
    'Seabed'                    :0
    'Diaenetic Event'           :1
    'Strachan Fan (top Breydon)':2
    'Intra Fan DHI (GWC)'       :3
    'Base Strachan Fan'         :4
    'Top High Amplitude Event'  :5
    'Balder Tuff'               :6
    'Top Lava'                  :7
    'Base Lava'                 :8
    'Total Depth'               :9


    Input:
         log      : Path/Name of PYD file of logs
         mark_ind : Index from events marks of match point of depth and time
    Output:
         log_intp_dict : Updated dictionary ofinterpolated logs
                    with equivalent two-way traveltimes for depths
    """
    log_intp_dict=interpolate(log)
    # Read log values from dictionary
    vars_dict=['tvdss','dtco']
    for name in vars_dict:
        exec('%s = log_intp_dict[name]' %name)  


    # Form two-way travel time  
    twt0=2*__np.concatenate([ [0],__np.diff(tvdss)*dtco[0:-1] ])
    twt = __np.cumsum(twt0)  #[sec]
    ind = mark_ind    #_event_marks['Strachan Fan (top Breydon)']
    ind2= __np.argmin(abs (tvdss - _tvd_marks[ind]) )
    twt += (_twt_marks[ind] - twt[ind2])
    log_intp_dict.update({'twt':twt})

    return log_intp_dict

def impedance(log=_df_pyd,domain='time'):
    """
    Calculate acoustic impedance [kg/m^2.s]  versus time or depth from
    Well logs.

    Input:
         log    : Path/Name of PYD file of logs
         domain : indicate calcuting versus ['time']/'depth' 
    Output:
         imp    : Imppedance with respect to time/deph 
         t      : Time/Depth based on domain argument
    """
    log_intp_dict=depth2time(log=_df_pyd)
    for (name,value) in log_intp_dict.items():
        exec('%s= value' %name)  # includes twt

    imp_z = (rhoz/dtco).transpose()   # [kg/m2.s]
    if domain.lower()=='depth':
        return imp_z,tvdss

    # Re-define time samples
    t   = __np.linspace( min(twt), max(twt), len(twt))
    op=__intp.interp1d(twt,imp_z,bounds_error=True,fill_value=_log_null)
    imp=op(t).transpose()
    return imp,t


def ref_coeff(z,vp,den,const=1):
    """
    Computes the reflection coefficients

    Input:
       z    : Depth
       vp   : P-wave velocity profile along z
       den  : Density profile along z
       const: Constant magnifying factor for coefficients

    Output:
       R    : Reflection coefficients 
    """

    N=len(vp)
    dz=__np.diff(z)
    R  = __np.zeros(N)
    T  = __np.zeros(N)
    AI=vp*den

    for i in range(0,N-1):
        R[i]  = (AI[i+1]-AI[i]) / ((AI[i+1]+AI[i])*dz[i])
        T[i]  = 2*AI[i+1] / (AI[i+1]+AI[i])
    R[-1]=R[-2]
    T[-1]=T[-2]
    return R*const,T



def layers(z,vp,den,Nt=2001,dt=0.004,const=1):
    """
    Computes the flux normalized response of an 1D acoustic medium

    Input:
       z    : Depth
       vp   : P-wave velocity profile along z
       den  : Density profile along z
       Nt   : Number of Time samples 
       dt   : Time sampling rate (If None, uses whole range of well logs)

    Output:
       Rd   : Downgoing reflection coefficients in the time domain (complex)  
       Ru   : Upgoing reflection coefficients in the time domain (complex)  
       T    : Downgoing transmission coefficients in the time domain (complex)
    """

    dz=__np.diff(z)
    N=len(vp)
    if Nt == None:
        Nt = (N+N%2)    # Make Nt even
    else:
        Nt = (Nt+Nt%2) 
   

    # Frequencies
    Nf = (0.5*Nt)+1
    w  = __np.arange(0,Nf) * (2*__np.pi/(Nt*dt))
    w  = (w+1j*_eps).reshape(Nf,1)
    # Initialise the GLOBAL quantities
    Rd  = __np.zeros((Nf,1),dtype='complex')
    Ru  = __np.zeros((Nf,1),dtype='complex')
    T   = __np.ones((Nf,1),dtype='complex')

    # Start the recursion loop over the N-1 layers
    for i in range(0,N-1):
        # Vertical slowness
        s1 = (vp[i]**-1)
        s2 = (vp[i+1]**-1)
        # Reflection coefficients
        r    = (den[i+1]*s1 - den[i]*s2) / (den[i+1]*s1+den[i] *s2)
        r    = __np.ones((Nf,1)) * r * const

        # Calculate the phase shift operator
        s   = __np.ones((Nf,1))*s1
        phi  = __np.exp(1j * s * w * dz[i])
        # Calculate the R downgoing & Upgoing + T downgoing
        Rd  += (T**2)*(phi**2) *r
        Ru   = -r + (phi**2)*Ru
        T   *= phi

    # Calculate the inverse fft's and apply the taper (laplace)
    T=__ft.irfft(__np.ravel(T).conj())
    Rd=__ft.irfft(__np.ravel(Rd).conj())
    Ru=__ft.irfft(__np.ravel(Ru).conj())

    #TT1=__misc.fillfreq(T,tnorm=0,stype='real')
    #TT2=    __np.real(__ft.ifft(__np.conj(TT1)))
    #RR1=__misc.fillfreq(Rd,tnorm=0,stype='real')
    #RR2=    __np.real(__ft.ifft(__np.conj(RR1)))
    
    return Rd,Ru,T





def layers_multiple(z,vp,den,Nt=2001,dt=0.004):
    """
    Computes the flux normalized response of an 1D acoustic medium

    Input:
       z    : Depth
       vp   : P-wave velocity profile along z
       den  : Density profile along z
       Nt   : Number of Time samples 
       dt   : Time sampling rate (If None, uses whole range of well logs)


    Output:
       Rd  : Downgoing reflection coefficients in the time domain (complex)  
       Ru  : Upgoing reflection coefficients in the time domain (complex)  
       T   : Downgoing transmission coefficients in the time domain (complex)
    """
    # Setting for p Values
    Np,dp,p0,ep = 1,0,0,_eps

    dz=__np.diff(z)
    N =  len(vp)
    if Nt == None:
        Nt = (N+N%2)    # Make Nt even
    else:
        Nt = (Nt+Nt%2) 

    # Frequencies
    Nf = (0.5*Nt)+1
    w  = __np.arange(0,Nf) * (2*__np.pi/(Nt*dt))
    w  = (w+1j*ep).reshape(Nf,1)

    # Horizontal slowness
    p   = __np.arange(0,Np)*dp + p0
    # Initialise the GLOBAL quantities
    Rd  = __np.zeros((Nf,Np),dtype='complex')
    Ru  = __np.zeros((Nf,Np),dtype='complex')
    T   = __np.ones((Nf,Np),dtype='complex')

    # Start the recursion loop over the N-1 layers
    for i in range(0,N-1):
        # Vertical slowness
        s1 = __np.sqrt(vp[i]**-2   - p**2)
        s2 = __np.sqrt(vp[i+1]**-2 - p**2)
        # Reflection coefficients
        r    = (den[i+1]*s1 - den[i]*s2) / (den[i+1]*s1 + 
                                              den[i] *s2)

        # Calculate the phase shift operator ,s [Nf*Np]
        s    = __np.ones((Nf,1))*s1
#         s    = __np.real(s) + 1j*(__np.sign(__np.real(w))*__np.ones((1,Np)))*__np.imag(s)
        s    = __np.dot( __np.real(s)+1j*__np.sign(__np.real(w))*__np.imag(s)
                         ,__np.ones((1,Np)))

        r    = __np.ones((Nf,1)) * r
        wvec = w * __np.ones((1,Np))
        phi  = __np.exp(1j * s * wvec * dz[i])

        # Calculate the R downgoing & Upgoing + T downgoing
        Rd  += (T**2)*(phi**2) *r
        Ru   = -r + (phi**2)*Ru
        T   *= phi

    # Calculate the inverse fft's and apply the taper (laplace)
    T=__ft.irfft(__np.ravel(T).conj())
    Rd=__ft.irfft(__np.ravel(Rd).conj())
    Ru=__ft.irfft(__np.ravel(Ru).conj())
    # TT1=__misc.fillfreq(T,tnorm=0,stype='real')
    # TT2=    __np.real(__ft.ifft(__np.conj(TT1)))
    # RR1=__misc.fillfreq(Rd,tnorm=0,stype='real')
    # RR2=    __np.real(__ft.ifft(__np.conj(RR1)))
    
    import cPickle as __cpk
    myfile=open('./datalayer.pyd','w+')
    __cpk.dump(dict(Rd=Rd,Ru=Ru,T=T,Nt=Nt,dt=dt),myfile)
    myfile.close()
    return Rd,Ru,T



# def match(seismic=_df_sfile,tracenum=870,log=_df_pyd):
#     """
#     (IN PROGRESS)
#     Input:
#        log : Path/Name of PYD file of logs

#     Output:

#     """

#     log_dict=depth2time(log=_df_pyd)
#     trace,hdr=__API.sftrace2d(infile=seismic,outfile=None,tracenum=tracenum)
#     twt= log_dict['twt']
#     dt=hdr['d1']
#     Nt=(twt[-1]-twt[0]) / dt

    

#     log_intp_dict=depth2time(log=log)
#     for (name,value) in log_intp_dict.items():
#         exec('%s= value' %name)
#     if dt == None:
#         dt =  (twt[-1]-twt[0])/Nt
#     Rd,Ru,T = layers(tvdss,1./dtco,rhoz,Nt=Nt,dt=dt)

#     __py.figure();__py.plot(Rd)
#     __py.figure();__py.plot(trace)
#     pass











# ************************************************************
#                Main Body (Running Script)
# ************************************************************
# if __name__ == '__main__':

#     tsize=__py.rcParams['axes.titlesize']
#     lsize=__py.rcParams['axes.labelsize']
#     # Compute Acoustic Impedance & Reflectivity
#     import Charm.Manifold as __manif
#     imp,t,dt = impedance(log=_df_pyd,domain='time')
#     r=__np.diff(imp)
#     rs=__manif.smooth(r,s=15)[0]
#     print len(rs)
#     if len(rs)%2 != 0:
#         rs=__np.concatenate([rs,[0]])
#     source=__manif.gmanf(len(rs),0,8,-2,0)
#     s=__ft.ifft(__ft.fft(rs)*source.data(domain='freq',tnorm=1)[0])

#     #================================
#     #         Plot Figures 
#     #================================
#     if __py=None:
#         raise ImportError,__Disp_Err
#     # Set path for saving figures
#     prefix=__join('demos','Seismic')
#     if not __os.access(prefix,__os.F_OK):
#         __os.makedirs(prefix)

#     saveflg=0
#     """
#     __py.figure()
#     __py.plot(tvdss,dtco_ma)
#     __py.title('Well Logs',size=tsize)
#     __py.xlabel('Depth[ft],  TVDSS',size=lsize)
#     __py.ylabel('Sonic Log [us/ft],   DTCO',size=lsize)
#     __py.show()
#     """

#     __py.figure()
#     __py.plot(t[0:len(rs)],rs)
#     __py.title('Smooth Reflectivity',size=tsize)
#     __py.xlabel('Time [sec]',size=lsize)
#     __py.ylabel('Smooth Reflectivity',size=lsize)
#     __py.show()

#     if saveflg:
#         __py.savefig(__join(prefix,"logs_%s.eps" %f))
#     saveflg=0
