#! /usr/bin/env python 
"""

       Main (reflector Characterization analysis)
    ================================================
    Main Body of the pakcage that deals with characterization analysis

    char       --  Wrapped Characterization Function
    generate   --  Serves 'char' function with providing input data



AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""

__all__=['char','generate']


from sys     import stderr as __stderr
from os      import path   as __path
from Charm   import _Null, _df_input,_df_rx,_df_xlabel
import Charm      as __chr
import Charm.Core as __core
import numpy      as __np
import cPickle    as __cpk
import time       as __tm
import os         as __os

try:
    import pylab as __py
    __py.close('all')
except ImportError:
    __py = None


# Default answers for user prompts
_df_est   = 'yes'
_df_sav   = 'no'
_df_fig   = 'no'

# Finding path to data files
data_path =__chr.__path__[2]        # Abs. path to data
pydata    =__chr.__path__[4]        # Abs. path to pydata


# ************************************************************
#                        Characterization Function
# ************************************************************
def char(type='rsf',input=_df_input,param=None,user=0,
         solver='BFGS',smooth=0,major=0.1,log=__stderr):

    """
    Characterize a single seismic trace

    Generate / Load seismic data(trace) for analysis
    
    Input:
    
        type  : Type of input data ['rsf']/'trace'/'new'/'loadreal'/'loadsynt'
                rsf      : get data from rsf files
                trace    : use givenarray/matrix by input argument as trace
                new      : Generate a new synthetic trace
                loadreal : Reload stored real model from last run
                loadsynt : Reload stored real model from last run 
                
        input : RSF Data file or Seismic trace (according to type argument)

        This may be required based on type argument :
        param : Trace number for rsf files and dictionary of synthesize 
                parameters for type='new' as below)

                N    : # Length of trace                
                K    : # Number of events
                Delta: # Minimum distance between events
                amp  : # Amp Ranges 
                rnga : # Alpha Ranges
                rngp : # Phase Range
                rngs : # Scale Range for gaussian trace generator
                scl  : # Scale for synthesis with splines and source wavelet
                regul: # Regularized spacing between events
                uniform: # Uniform amplitude for all events

                
        user  : Promp user for performing estimation and other steps if one
        solver: Parametric inversion method in estimation ['BFGS']/'LM'/'LS'
        smooth: Smooth factor for input trace
        major : Threshold level (pick up events with amp. percentage above this)
        log   : File name to print logs
    Output:
        Characterization results
    """

    # Set answers for no user case
    ans_est=_df_est
    ans_sav=_df_sav
    ans_fig=_df_fig

    #==================================================================
    #             ============== Initializing  ===============
    #==================================================================
    model_dict,wave_dict,attrib_org,attrib_vect,tm_syn=generate(
        type=type,input=input,param=param,log=log)
    for (name, value) in model_dict.items():
        exec('%s = value' %name)    # Includes Abstract , N

    trace=__core.Manifold.smooth(trace,s=smooth,fi=0,tnorm=0)[0]
    #==================================================================
    #          ==================<< Parameters >>================
    #==================================================================
    #---------- Multi-scale detection by MML --------
    "#     scl_rng=(-2,2)  # Scales for Morlet WT"
    scl_rng=(-0.5,3.5) #(-1,3.5)
    if not abstract:
        scl_rng=(0,4)
    scale_match=1.35

    """
    scl        : # scale range for wavelet transform 1(1.5),4
    fineloc    : # Put fine-scale location of MML as event location
    scale_match: # Match the scale of CWT to real scale   
    major      : # Pick events with CWT coefficient down to this percent of 
                   global maximum of CWT coefficients (from input)

    con     : # Connectivity Factor for MML
    rad     : # topological radius for extrema>=2 is  advised
    acc     : # Required accuracy for Extrema 
    search  : # Portion of scales to be searched for start point of MMLs
    negline : # Include Negative MMLLs
    mindist : # Min. Distance of MMLs to be detected as separate lines
    chkdist : # How far (# of offsets) from MML location, look for local maxima
    """
    cwt_dict={'con':__np.sqrt(12),'rad':2,'acc':1e-6,'search':0.25,
              'negline':0,'mindist':10,'chkdist':3,'fineloc':0,'log':log,  
              'scales':2**__np.linspace(scl_rng[0],scl_rng[1],60)}


    #---------- Estimation Parameters --------
    """
    estimate_type  : # Type of optimization for esimation([BFGS]/ LM /LS/LSQR)
                       Use BFGS or LM for best results
    _df_Ai         : # Default initial guess for alpha
    est1_dict      : # Dictionay of settings for first estimtion
    est2_dict      : # Dictionay of settings for second estimtion

    prm_ind: # index of parameters to inverse for 
    niter  : # Number of iteration
    smth   : # Smoothing factor for events
    intg   : # Fractional order integration for events
    dphi   : # Constant value for derivative of manifolds w.r.t. phase
    flatphi: # Make phase between 0,2*pi
    tol    : # Error tolerance to stop estimation
    verb   : # Verbosity of event estimation
    show   : # Show graphical results of event estimation
    """
    _df_Ai= -1     # -1 for LM & -0.2 for LS  
    estimate_type =solver  # "BFGS","LM":Levenberg-Marquardt,"LS":Line-Search
    est1_dict={'dphi':-1j,'flatphi':1,'tol':1e-7,'verb':0,'show':0}
    est2_dict=est1_dict.copy()
    est1_dict.update({'prm_ind':[1,2],'niter':500,'smth':0,'intg':0,'log':log})
    est2_dict.update({'prm_ind':[1,2,3],'niter':50,'smth':0,'intg':0,'log':log})



    #==================================================================
    #  ====== Detection & Numerical Results + Segmentation ======
    #==================================================================
    tm_det=__tm.time()
    print >>log,'\n   =============== Detection & Segmenting ==============='
    wr=__core.Cwt.wavelet('gaus',wtype='complex',param={'p':2})
    # wr=__core.Cwt.wavelet('morl',wtype='complex',param={'w0':1.5*__np.pi})
    (attrib_det,mml,cw,lnpos,ext,max_points)=__core.Steps.detect(
    trace,wr=wr,scl=cwt_dict['scales'],major=major,**cwt_dict)
    print >>log, "major=",major
    acw=abs(cw)
    wave_dict.update({'wr':wr,'cw':cw})
    Kdet=attrib_det.shape[0]

    """# ----------------- Segmenting ------------------
    # The window setting based on raw attributes 
    (Do segmenting first and then update attributes)"""
    events_det,masks=__core.Steps.segment(trace,attrib_det)
    trace_det=__np.sum(events_det,axis=0)

    if attrib_det.shape[0] !=0:
        # Consider one order difference in phase values (originial & Detected)
        attrib_det[:,1]*=scale_match
        "if wavelet is not zero-phase subtract its phase e.g. -1 for Ricker"
        attrib_det[:,3]=(attrib_det[:,3])%2


        "  ======== Results ========="
        if not abstract:
            print >>log, 'Original Events :',attrib_org[:,0]
        print >>log, 'MML Detection   :',attrib_det[:,0]

        if not abstract:
            print >>log, '\nOrig. Scales        :',attrib_org[:,1]
        print >>log, 'Detected Scales*',scale_match,':',attrib_det[:,1]


        if not abstract:
            # Source wavelet causes adding -2 to alpha
            print >>log, '\nOriginal Phase  :',attrib_org[:,3]
        print >>log, 'Detected Phases :',attrib_det[:,3]


    tm_det=__tm.time()-tm_det
    print >>log, '\n==> Done in',str(float("% .3f"%tm_det)),'second(s).'

    #==================================================================
    #        ================= Estimation =================
    #==================================================================
    print >>log, '\n   ============== Estimation ==============\n'
    if user:
        ans_est=raw_input ('Do you want to do the estimation part [yes]/no? ')
    tm_est=__tm.time()
    if ans_est.lower() not in ['n','no']:

        attrib_est=__np.zeros( (Kdet,4) )
        err=range( Kdet )
        events_est=list()
        attrib_amp=list()

        if estimate_type.upper() not in ['LSQR','LS','LM','BFGS']:
            estimate_type='BFGS'
        estimate=eval('__core.Steps.'+estimate_type.upper()+'estimate')    

        for ind in range(Kdet):
            guess=attrib_det[ind,:]
            guess[2]=_df_Ai
            thr_scl=30
            thr_alp=6

            OUT= estimate(events_det,guess,ind=ind,**est1_dict)
            if (OUT[0][1::] !=  3*[_Null]).all():
                if (abs(OUT[0][1])<thr_scl and abs(OUT[0][2])<thr_alp):
                    guess=OUT[0]
                OUT= estimate(events_det,guess,ind=ind,**est2_dict)
                if (abs(OUT[0][1])>thr_scl or abs(OUT[0][2])>2*thr_alp):
                    OUT=list(OUT) # In case of divergence reload initial guess
                    OUT[0]=__np.array([attrib_det[ind,0],_Null,_Null,_Null])

            attrib_est[ind,:],err[ind],f,ind=OUT
            err[ind]=err[ind][-1]
            # Form event amplitude ...
            if (OUT[0][1::] !=  3*[_Null]).all():
                f.update(OUT[0])
                event=f.data(tshift=0,tnorm=1)[0]
            else:
                event=__np.zeros(N)
            events_est.append(__core.Misc.norm(events_det[ind,:])*event)
            attrib_amp.append(__core.Misc.norm(events_det[ind,:]) *max(abs(event)))


        if not abstract:
            # Source wavelet causes adding -2 to alpha 
            print >>log, '\nOrignal Scales & Orders(source included) : ',
            print >>log, '\n',attrib_org[:,1],'\n',attrib_org[:,2]

        print >>log, '\nEstimated Scales & Orders  :   ',
        print >>log, '\n',attrib_est[:,1],'\n',attrib_est[:,2]
        
        
        attrib_amp=__np.array(attrib_amp)
        if Kdet > 0:
            events_est=__np.array(events_est)
            trace_est=__np.sum(events_est,axis=0)
        else:
            events_est=__np.zeros(N)
            trace_est=events_est

        print >>log, '\n\nError      : ',err
        attrib_vect_est = _Null*__np.ones((4,N))
        ind = __np.int32(attrib_est[:,0])
        attrib_vect_est[0,ind]=attrib_amp
        attrib_vect_est[1,ind]=attrib_est[:,1]
        attrib_vect_est[2,ind]=attrib_est[:,2]
        attrib_vect_est[3,ind]=attrib_est[:,3] %2

    tm_est=__tm.time()-tm_est
    print >>log, '\n==> Done in',str(float("% .3f"%tm_est)),'second(s).'


    #==================================================================
    # ====================== Save data files ========================
    #==================================================================
    if not __os.access(pydata,__os.F_OK):
        __os.makedirs(pydata)
    print >>log, '\n   =========== Save data files ============\n'
    # Prompt user for saving data to file
    if user:
        ans_sav=raw_input (
            'Do you want to save this reflectivity model to file [yes]/no? ')
    if ans_sav.lower() not in ['n','no']:
        if abstract:
            model_file='real_model.pyd'
        else:
            model_file='synt_model.pyd'
        mymodel=open(__path.join(pydata,model_file),'w+') 
        # Seismic reflectivity model
        mywav=open(__path.join(pydata,'waves.pyd'),'w+')  
        # Source wavelet and CWT 
        __cpk.dump(wave_dict,mywav)
        __cpk.dump(model_dict,mymodel)
        mymodel.close()
        mywav.close()
        events_det.dump(__path.join(pydata,'events_det.pyd'))
        attrib_det.dump(__path.join(pydata,'attrib_det.pyd'))

        if ans_est.lower() not in ['n','no']:
            attrib_est.dump(__path.join(pydata,'attrib_est.pyd'))
            attrib_vect_est.dump(__path.join(pydata,'attrib_vect_est.pyd'))
            events_est.dump(__path.join(pydata,'events_est.pyd'))
            trace_est.dump(__path.join(pydata,'trace_est.pyd'))

    print >>log, '\n==> Done.\n'

    #==================================================================
    #        =============   Graphical Results   =============
    #==================================================================
    print >>log, '\n   =========== Graphical Results  ============\n'
    if __py == None:
        print >>log, __chr.__Disp_Err

    # Prompt user for showing results as figures
    if user and (__py != None):
        ans_fig=raw_input (
            'Do you want to See the results as figures yes/[no]? ')
    if ans_fig.lower() in ['y','yes'] and (__py != None) :
        print >>log, 'Close the figure to continue (If others are not shown)'
        _tsize=__py.rcParams['axes.titlesize']
        _lsize=__py.rcParams['axes.labelsize']
        if not abstract:
            #__py.figure();__py.plot(attrib_vect[0,:],'k',lw=1);
            #__py.title('Location and amplitude of events',size=_tsize);
            #__py.xlabel(_df_xlabel,size=_lsize)
            #__py.ylabel('Amplitude',size=_lsize)


            __py.figure();__py.title('Actual events',size=_tsize)
            numevents=events.shape[0]
            ev=list(range(numevents))
            for k in __np.arange(numevents):
                x=events[k,:]
                ev[k],=__py.plot(x,label='event'+str(k),lw=1)
                __tm.sleep(0.01)
            __py.xlabel(_df_xlabel,size=_lsize)
            __py.ylabel('Amplitude',size=_lsize)


        __core.Cwt.mmlimage(max_points,mml,acw,
                           scatter=(attrib_det.shape[0] !=0))

        __py.figure();__py.title('Segmentation of detected events',size=_tsize)
        wev=list(range(Kdet))
        for k in __np.arange(Kdet):
            x=events_det[k,:]
            m=masks[k,:]*max(trace)
            if abstract or len(ev)-1<k:
                wev[k],=__py.plot(x,label='event'+str(k),lw=1)
            else:
                wev[k],=__py.plot(x,c=ev[k].get_color(),label='event'+str(k),
                                lw=1)
            __py.plot(m,'--',c=wev[k].get_color(),label='window'+str(k),lw=1)
            __tm.sleep(0.01)      
        __py.xlabel(_df_xlabel,size=_lsize)
        __py.ylabel('Amplitude',size=_lsize)



        if ans_est.lower() not in ['n','no'] and attrib_det.shape[0] !=0:
            tmp_alpha=0*_Null * __np.ones(N)
            tmp_alpha[__np.int32(attrib_est[:,0])]=attrib_est[:,2]
            tmp_alpha[__np.logical_and(tmp_alpha<-10,tmp_alpha!=_Null)]=-10
            rng=range(N)
            __py.figure()
            legloc=2
            lsize=14
            __py.plot(rng,tmp_alpha,'-or')
            if not abstract:
                tmp_alp_vect=attrib_vect[2,:]
                tmp_alp_vect[tmp_alp_vect==_Null]=0
                __py.plot(rng,tmp_alp_vect,'-sb',alpha=0.65)

            __py.legend(['Estimated alpha','Real alpha'],loc=4)
            __py.title("Error in Estimation of Alphas",size=_tsize)
            __py.ylabel('Fractional Order',size=_lsize)
            __py.xlabel(_df_xlabel,size=_lsize)


        __py.figure()
        __py.plot(trace,lw=1)
        __py.plot(trace_det,'-.g',lw=1)
        if ans_est.lower() not in ['n','no'] and attrib_det.shape[0] !=0:
            __py.plot(trace_est,'--r',lw=1)
        __py.legend(['Actual','Detected','Estimated'])
        __py.title("Actual and estimated trace",size=_tsize)
        __py.ylabel('Amplitude',size=_lsize)
        __py.xlabel(_df_xlabel,size=_lsize)
        __py.show()


    print >>log, '\n==> Done.'


    print >>log, '\n\n        Summary'
    print >>log, ' ====================\n'
    if not abstract:
        print >>log, 'Number of original events : ',attrib_org.shape[0]
    print >>log, 'Number of detected events : ',attrib_det.shape[0],'\n'
    print >>log, 'Synthesiz/Load  Time (s) : ',str(float("% .3f"%tm_syn))
    print >>log, 'Detection       Time (s) : ',str(float("% .3f"%tm_det))
    print >>log, 'Estimation      Time (s) : ',str(float("% .3f"%tm_est))
    print >>log, '\n==> Total Processing Time (s) : ',str(
        float("% .3f"%(tm_est+tm_det+tm_syn)))

    DICT=locals()
    if Kdet == 0:
        print >>__stderr," ===> Warning: No event was detected. The trace is either too smooth or has no events."
    return DICT



# ************************************************************
#                      Generation Function
# ************************************************************
def generate(type='rsf',input=_df_input,param=None,log=__stderr):
    """
    Generate / Load seismic data(trace) for analysis
    
    Input:
    
        type  : Type of input data ['rsf']/'trace'/'new'/'loadreal'/'loadsynt'
                rsf      : get data from rsf binary files
                trace    : use givenarray/matrix by input argument as trace
                new      : Generate a new synthetic trace
                loadreal : Reload stored real model from last run
                loadsynt : Reload stored real model from last run 
                
        input : RSF Data file or Seismic trace (according to type argument)

        This may be required based on type argument :
        param : Trace number for rsf files and dictionary of synthesize 
                parameters for type='new' as below)

                N    : # Length of trace                
                K    : # Number of events
                Delta: # Minimum distance between events
                amp  : # Amp Ranges 
                rnga : # Alpha Ranges
                rngp : # Phase Range
                rngs : # Scale Range for gaussian trace generator
                scl  : # Scale for synthesis with splines and source wavelet
                regul: # Regularized spacing between events
                uniform: # Uniform amplitude for all events

    Output:
        Seismic Trace and Attributes if applicable (None otherwise)
        model_dict  : Dictionary of trace and its model parameters
        wave_dict   : Dictionary of wavelets and their scales
        attrib_org  : Attributes of original events (Synthetic data)
        attrib_vect : Vector of Attribute (Synthetic data)
        tm_syn      : Elapsed time (s) for data loading/generation
    """

    #---------$ Synthetic Data Defaul Parameters $---------
    synt_dict={'scl':8,'rngs':[8,12],'rnga':__np.array([-1,1]),'rngp':[0,2],
               'amp':[4,10],'regul':0,'uniform':0}

    engine='gauss'     # Generate synth data from 'spline' or 'gauss'
    manaul=0           # Synthesize with pre-defined values


    N=1000;K=11;Delta=65
    N+=(N%2)
    #---------$  MANAUL  CONTROL  $---------
    "This is used to generate event(s) with certain attributes"
    if manaul:
        N=512
        loc=44 + N/2
        amp=4.25
        sigma=19
        alpha=0.78
        phi=0.23

   
    #---------$  Create Seismic trace   $---------
    print >>log, '\n\n   ====== Generate/Load Seismic data ======\n'
    tm_syn=__tm.time()
    wave_dict={}
    attrib_org=None
    attrib_vect=None

    # === Get real seismic trace from RSF file ===
    if type.lower()=='rsf':

        if not isinstance(input,str):
            input=_df_input
        if not isinstance(param,int):
            param=_df_rx


        inputf=__core.Misc.abs_file(input, data_path)
        if inputf==None:
            raise IOError,'RSF file not found in data path: '+input
        elif not __os.access(inputf,__os.F_OK):
            raise IOError,'RSF file not found in data path: '+input

        try:
            print >>log,'Reading trace #'+str(param)+' of data from rsf to python...'
            trace,hdr=__core.API.sftrace2d(inputf,tracenum=param)	
            if len(trace)%2==1:
                trace=__np.concatenate([trace,[0]])
            N=len(trace)
            abstract=1
        except IOError or EOFError:
            raise IOError,'The File is either empty or not existing !'

        model_dict={'trace':trace,'N':N,'abstract':abstract}


    # ===  Set variables for input trace for characterization ===
    elif type.lower()=='trace':
        if  not issubclass(input.__class__,(__np.ndarray,__np.matrix)):
            raise ValueError, "Trace Mode : Input argument should be array or matrix of seismic trace !"
        trace=__np.array(input).flatten()
        if len(trace)%2==1:
            trace=__np.concatenate([trace,[0]])
        N=len(trace)
        abstract=1
        model_dict={'trace':trace,'N':N,'abstract':abstract}


    # === Generating synthetic data by splines or Gaussian waves ===
    elif type.lower()=='new':
        abstract=0
        src=__core.Cwt.wavelet('gaus',wtype='real',param={'p':2})
        wave_dict={'srcscl':synt_dict['scl'],'srcwavelet':src}
        if isinstance(param,dict):
            synt_dict.update(param)
    

        if engine.lower()=='gauss':
            # Synthetic trace with Gaussian Manifold
            print >>log, "\n Synthetic data by using Gaussian manifold(alpha-2) \n"
            synt_dict['rnga']-=2
            synt_dict.pop('scl')
            if manaul:
               OUT=__core.Synthesize.det_gausstrace(N,loc,amp,sigma,alpha-2,phi)
            else:
               OUT=__core.Synthesize.rand_gausstrace(N,K,Delta,**synt_dict)
            [trace,events,attrib_org,attrib_vect,
             k,delta,regul,model_dict]=OUT

        else: # 'splines'
            # Synthetic trace with Splines
            print >>log, "\n Synthetic data by using splines and source convolution \n"
            synt_dict.pop('rngs')
            if manaul:
                OUT=__core.Synthesize.det_trace(N,loc,amp,s0,alpha,phi,src)
            else:
                OUT=__core.Synthesize.rand_trace(N,K,Delta,**synt_dict)
            [trace,events,attrib_org,attrib_vect,
             k,delta,regul,model_dict]=OUT

            ind=__np.int32(attrib_org[:,0])
            # Match scale & Other params
            attrib_org[:,1]*=1.4
            attrib_org[:,2]-=2
            attrib_vect[1,ind]= attrib_org[:,1]
            attrib_vect[2,ind]= attrib_org[:,2]

        model_dict.update({'N':N,'K':K,'Delta':Delta,'attrib_org':attrib_org
                           ,'attrib_vect':attrib_vect,'abstract':abstract})


    # === Reload previously used synthetic trace from data files ===
    elif type.lower()=='loadsynt':
        try: # Synthetic  model
            mymodel=open(__path.join(pydata,'synt_model.pyd'),'r+')
            model_dict=__cpk.load(mymodel)
            mymodel.close()
        except IOError or EOFError:
            raise IOError,'The File is either empty or not existing !'

        try:
            mywave=open(__path.join(pydata,'waves.pyd'),'r+') # Wavelets and ...
            wave_dict=__cpk.load(mywave)
            mywave.close()
        except EOFError or IOError:
            wave_dict={}

        for (name, value) in model_dict.items():
            exec('%s = value' %name)

    # === Reload previously used real trace from data files ===
    elif type.lower()=='loadreal':
        try:
            mymodel=open(__path.join(pydata,'real_model.pyd'),'r+') # Real model
            model_dict=__cpk.load(mymodel)
            mymodel.close()
            trace=model_dict['trace']
            N=model_dict['N']
            abstract=1
        except IOError or EOFError:
            raise IOError,'The File real_model.pyd is not existing!'

    else:
        raise ValueError," Wrong Type : The type of input argument is not any of 'rsf', 'trace', 'new', 'loadreal', or 'loadsynt' "

    tm_syn=__tm.time()-tm_syn    
    print >>log, '\n==> Done in',str(float("% .3f"%tm_syn)),'second(s).'
    return model_dict,wave_dict,attrib_org,attrib_vect,tm_syn




# ************************************************************
#                Main Body (Running Script)
# ************************************************************
if __name__ == '__main__':
    DICT = char(type='rsf',input=_df_input,param=_df_rx,user=0)
    for (name, value) in DICT.items():
        exec('%s = value' %name)

