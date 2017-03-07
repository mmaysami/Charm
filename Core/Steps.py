"""

           Characterization Steps
    ===================================
    Functions for different steps of characterization e.g. detection, 
    segmentation, and estimation

    Functions
    ----------
    detect        --  Detection of  events based on CWT and MML
    detect_alt    --  Alternative detection (Obsolete)

    segment       --  Segmentation of events

    LSQRestimate  --  Attribute estimation with LSQR(Obsolete)
    LSestimate    --  Attribute estimation with Line Search (Semi-obsolete)
    LMestimate    --  Attribute estimation with Levenberg-Marquardt
    MSLMestimate  --  Attribute estimation with Multi-scale Levenberg-Marquardt


 
AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""
# __all__=['detect','detect_alt','segment','LSQRestimate','LSestimate','LMestimate','MSLMestimate','BFGSestimate']

import numpy as __np
from   sys   import stderr as __stderr
from   scipy import optimize as __opt
from   Charm import _Null,__Disp_Err,_df_xlabel
from   Charm.Core import Misc as __misc, Window as __wnd, Manifold as __manif
from   Charm.Core.Cwt import wavelet as __wavelet

try:
    import pylab as __py
except ImportError:
    from Charm import __Disp_Err
    __py=None

_df_wr=__wavelet('gaus',wtype='complex',param={'p':2})
_df_scl=2**__np.linspace(1,4,60)
_df_width = 50
_df_ord = 4
_df_ind = [1,2]
_df_TRiter=20
_df_log=__stderr
_df_dphi=-1j
#=====================================
#              Detection
#=====================================
def detect(trace,wr=_df_wr ,scl=_df_scl,major=0.15,fineloc=0,*args,**kwargs):
    """
    Multi-scale event detection by CWT

    Input:
        trace   : 1-D seismic signal which is going to be characterized
        wr      : Wavelet object to be used for Continous wavelet transform 
                  [See cwt Module]
        scl     : Scale range to be used for CWT [See cwt Module]
        major   : Pick events with CWT coefficient down to this percent of 
                  global maximum of CWT coefficients 
        fineloc : Put fine-scale location as event location if 1; otherwise 
                  location of local maxima on MML will be counted as location
 
        Other parameters that should be passed as keyword arguments are:
        con     : Connectivity Factor for MMLs [See cwt.mml(e)]
        negline : Include Negative MMLs if 1
        rad     : Topological radius for forming MMLs [See cwt.mml(e)]
        acc     : Accuracy Factor for MMLs [See cwt.mml(e)]
        search  : Portion of scales to be searched for start point of MMLs 
                  [See cwt.mml(e)] 
        mindist : Min. Distance of 2 MMLs to be deteced as separate lines
                  [See cwt.mml(e)] 
        chkdist : How far (# of offsets) from MML location, 
                  look for local maxima
        log     : stderr or opened file to print the log if there is any


    Output:      
        attrib_det : 2-D array of (location,scale,0,phi%2) pairs for 
                     detected events sorted with their location
        mml        : Modulus maxima lines (MML)
        cw         : Continuous wavelet transform (CWT)
        lnpos      : Positive MML as list (To separately access to each line)
        ext        : Extrema points of CWT
        max_points : Position & scale of maximum points along MMLs in CWT plane
    """
    if not kwargs.has_key('log'):
        kwargs.update({'log':_df_log})

    # Computing Continuous Wavelet transform 
    (mml,cw,lnpos,ext)=wr.mmle(trace,s=scl,**kwargs)
    if lnpos==[]:
        print >>kwargs['log']," Trace is either too smooth or has no events !"
        return __np.array([]),mml,cw,lnpos,ext,__np.array([])


    acw=abs(cw)
    MAXacw,(MAXtau,MAXscale),flag=__misc.max2d(acw)


    # Get location and scales from MMLs
    attrib_det=list()
    nbposl=len(lnpos)
    mindist=15
    chkdist=5
    
    if kwargs.has_key('chkdist'):
        chkdist=kwargs['chkdist']
    if kwargs.has_key('mindist'):
            mindist=kwargs['mindist']

    for j in range(nbposl):

        pl0=lnpos[j][:,0]
        pl1=lnpos[j][:,1]
        tau=pl0[0]
        N=acw.shape[0]
        acwline=__np.zeros(acw.shape)
        # Check neighborhood of MML for any local maxima
        for i in range(max(-chkdist,-min(pl0)),min(chkdist+1,N-max(pl0)-1)):
            acwline[pl0+i,pl1]=acw[pl0+i,pl1]
        mx,(tau_mx,scale_mx),flag=__misc.max2d(acwline)
        scale=scl[scale_mx]
        phi=(__np.angle( cw[tau_mx,scale_mx])/__np.pi)%2
        "#Detected phase includes the 'pi':e.g. phi*pi"
        if fineloc==0:
            tau=tau_mx

        # Adding event's attribute to the list
        if mx >major*MAXacw:
            if len(attrib_det)==0:
                '# Add the first event'
                attrib_det.append((tau,scale,__np.NaN,phi,scale_mx))
                prv_mx=acw[tau_mx,scale_mx]
            elif  (mindist <  abs(tau-attrib_det[-1][0])):
                # min( abs(pl0[0]-lnpos[j-1][0,0]),)
                '# Add a new event'
                attrib_det.append((tau,scale,__np.NaN,phi,scale_mx))
                prv_mx=acw[tau_mx,scale_mx]
            elif mx > prv_mx:
                '# Fine tune the event location' 
                attrib_det[-1]=(tau,scale,__np.NaN,phi,scale_mx)
                prv_mx=acw[tau_mx,scale_mx]


    # Changing the lists to 2-D arrays
    attrib_det=__np.array(attrib_det)
    max_points=__np.array([]) 
    if attrib_det.shape[0] != 0:
        max_points=__np.concatenate([[attrib_det[:,0]],
                                     [attrib_det[:,4]]]).transpose()
        attrib_det=attrib_det[:,0:4]    
    if attrib_det.shape[0] > 1:
        perm=__np.argsort(attrib_det[:,0])
        attrib_det=attrib_det[perm,:]

    return attrib_det,mml,cw,lnpos,ext,max_points



def detect_alt(trace,maxratio = 0.1,sclratio = 0.4,wr=_df_wr ,scl=_df_scl,
               *args,**kwargs):
    """
    Alternative method for  multi-scale event detection by CWT

    Input:
        trace    : 1-D seismic signal which is going to be characterized
        maxratio : Percentage of max. value neighborhood to form active sets
        sclratio : Percentage of fine scales to look for local maxima
	wr       : Wavelet object to be used for Continous wavelet transform 
                   [See cwt Module]
	scl      : Scale range to be used for CWT [See cwt Module]

    Output:      
        attrib_det : 2-D array of (location,scale,0,phi%2) pairs for 
                     detected events sorted with their location
        max_points : Position & scale of maximum points along MMLs in CWT plane
    """
    N=len(trace)
    # Computing Continuous Wavelet transform 
    cw=wr.cwt(trace,s=scl,out='complex')
    acw=abs(cw)

    # Computing Abs. values and 1st & 2nd order gradients
    g=__np.gradient(acw)
    h0=__np.gradient(g[0])
    h1=__np.gradient(g[1])
    d1=__np.diff(acw,n=1)
    d2=__np.diff(acw,n=2)

    eps=1e-12
    max_acw=__misc.max2d(acw)[0]
    min_g=__misc.min2d(abs(g[0]))[0]
    max_g=__misc.max2d(abs(g[0]))[0]
    thr=0.5*(min_g+max_g)

    # Forming active sets where local maxima could be found
    # maxratio = 1/10.      # ==>  ===<IMPORTANT>===
    actset0=(acw>maxratio*max_acw)
    actset1=__np.logical_and(h0[0]<eps,actset0)
    actset2=__np.logical_and( abs(g[0])<thr,actset1)


    # Picking events of trace based on local maxima
    # sclratio = 4/10.      # ==>  ===<IMPORTANT>===
    chunk= actset2[:,0:actset2.shape[1]*sclratio].any(axis=1)
    locs = __np.nonzero(chunk)[0]
    edge = __np.concatenate(( (__np.diff(locs)!=1),[0] ))

    # Determining location and scale of all singularities
    st=0
    attrib_altdet=list()
    for i in edge.nonzero()[0]:
        if st != i:
            mx,(ind_mx,scl_mx),flag = __misc.max2d( acw[locs[st]:locs[i],:] )
        else:
            scl_mx=__np.argmax(acw[locs[i],:])
            ind_mx=0
        tau=locs[st]+ind_mx
        scale=scl[scl_mx]
        # Detected phase include 'pi': e.g. (phi*pi)
        phi=( __np.angle(cw[tau,scl_mx])/__np.pi)%2
        acw_max=__misc.max2d(acw)[0]
        # Limit on scale as well max value is possible
        if not(tau==0):
            attrib_altdet.append((tau,scale,__np.NaN,phi,scl_mx))
        st=i+1

    # Picking last event if applicable        
    mx,(ind_mx,scl_mx),flag = __misc.max2d( acw[locs[st]:,:] )
    tau=locs[st]+ind_mx
    scale=scl[scl_mx]
    phi=(__np.angle(cw[tau,scl_mx])/__np.pi)%2
    if not(tau==N-1):
        attrib_altdet.append((tau,scale,__np.NaN,phi,scl_mx))

    attrib_altdet=__np.array(attrib_altdet)    
    max_points=__np.array([])
    if attrib_det.shape[0] != 0:
        max_points=__np.concatenate([[attrib_altdet[:,0]],
                                     [attrib_altdet[:,4]]]).transpose()
        attrib_altdet=attrib_altdet[:,0:4]
    return attrib_altdet,max_points

#=====================================
#            Segmentation
#=====================================
def segment(trace,attrib_det,*args,**kwargs):
    """
    Segmenting seismic trace to detected events (atoms)

    Input:
        trace      : 1-D seismic signal which is going to be characterized
        attrib_det : 2-D array of (location,scale) pairs for detected events 
                     [See detect function]

    Output:      
        events_det : 2-D array containing segmented events in rows
        masks      : 2-D array of window function (in rows) for each event
    """

    N=len(trace)
    # Number of detected events
    numevents=attrib_det.shape[0]
    # if no or one detected event, pass trace
    if numevents < 2:
        return __misc.vector(trace,'row'),__np.ones((1,N))

    # Average Sigma(Scale) of events (2 or more events)
    scl_avg=sum(attrib_det[:,1])/numevents
    rng=__np.arange(0,N)
    events_det=list()
    masks=list()

    # Loop over events for windowing each of them
    for k in __np.arange(numevents):
        loc=attrib_det[k,0]
        scl=(4 *attrib_det[k,1] + 1*scl_avg)/5.
        
        if k==0:
            " The most left event"
            #nextloc=attrib_det[k+1,0]
            #midr=0.5*(nextloc+loc)
            window=__wnd.wbutterworth(trace,attrib_det[k,0],scl,
                                     width=_df_width,ordr=4)

        elif k!=numevents-1:
            " The intermediate events "
            if numevents > 2:
                #nextloc=attrib_det[k+1,0]
                #midr=0.5*(nextloc+loc)
                #prevloc=attrib_det[k-1,0]
                #midl=0.5*(prevloc+loc)
                window=__wnd.wbutterworth(trace,attrib_det[k,0],scl,
                                         width=_df_width,ordl=4,ordr=4)

        else:
            # The most right event         #attrib_det[k,1]
            #prevloc=attrib_det[k-1,0]
            #midl=0.5*(prevloc+loc)
            window=__wnd.wbutterworth(trace,attrib_det[k,0],scl,
                                     width=_df_width,ordl=4)

        events_det.append(window[0])
        masks.append(window[1])   


        # Alternative Windowing Method (Not Advised)
        """
        temp=trace.copy()
        width=2.5*s0*attrib_det[k,1]
        cond=(abs(rng-attrib_det[k,0])> width)
        if k!=numevents-1:
            cond1=rng > ( 3*attrib_det[k+1,0]+2*attrib_det[k,0] )/5
            cond= __np.logical_or(cond,cond1)
        if k!=0:
            cond2=rng < ( 3*attrib_det[k-1,0]+2*attrib_det[k,0] )/5
            cond= __np.logical_or(cond,cond2)
        __np.putmask(temp , cond , 0)

        masks.append(~cond)
        events_det.append(temp)
        """
    events_det=__np.array(events_det)
    masks=__np.array(masks)
    return events_det,masks




#=====================================
#            Estimation
#=====================================
_df_solve=' *'

def LSQRestimate(events_det,guess,ind=0,prm_ind=_df_ind,dphi=_df_dphi,flatphi=1
                 ,niter=50,tol=3e-6,verb=1,show=0,log=_df_log,**kwargs):
    """
    (OBSOLETE, Use LM or BFGS esstimation for best results)
    LSQR Estimating attributes of single events with a Gaussian waveform
    
    Input:
        events_det : 2-D array containing segmented events in rows 
                     (from segmentation)
        guess      : initial guess for all parameters
        ind        : Index of the event to be characterized
        prm_ind    : Index of parameter to be inverted for (tau,sigma,alpha,phi)
        dphi       : Constant for d3 (w.r.t phi) [+-1j + 0 or 0.1]
        flatphi    : Makes phi in the range of [0,2)
        niter      : Number of iteration for solver
        tol        : Error tolerance to stop 
        verb       : Show the numerical results if is 1
        show       : Show the graphical results if is 1

    Output:      
        attrib_est : Estimated values for corresponding attributes of the event

        err        : Initial & Final Error value
        f          : Gaussian manifold which estimates windowed event
        ind        : Index of the event (Useful when input ind=None)
    """
    # Initializing
    print >>log, _df_solve,    
    if __py==None and show:
        print >>log, __Disp_Err
        show= show and 0
    NRM=1     # Use normalized Gaussian waveform 
    dNRM=1    # Use Normalized derivatives of Gaussian

    # Check if there is only one event
    if len(events_det.shape)==1:
        events_det=__misc.vector(events_det,type='row')
    N=events_det.shape[1]
    attrib_est=guess.copy()
    err=list()
    iter=0

    # Form Gaussian to estimate event (sigma ~ 1.4 * wavelet scale)
    f = __manif.gmanf(n=N,s=0,fi=0); f.update(guess)
    # Pick segmented event to be characterized
    In= events_det[ind,:].copy()     
    if NRM:
        In/=__misc.norm(In)
    res = In - f.data(tshift=0,tnorm=NRM)[0]
    err.append( __misc.mse(res) )        

    G=f.grad(tshift=0,tnorm=dNRM,dphi=dphi)[0]
    G=__np.matrix( G.take(prm_ind,axis=1) )
    A=(G.T*G)
    b=(G.T*__np.matrix(res).T)

    dTheta,res,rank,s=__np.linalg.lstsq(A,b)
    # dTheta,res=ss.solvers.solvelsqr(A, b, x0=None, maxiter=niter, tol=tol)

    attrib_est.put(prm_ind,guess.take(prm_ind)+dTheta)
    f.update(attrib_est)
    err.append( __misc.mse(f.data(tshift=0,tnorm=NRM)[0] , In) )
    __est_show("LSQR",ind,prm_ind,"?",guess,attrib_est,
               In,f,err,NRM,verb=verb,show=show,log=log)
    print >>log, '\n'
    return attrib_est,err,f,ind



def LSestimate(events_det,guess,ind=0,prm_ind=_df_ind,smth=0,intg=0,
               dphi=_df_dphi,flatphi=1,niter=50,tol=1e-5,verb=0,show=0,
               log=_df_log,**kwargs):
    """
    (OBSOLETE, Use LM or BFGS esstimation for best results)
    Line-Search Estimating of one/some attribute(s) of single events with a 
    Gaussian waveform keeping the rest of parameters fixed

    Input:
        events_det : 2-D array containing segmented events in rows 
                     (from segmentation)
        guess      : Initial guess/values for all parameters as a list/tuple
        ind        : Index of the event to be characterized
        prm_ind    : Index of parameter to be inverted for(tau,sigma,alpha,phi)
        smth       : Scale of Gaussian smoother (for both event and manifold)
        intg       : Fractional integration order (for both event and manifold)
        dphi       : Constant for d3 (w.r.t phi) [+-1j + 0 or 0.1]
        flatphi    : Makes phi in the range of [0,2)
        niter      : Number of iteration for solver
        tol        : Error tolerance to stop 
        verb       : Show the numerical results if is 1
        show       : Show the graphical results if is 1

    Output:      
        attrib_est : Estimated values for corresponding attributes of the event
        err        : Initial & Final Error value
        f          : Gaussian manifold which estimates windowed event
        ind        : Index of the event (Useful when input ind=None)
    """

    # Initializing
    print >>log, _df_solve,
    if __py==None and show:
        print >>log, __Disp_Err
        show=show and 0
    NRM=1     # Use normalized waveform for comparing Gaussian & seismic event
    dNRM=1    # Use Normalized derivatives of Gaussian for parametric inversion
    # Check if there is only one event
    if len(events_det.shape)==1:
        events_det=__misc.vector(events_det,type='row')
    N=events_det.shape[1]
    err   = list()
    Theta = list()
    iter=0
    s0=smth
    fi0=intg

    # Form Gaussian to estimate event (sigma ~ 1.4 * wavelet scale)
    param=__np.array(guess)
    param[3]%=2
    [ti,si,ai,pi]=param.copy()
    Theta.append(param.take(prm_ind))
    f = __manif.gmanf(n=N,tau=ti,sigma=si,alpha=ai,phi=pi,s=s0,fi=fi0)
    # Pick segmented event to be characterized
    In= events_det[ind,:].copy()     
    if NRM:
        In/= __misc.norm(In)
    In=__manif.smooth(In,s=s0,fi=fi0,tnorm=NRM)[0]
    err.append(__misc.mse(f.data(tshift=0,tnorm=NRM)[0] , In))
   
    # Loop over iterations for updating estimation values
    indmin= -1
    while iter < niter and err[iter] > tol:
        iter+=1
        # Computing Differences and Error
        E= __np.matrix(f.data(tshift=0,tnorm=NRM)[0]-In).T
        G=f.grad(tshift=0,tnorm=dNRM,dphi=dphi)[0]
        G=__np.matrix( G.take(prm_ind,axis=1) )
        # Computing Jacobian & Hessian
        J=(G.T * E)
        H=__np.matrix (__np.diag(__np.diag( G.T * G  )))
        #         H=__np.matrix (( G.T * G  ))
        # Determining new parameter values (Direction:d & Step:a0)
        try:
            I=H.I
        except:
            attrib_est=__np.array([guess[0],_Null,_Null,_Null])
            f.update(attrib_est)
            print >>log,"Warning :Singular Hessian matrix in estimation"
            return attrib_est,__np.inf,f,ind
        d=__np.matrix(I * J)
        mu=0.1
        a0=1
        param.put(prm_ind,Theta[iter-1]-a0*d)
        if flatphi:
            param[3]%=2
        f.update(param)
        E_new=f.data(tshift=0,tnorm=NRM)[0]-In
        # Find step size along direction d
        while __misc.norm(E_new) > __misc.norm(E+mu*a0*E.T*G*d) and a0>1e-3:
            a0 *= 0.50
            param.put(prm_ind,Theta[iter-1]-a0*d)
            if flatphi:
                param[3]%=2
            # Print & Show results       
            f.update(param)
            E_new= f.data(tshift=0,tnorm=NRM)[0]-In
            
        Theta.append( Theta[iter-1] - a0*d.ravel() )
        err.append(__misc.mse(E_new))

 
    param.put(prm_ind,Theta[-1])
    param[3]%=2
    f.update(param)
    attrib_est=__np.array(param)
    # Show results
    #print >>log, 'Stopped after ',str(iter),' iterations.'
    __est_show("LINE SEARCH",ind,prm_ind,iter,guess,attrib_est,
               In,f,err,NRM,verb=verb,show=show,log=log) 
    print >>log, '\n'
    return attrib_est,__np.array([err[0],err[-1]]),f,ind


def LMestimate(events_det,guess,ind=0,prm_ind=_df_ind,smth=0,intg=0,
               dphi=_df_dphi,flatphi=1,niter=100,tol=1e-5,verb=1,show=0,
               log=_df_log,**kwargs):
    """
    Levenberg-Marquardt Estimating one/some attribute(s) of single events with 
    a Gaussian waveform keeping the rest of parameters fixed (TOLERANCE is set 
    to MSE of difference)
    NOTE: Holds the best results

    Input:
        events_det : 2-D array containing segmented events in rows 
                     (from segmentation, It will normalize events)
        guess      : Initial guess/values for all parameters as a list/tuple
        ind        : Index of the event to be characterized
        prm_ind    : Index of parameter to be inverted for (tau,sigma,alpha,phi
        smth       : Scale of Gaussian smoother (for both event and manifold)
        intg       : Fractional integration order (for both event and manifold)
        dphi       : Constant for d3 (w.r.t phi) [+-1j + 0 or 0.1]
        flatphi    : Makes phi in the range of [0,2)
        niter      : Number of iteration for solver
        tol        : Error tolerance to stop 
        verb       : Show the numerical results if is 1
        show       : Show the graphical results if is 1

    Output:      
        attrib_est : Estimated values for corresponding attributes of the event
        err        : Initial & Final Error value
        f          : Gaussian manifold which estimates windowed event
        ind        : Index of the event (Useful when input ind=None)
    """
    # Initializing
    if __py==None and show:
        print >>log, __Disp_Err
        show=show and 0
    print >>log, _df_solve,
    NRM,dNRM=1,1     # Use normalized gaussian waveforms & derivative
    s0,fi0=smth,intg
    err=list()
    # Check if there is only one event
    if len(events_det.shape)==1:
        events_det=__misc.vector(events_det,type='row')
    N=events_det.shape[1]  # Length of signal
    M=len(prm_ind)         # Number of attributes for inversion

    # Form segmented event and Gaussian estimator
    param=__np.array(guess).ravel(); param[3]%=2
    f = __manif.gmanf(n=N,s=s0,fi=fi0);
    In= events_det[ind,:].copy()     
    if NRM:
        In/= __misc.norm(In)
    if __np.logical_or(s0 !=0 , fi0!=0):
        In=__manif.smooth(In,s=s0,fi=fi0,tnorm=NRM)[0]

    # Levenberg-Marquardt Algorithm
    iterc=0; xc=param.take(prm_ind); nu0=.0010
    f.update(param)
    fc=f.data(tnorm=NRM)[0]-In
    err.append( __misc.mse(fc) )

    jac=f.grad(tnorm=dNRM, tshift=0, dphi=dphi)[0]
    jac=jac.take(prm_ind,axis=1)
    gc=__np.dot(jac.transpose(),fc)
    fc=0.5*__np.dot(fc,fc)


    numf,numg,numh=1,1,0
    nu=__misc.norm(gc)
    ithist=__np.zeros((niter+1,4))
    ithist[0,:]=__np.array([err[-1],fc,0,iterc])

    while __np.logical_and(err[-1] >tol , iterc <niter):
        iterc += 1
        # Compute Trial point xt
        hc=__np.dot(jac.transpose(),jac)+ nu*__np.eye(M);
        xt=xc-__np.dot( __np.linalg.inv(hc),gc)
        # Trust Region Test
        [xp,nup,TRiter]=__trtestlm(
            f,xc,xt,fc,gc,jac,nu,nu0,In,param,prm_ind,NRM)
        if TRiter > _df_TRiter:
            #print >>log,"  WARNING: Too many iterations in TR solve."
            param.put(prm_ind,xc)
            f.update(param)
            __est_show("Levenberg-Marquardt",ind,prm_ind, iterc,guess,param, 
                       In,f,err,NRM,verb=verb,show=show,log=log)
            return param,__np.array([err[0],err[-1]]),f,ind
        # Update Current Point (xc), nu
        xc,nu=xp,nup; numf += TRiter

        # Check TR results
        if TRiter > 1:
            param.put(prm_ind,xc)
            f.update(param)
            fc=f.data(tnorm=NRM)[0]-In
            err.append( __misc.mse(fc) )
            jac=f.grad(tnorm=dNRM, tshift=0, dphi=dphi)[0]
            jac=jac.take(prm_ind,axis=1)
            gc=__np.dot(jac.transpose(),fc)
            fc=0.5*__np.dot(fc,fc)
            numf += 1; numg += 1
        ithist[iterc,:]=__np.array([err[-1],fc,TRiter,iterc])

    costdata=[numf, numg, numh]
    param.put(prm_ind,xc)
    f.update(param)
    attrib_est=param
    __est_show("Levenberg-Marquardt",ind,prm_ind,iterc,guess,attrib_est,In,f,
               err,NRM,verb=verb,show=show,log=log)
    print >>log, '\n'
    return attrib_est,__np.array([err[0],err[-1]]),f,ind



def MSLMestimate(events_det,guess,ind=0,smth=[12,9,6,3,0],tol=[1e-5,5e-6,1e-6,
                                                               5e-7,1e-7],
                 verb=1,show=0,log=_df_log,LM_dict={},**kwargs):
    """
    Multi-Scale Levenberg-Marquardt Estimating one/some attribute(s) of single 
    events with a Gaussian waveform keeping the rest of parameters fixed 
    (TOLERANCE is set to MSE of difference)

    Input:
        events_det : 2-D array containing segmented events in rows 
                     (from segmentation)
        guess      : Initial guess/values for all parameters as a list/tuple
        ind        : Index of the event to be characterized
        smth       : Scale of Gaussian smoother (for both event and manifold)
        tol        : Error tolerance to stop 
        verb       : Show the numerical results if is 1
        show       : Show the graphical results if is 1
        LM_dict    : Dictionary of arguments for LMestimate function 
                     (See it's help for details)

    Output:      
        attrib_est : Estimated values for corresponding attributes of the event
        err        : Initial & Final Error value
        f          : Gaussian manifold which estimates windowed event
        ind        : Index of the event (Useful when input ind=None)
    """
    print >>log, _df_solve,
    NRM=1
    guess0=guess.copy()
    err=range(len(smth))
    # Check if there is only one event
    if len(events_det.shape)==1:
        events_det=__misc.vector(events_det,type='row')
    LM_dict.update({'ind':ind,'verb':0,'show':0,'log':log})

    for i in range(len(smth)):
        LM_dict.update({'smth':smth[i],'tol':tol[i]})
        attrib_est,err[i],f,ind=LMestimate(events_det,guess,**LM_dict)
        guess=attrib_est.copy()


    In=events_det[ind,:]/__misc.norm(events_det[ind,:])
    In=__manif.smooth(In,s=smth[i],fi=0,tnorm=NRM)[0]

    err=__np.array([err[0][0],err[-1][-1]])
    __est_show("MS Levenberg-Marquardt",ind,LM_dict['prm_ind'],"?",guess0,
               attrib_est,In,f,err,NRM=1,verb=verb,show=show,log=log)
    print >>log, '\n'
    return attrib_est,err,f,ind





def BFGSestimate(events_det,guess,ind=0,prm_ind=_df_ind,smth=0,intg=0,
               dphi=_df_dphi,flatphi=1,niter=100,tol=1e-5,verb=1,show=0,
               log=_df_log,**kwargs):
    """
    BFGS Estimating one/some attribute(s) of single events with 
    a Gaussian waveform keeping the rest of parameters fixed (TOLERANCE is set 
    to MSE of difference)
    NOTE: Holds the best results with fastest computation

    Input:
        events_det : 2-D array containing segmented events in rows 
                     (from segmentation, It will normalize events)
        guess      : Initial guess/values for all parameters as a list/tuple
        ind        : Index of the event to be characterized
        prm_ind    : Index of parameter to be inverted for (tau,sigma,alpha,phi
        smth       : Scale of Gaussian smoother (for both event and manifold)
        intg       : Fractional integration order (for both event and manifold)
        niter      : Number of iteration for solver
        tol        : Error tolerance to stop 
        verb       : Show the numerical results if is 1
        show       : Show the graphical results if is 1

    Output:      
        attrib_est : Estimated values for corresponding attributes of the event
        err        : Initial & Final Error value
        f          : Gaussian manifold which estimates windowed event
        ind        : Index of the event (Useful when input ind=None)
    """
    # Initializing
    if __py==None and show:
        print >>log, __Disp_Err
        show=show and 0
    print >>log, _df_solve,
    
    NRM=1     # Use normalized gaussian waveforms
    s0,fi0=smth,intg
    err=list()
    # Check if there is only one event
    if len(events_det.shape)==1:
        events_det=__misc.vector(events_det,type='row')
    N=events_det.shape[1]  # Length of signal

    # Form segmented event and Gaussian estimator
    param=__np.array(guess,'float64').ravel(); param[3]%=2
    f = __manif.gmanf(n=N,s=s0,fi=fi0)
    f.update(param)
    In= events_det[ind,:].copy()     
    if NRM:
        In/= __misc.norm(In)
    if __np.logical_or(s0 !=0 , fi0!=0):
        In=__manif.smooth(In,s=s0,fi=fi0,tnorm=NRM)[0]

    err.append(__residue(param.take(prm_ind),param,prm_ind,f,In,NRM))
    xc,res,dict=__opt.fmin_l_bfgs_b(__residue,param.take(prm_ind),fprime=None
                                    ,args=(param,prm_ind,f,In,NRM),approx_grad
                                    =1,factr=tol/1e-16,maxfun=3*niter)


    # Alternatively BFGS can be used, holds quite similar results
    #     (xc,fc,gc,Hc,funcalls,gradcalls,warnflag)=__opt.fmin_bfgs(__residue,param.take(prm_ind),fprime=None,args=(param,prm_ind,f,In,NRM),gtol=1e-5,norm=__np.inf,epsilon=1.5e-8,maxiter=3*niter,full_output=1,disp=0,retall=0,callback=None)


    err.append(res)
    param.put(prm_ind,xc)
    f.update(param)
    attrib_est=param
    __est_show("BFGS",ind,prm_ind,dict['funcalls']/3,guess,attrib_est,In,f,
               err,NRM,verb=verb,show=show,log=log)
    print >>log, '\n'
    return attrib_est,__np.array(err),f,ind



#============================================
#   Internal Functions for Estimation Part
#============================================
def __residue(x,param,prm_ind,f,In,NRM=1):
    param.put(prm_ind,x)
    f.update(param)
    fc=f.data(tnorm=NRM)[0]-In
    return __misc.mse(fc)

def __trtestlm(f,xc,xt,fc,gc,jac,nu,nu0,In,param,prm_ind,NRM=1):
    """
    Trust Region test for Levenberg-Marquardt Method
    
    Input:
        f      : Function to be minimized
        xc     : Curent value of function parameters
        xt     : Trial point value of function parameters
        fc     : Funcation value at current point
        gc     : Gradient of function at current point
        jac    : Jacobian of function at current point
        nu     : Value of nu parameter
        nu0    : Value of nu0 parameter
        In     : Segmented event (waveform)
        NRM    : Use normalized Gaussian waveform 

    Output:      
        xp     : Predicted parameteres (x)
        nup    : Predicted nu value
        TRiter : Number of iteration for TR

    """

    mu0=.1; mulow=.25; muhigh=.75; wup=2.0; wdn=.50
    z=xc; N=len(z); iter=0

    while __np.logical_and( (z==xc).all() , iter <= _df_TRiter):
        param.put(prm_ind,xt)
        f.update(param)
        ft=f.data(tnorm=NRM)[0]-In
        ft=0.5*__np.dot(ft,ft)

        iter += 1
        s=xt-xc; 
        ared = fc-ft; 
        pred = -0.5*__np.dot(gc.transpose(),s)
        ratio  = ared/pred
            
        if ratio < mu0:
            nu=max(nu*wup,nu0)
            hc=__np.dot(jac.transpose(),jac)+ nu*__np.eye(N)
            xt=xc-__np.dot(__np.linalg.inv(hc),gc)
        elif ratio < mulow:
            nu=max(nu*wup,nu0)
            z=xt
        else:
            z=xt
            if ratio > muhigh:
                nu=wdn*nu
                if nu < nu0:
                    nu=0.0

    nup=nu; xp=z; TRiter=iter;
    return (xp,nup,TRiter)





def __est_show(method,ind,prm_ind,iter,guess,param,In,f,err,NRM=1,
               verb=1,show=0,sepfig=0,log=_df_log):
    """
    INTERNAL USE Only for printing and showing estimation results
    Used in all **estimate functions

    """
    _tsize=__py.rcParams['axes.titlesize']
    _lsize=__py.rcParams['axes.labelsize']
    __py.rc('legend',fontsize=12)
    # Print  initial numerical results
    if verb:
        print >>log, "\n___________________________________"
        print >>log, "         "+str(method.upper())+"           \n"
        print >>log, '    < INITIAL Guess >'
        print >>log, '  Event # : ',ind
        print >>log, '  Parameter index = ',prm_ind
        print >>log, '  [ti,si,ai,pi]   = ',__np.array(guess)
        print >>log, '  Error Value     = ',err[0],'\n'

    # Show initial graphical results
    if show:
        g=f.copy()
        g.update(guess)
        # Legend position & text size
        legloc=0 #2:top-left, 0:best
        nfig=__py.figure()
        nfig=nfig.number
        __py.clf()
        __py.plot(In,'--b',g.data(tshift=0,tnorm=NRM)[0],'r',lw=1)
        __py.legend(['Windowed Event','Initial guess'],loc=legloc)
        __py.title('Parameter pre_estimation for '+str(ind+1)+
                 'th event, Initial guess',size=_tsize)
        __py.ylabel('Normalized amplitude',size=_lsize)
        __py.xlabel(_df_xlabel,size=_lsize)
        #__py.show()

    # Print numeric results of final iteration
    if verb:
        print >>log, '   < FINAL Values > ' 
        print >>log, '  Iteration # ',iter
        print >>log, '  [ti,si,ai,pi]   = ',param
        print >>log, '  Error Value     = ',err[-1]
            
    # Show graphical results of final iteration
    if show:
        g.update(param)
        if sepfig:
            __py.figure(nfig+1)
        else:
            __py.figure(nfig)
        __py.clf()
        __py.plot(In,'--b',g.data(tshift=0,tnorm=NRM)[0],'r',lw=1)
        __py.legend(['Windowed event','Estimation'],loc=legloc)
        __py.title('Parameter estimation for '+str(ind+1)+
                 'th event, Final results',size=_tsize)
        __py.ylabel('Normalized amplitude',size=_lsize)
        __py.xlabel(_df_xlabel,size=_lsize)
        #__py.show()
    pass
