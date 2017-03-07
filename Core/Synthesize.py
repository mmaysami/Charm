"""

            Synthetize
    ==========================
    Generation of synthetic 1-D reflectivity profiles with fractional 
    order transitions

    Functions
    ----------
    spline           -- Fractional spline 
    fspline          -- Fractional splines in frequency domain
    gfspline         -- Generalized Fractional splines in frequency domain

    spiketrain       -- Sike train
    fsplinetrain     -- Train of fractional splines with different orders
    src_conv         -- Convolution of reflectivity with seismic source to 
                        get seismic trace
    
    det_trace        -- Deterministic Trace with fractional splines
    rand_trace       -- Random Trace with fractional splines
    det_gausstrace   -- Deterministic Trace with Gaussian manifolds
    rand_gausstrace  -- Random Trace with Gaussian manifolds




AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""
# __all__=['spline','fspline','gfspline','spiketrain','fsplinetrain','src_conv','det_trace','rand_trace','det_gausstrace','rand_gausstrace']

from Charm.Core import Cwt as __wav, Misc as __misc, Manifold as __manf
from Charm      import _eps, _Null
from sys        import stderr as __stderr
import numpy     as __np
import numpy.fft as __ft
import time      as __tm


_df_amp=(2,10)
_df_scl=10
_df_rnga=(-0.85,-0.05)
_df_rngs=(8,12)
_df_rngp=(0,2)
_df_src=__wav.wavelet('gaus',wtype='real',param={'p':2})

#=====================================
#  Define Spline Functions
#=====================================
def spline(N,alpha,loc=None,caus=1):
    """
    Generates a Single spline with certain order and location (Obsolete)
    The spline shows (alpha-1) order.

    Input:
        N     : Number of samples
	alpha : >0,Singularity order of spline (in [0 .. 1] for our purpose)
	loc   : Location of spline w.r.t. origin [0..N]
	caus  : Causal flag, returns Casual splines 

    Output:
        Spline Function (Returns Zero for Alpha <0)
    """
    if N%2 == 1:
	N+=1
    #     if (alpha <-1):
    #         return __np.zeros(N)
    if (loc > N)| (loc==None):
	loc=N/2+1

    if caus==1:
	t = __np.arange(_eps,N-loc+1,dtype='float')
	sp = t**(alpha)
        if (alpha == -1) or (alpha==- __np.inf):
            sp[0] = 0
        if loc==0:
            Z=__np.array([])
        else:
            Z=__np.zeros(loc-1)
    	sp = __np.concatenate((Z,sp))
	sp=__np.diff(sp)
	sp =__np.concatenate([[0],sp])
        sp /= max(abs(sp))+1e-16*(max(sp)==0)
    else:
	t = __np.arange(_eps,loc,dtype='float')
	sp = t**(alpha)
	sp = __np.concatenate((sp[loc-1::-1],__np.zeros(N-loc)))
	sp = __np.diff(- sp)
	sp =__np.concatenate([sp,[0]])
        sp /= max(abs(sp))+1e-16*(max(sp)==0)    
    return sp

def  fspline(N,alpha,phi=None,loc=None,tnorm=1,domain=None):
    """
    (alpha, tau) Fractional Spline with certain order, phase, and 
    location in Frequency Domain

    Note: tau may be shift(integer values) or phase rotation
    Note: Due to sampling problems, it might be better to be used 
          with some source wavelets

    Input:
        N      : Number of samples
	alpha  : Singularity order of spline (in [-1 .. _eps] for our purpose)
        phi    : Phase rotation component. {[a,a+2] | a in R} covers all shifts
	loc    : Location of spline w.r.t. origin [0..N]
        tnorm  : L2 time normalized output if True
        domain : domain of output data ['time']/'freq'/None (half freq. range)

    Output:
        Fourier transform of Spline over non negative frequency range
    """
    if N%2 == 1:
	N+=1
    if (loc > N) | (loc==None):
	loc=N/2+1

    alpha = -(alpha)
    if phi==None:
        phi=0.5*alpha

    w = 2*__np.pi/N * __np.arange(_eps,N/2+1)
    #To avoid NaN error :  (-1j*w)**(0.5*alpha-phi)  *  (1j*w)**(0.5*alpha+phi)
    spf = __np.exp(-1j*w*loc) * (w**alpha) * (-1+0j)**phi

#     if alpha <= 0:
    spf[0]=__np.exp(-1j*w[0]*loc)
    spf=__misc.outtype( spf, domain=domain, mode='complex',tshift=0, 
                        tnorm=tnorm, stype='real')[0]
    return spf

def gfspline(N,alpha,phi=None,loc=None,tnorm=1,domain=None):
    """
    Generalized (alpha,tau) Fractional Spline with certain order, phase, 
    and location in Frequency Domain
    
    Note: tau may be shift(integer values) or phase rotation
    Note: Due to sampling problems, it might be better to be used 
          with some source wavelets

    Input:
        N      : Number of samples
	alpha  : Singularity order of spline (in [-1 .. _eps] for our purpose)
        phi    : Shift/phase rotation component (shift for integers) 
	loc    : Location of spline w.r.t. origin [0..N]
        tnorm  : L2 time normalized output if True
        domain : domain of output data ['time']/'freq'/None (half freq. range)

    Output:
        Fourier transform of Spline over non negative frequency range
    """
    if N%2 == 1:
	N+=1
    if (loc > N) | (loc==None):
	loc=N/2+1

    alpha = -(alpha)
    if phi==None:
        phi=0.5*alpha


    w = 2*__np.pi/N * __np.arange(_eps,N/2+1)

    trm1=(__np.exp(1j*w)-1)/(1j*w)
    trm2=(1-__np.exp(-1j*w))/(1j*w)
    spf =trm1**(0.5*(alpha)-phi) * trm2**(0.5*(alpha)+phi)  * __np.exp(-1j*w*loc)
    # Limit of spline in w=0 is 1
    spf[0]=__np.exp(-1j*w[0]*loc)

    spf=__misc.outtype( spf, domain=domain, mode='complex',tshift=0, 
                        tnorm=tnorm, stype='real')[0]
    return spf


#=====================================
#  Define Spike/Spline Trains
#=====================================
def spiketrain(N,k=None,delta=None,amp=_df_amp,regul=0,uniform=0):
    """
    Generates spike train with certain features
    
    Input:
        N       : Number of samples
        k       : Number of spikes
        delta   : Minimum distance between spikes (delta >= 0.8*N/k leads 
                  to regular spacing)
        amp     : abs(amplitude) range for spikes, amplitude range 
                  in +/-(minamp,maxamp)
        regul   : Flag for regular spacing between spikes
        uniform : Flag for uniform amplitude of spikes

    Output:
        Tuple of output arguments
        r     : Spike train (Spiky Reflectivity)
        k     : Number of spikes (useful when k=None as input)
        delta : Actual minimum distance between spikes
        regul : Regular spacing status flag 
    """    
    if delta==None:
        delta=1
    if k==None:
        k=round(N/10 * __np.random.rand())
    elif k==1:
        delta=1
    minamp=min(abs(__np.array(amp)))
    maxamp=max(abs(__np.array(amp)))

    regflag=0
    r=__np.zeros(N)
    if (N/k < delta) & (regul ==0):
        print 'Warning: Invalid delta > N/k: Min. Distance is not feasible!'
        return (r,k,delta,regul)
    elif (delta < __np.floor(0.8 * N/k)) & (delta > 0) & (regul ==0):
        tm0=__tm.time()
        for p in __np.arange(k):
            assignflag = 0
            while assignflag == 0:
                margin=N/15.
                n = __np.floor(__np.random.rand() * (N-2*margin))+margin
                nnz=r.nonzero()[0]
                if  (sum(abs(nnz-n)>delta) == p) | (p==0):
                    r[n] = 1
                    assignflag = 1
                elif __tm.time()-tm0 >5:
                    assignflag = 1
                    k-=1
                    print >>__stderr,"1 Spike is skipped !"
        nnz=r.nonzero()[0]

    else: #[regul=1 | (N/k>=delta)] & [(N/k<delta+1)|(floor(0.8*N/k)>=delta)]
        if not regul:
            print 'Warning: Regular Spacing is chosen automatically for generated data!'
            regul=1
        if (N/k < delta):
            delta=__np.floor(N/k)
            print 'Invalid delta > N/k: Min. Distance is not feasible.'
            print 'Delta will be changed to the maximum value for regular spacing.'
        offset=__np.floor((N-k*delta)/2);
        nnz=__np.floor(__np.linspace(offset,N-1-offset,k)).astype('int')


    print len(nnz),k
    if not uniform:
        r[nnz]=__np.sign(2*__np.random.rand(k)-1)*(__np.random.rand(k) * 
                                                   (maxamp-minamp)+minamp) 
    else:
        amp=__np.random.rand() * (maxamp-minamp)+minamp
        r[nnz] = __np.sign(__np.random.rand(k)-0.5)*amp

    if k>1:
        delta=min(__np.diff(r.nonzero()[0]))
    else:
        delta=N
    r=abs(r)        # %%%%%%%%%%%%%%%%% ABS (SPIKES) %%%%%%%%%%%%%%%%%%%
    return (r,k,delta,regul)



def fsplinetrain(N,k=None,delta=None,rnga=_df_rnga,rngp=_df_rngp,amp=_df_amp
                 ,regul=0,uniform=0,caus=1):
    """
    Generates fractional spline train with certain features 
    (Use tracegen instead)
    
    Input:
        N       : Number of samples
        k       : Number of spikes
        delta   : Minimum distance between spikes (delta >= 0.8*N/k leads 
                  to regular spacing)
        rnga    : range of spline order values
        rngp    : range of instantaneous phase component values
        amp     : abs(amplitude) range for spikes, amplitude range 
                  in +/-(minamp,maxamp)
        regul   : Flag for regular spacing between spikes
        uniform : Flag for uniform amplitude of spikes
	caus    : Causal flag, Casual splines if 1 

    Output:
        Tuple of output arguments
        r     : Spline train (Reflectivity)
        r_spk : Equivalent Spike train (Spiky Reflectivity)
        alpha : Fractional orders of events (zero elsewhere)
        phi   : Phase shift of events (zero elsewhere)
        k     : Number of spikes (useful when k=None as input)
        delta : Actual minimum distance between spikes
        regul : Regular spacing status flag 
    """
    if N%2 == 1:
	N+=1
    (r_spk,k,delta,regul)=spiketrain(N,k,delta,amp,regul,uniform)
    
    alpha=None
    lb=__np.min(__np.array(rnga));ub=__np.max(__np.array(rnga))
    a=(lb-ub)*__np.random.rand(k)+ub
    phi=None
    lb=__np.min(__np.array(rngp));ub=__np.max(__np.array(rngp))
    p=((lb-ub)*__np.random.rand(k)+ub)%2 #-(a+1)/2.
    if r_spk != None :
        ind = r_spk.nonzero()[0]
        alpha= __np.zeros(N)
        phi  = __np.zeros(N)
        alpha[ind]=a
        phi[ind]=p

        r=__np.zeros(N/2+1)
        for i in range(k):
            r=r+r_spk[ind[i]]* fspline(N,a[i],p[i],ind[i],tnorm=0,domain=None)  
        r=__misc.outtype(r,domain='time',mode='complex',tshift=0,tnorm=0,
                         stype='real')[0]
        return (r,r_spk,alpha,phi,k,delta,regul)
    else:
        return (None,None,None,None,k,delta,regul)


#=====================================
#  Define Seismic Trace Generators
#=====================================
_doc_out="""

    Output:
        Tuple of output arguments
        trace : Seismic trace
        events: 2D Array having single events in columns
        attrib_org  : Attributes of single events (location,scale,order,phase) 
                      without source effect (events attribs in rows)
        attrib_vect : Matrix with vector of attributes 
                     (amplitude,scale,order,phase) with zero elsewhere 
                     (without source effect) in columns
        k     : Number of spikes (useful when k=None as input)
        delta : Actual minimum distance between spikes
        regul : Regular spacing status flag 
        model_dict : Dictionary of model parameters
"""

def rand_trace(N,k=None,delta=None,scl=_df_scl,rnga=_df_rnga,
               rngp=_df_rngp,src=_df_src,amp=_df_amp,regul=0,uniform=0):
    """
    Generates semi-random seismic trace  with certain features using fractional 
    splines.Seismic source waveform is positive i.e.           __  
                                            (@ zero has max \_/  \_/ )
    
    Input:
        N       : Number of samples
        k       : Number of spikes
        delta   : Minimum distance between spikes (delta >= 0.8*N/k leads 
                  to regular spacing)
        scl     : scale of source signature to be used
        rnga    : range of spline order values
        rngp    : range of instantaneous phase component values
        amp     : abs(amplitude) range of singularities, amplitude range 
                  in +/-(min,max) amp
        src     : a wavelet class object to define source signature 
        regul   : Flag for regular spacing between spikes
        uniform : Flag for uniform amplitude of spikes
    """    
    if N%2 == 1:
	N+=1
    (r_spk,k,delta,regul)=spiketrain(N,k,delta,amp,regul,uniform)
    model_dict={'k_act':k,'delta_act':delta,'regul':regul}

    ###############   r_spk : POSITIVE +++++

    alpha=None
    lb=__np.min(__np.array(rnga));ub=__np.max(__np.array(rnga))
    a=(lb-ub)*__np.random.rand(k)+ub
    phi=None
    lb=__np.min(__np.array(rngp));ub=__np.max(__np.array(rngp))
    p=((lb-ub)*__np.random.rand(k)+ub)%2 #-(a+1)/2.



    if r_spk != None :
        attrib_org=__np.array([r_spk.nonzero()[0],scl*__np.ones(k),a,p]
                              ).transpose()
        attrib_vect=_Null * __np.ones((4,N))
        ind = r_spk.nonzero()[0]
        attrib_vect[0,ind]=r_spk[ind]
        attrib_vect[1,ind]=scl
        attrib_vect[2,ind]=a
        attrib_vect[3,ind]=p

        trace=__np.zeros(N)
        events=__np.zeros((k,N))
        for i in range(k):
            # Number of generated events
            Fr=fspline(N,a[i],p[i],ind[i],tnorm=0,domain=None)
            wav = src.wavf(N,s=scl,tnorm=0,domain=None)[0]
            Ftrace = __misc.outtype(Fr*wav,domain='time',mode='complex',
                                    tshift=0,tnorm=1,stype='real',axis=0)[0]
            events[i,:]=r_spk[ind[i]]*Ftrace /max(abs(Ftrace)) 
            # for Amp. same as spikes
            trace = trace+events[i,:]
        print 'All',k,' events are generated'


        model_dict.update({'trace':trace,'events':events,'attrib_org':
                           attrib_org,'attrib_vect':attrib_vect})
        return [trace,events,attrib_org,attrib_vect,k,delta,regul,model_dict]
    else:
        return [None,None,None,None,k,delta,regul,model_dict]



def det_trace(N,loc,amp,scl=_df_scl,alpha=0,phi=0,src=_df_src):
    """
    Generates seismic trace  with certain features using fractional splines
    Seismic source waveform is positive i.e.                   __  
                                            (@ zero has max \_/  \_/ )
    
    Input:
        N       : Number of samples
        loc     : Location of events (assumed as only nonzero points)
        amp     : Maximum amplitude for spikes, amplitude range in [-amp,+amp]
        scl     : scale of source signature to be used
        alpha   : Singularity order coressponding to each event
        phi     : Phase rotation component: [0,2)
        src     : a wavelet class object to define source signature 
    """    
    if N%2 == 1:
	N+=1

    loc = __np.array([loc]).ravel()
    amp = __np.array([amp]).ravel()
    alpha = __np.array([alpha]).ravel()
    phi = __np.array([phi]).ravel()
    k=loc.shape[0]


    delta=N
    attrib_org=__np.array([loc,scl*__np.ones(k),alpha,phi]).transpose()

    if k>1:
        perm=__np.argsort(attrib_org[:,0],axis=0)
        attrib_org=attrib_org[perm,:]
        loc=loc[perm]
        delta=__np.diff(loc)
        if k>2:
            delta=min(delta)
        if amp.shape[0]==1:
            amp=amp *__np.random.rand(k)
        else:
            amp=amp[perm]

            
    attrib_vect=_Null*__np.ones((4,N))
    attrib_vect[0,loc]=amp
    attrib_vect[1,loc]=scl
    attrib_vect[2,loc]=attrib_org[:,2]
    attrib_vect[3,loc]=attrib_org[:,3]
    
    trace=__np.zeros(N)
    events=__np.zeros((k,N))
    for i in range(k):
        Fr=fspline(N,attrib_org[i,2],attrib_org[i,3],loc=attrib_org[i,0],
                   tnorm=0)          
        wav = src.wavf(N,s=scl,tnorm=0,domain=None)[0]
        Ftrace = __misc.outtype(Fr*wav,domain='time',mode='complex',tshift=0,
                                tnorm=1,stype='real',axis=0)[0]
        events[i,:]=amp[i]*Ftrace/max(abs(Ftrace))
        trace = trace+events[i,:]
    
    regul=None
    model_dict={'trace':trace,'events':events,'attrib_org':attrib_org,
                'attrib_vect':attrib_vect, 'k_act':k,'delta_act':delta,
                'regul':regul}
    return [trace,events,attrib_org,attrib_vect,k,delta,regul,model_dict]




def rand_gausstrace(N,k=None,delta=None,rngs=_df_rngs,rnga=_df_rnga,
                    rngp=_df_rngp,amp=_df_amp,regul=0,uniform=0):
    """
    Generates semi-random seismic trace  with certain features (by Gaussian 
    waves). Seismic source waveform is positive i.e.           __  
                                            (@ zero has max \_/  \_/ )

    Input:
        N       : Number of samples
        k       : Number of spikes
        delta   : Minimum distance between spikes (delta >= 0.8*N/k leads 
                                                   to regular spacing)
        rngs    : scale of source signature to be used
        rnga    : range of spline order values
        rngp    : range of instantaneous phase component values
        src     : a wavelet class object to define source signature 
        amp     : abs(amplitude) range for singularities, amplitude range 
                  in +/-(minamp,maxamp)
        regul   : Flag for regular spacing between spikes
        uniform : Flag for uniform amplitude of spikes
    """    
    if N%2 == 1:
	N+=1
    (r_spk,k,delta,regul)=spiketrain(N,k,delta,amp,regul,uniform)
    model_dict={'k_act':k,'delta_act':delta,'regul':regul}
    ###############   r_spk : POSITIVE +++++


    sigma=None
    lb=__np.min(__np.array(rngs));ub=__np.max(__np.array(rngs))
    s=(lb-ub)*__np.random.rand(k)+ub
    alpha=None
    lb=__np.min(__np.array(rnga));ub=__np.max(__np.array(rnga))
    a=(lb-ub)*__np.random.rand(k)+ub
    phi=None
    lb=__np.min(__np.array(rngp));ub=__np.max(__np.array(rngp))
    p=((lb-ub)*__np.random.rand(k)+ub) %2 #-(a+1)/2.


    if r_spk != None :

        attrib_org=__np.array([r_spk.nonzero()[0],s,a,p]).transpose()
        attrib_vect=_Null * __np.ones((4,N))
        ind = r_spk.nonzero()[0]
        attrib_vect[0,ind]=r_spk[ind]
        attrib_vect[1,ind]=s
        attrib_vect[2,ind]=a
        attrib_vect[3,ind]=p

        trace=__np.zeros(N)
        events=__np.zeros((k,N))
        for i in range(k):
            # Number of generated events
            f=__manf.gmanf(N,ind[i],s[i],a[i],p[i])
            events[i,:]=r_spk[ind[i]]*f.data(tnorm=1)[0] /max(abs(
                    f.data(tnorm=1)[0] )) # for Amp. same as spikes
            trace = trace+events[i,:]
        print 'All',k,' events are generated'

        model_dict.update({'trace':trace,'events':events,'attrib_org':
                           attrib_org,'attrib_vect':attrib_vect})
        return [trace,events,attrib_org,attrib_vect,k,delta,regul,model_dict]
    else:
        return [None,None,None,None,k,delta,regul,model_dict]


def det_gausstrace(N,loc,amp,sigma,alpha,phi):
    """
    Generates seismic trace  with certain features(by Gaussian waveforms)
    Seismic source waveform is positive i.e.                   __  
                                            (@ zero has max \_/  \_/ )
    Input:
        N       : Number of samples
        loc     : Location of events (assumed as only nonzero points)
        amp     : Maximum amplitude for spikes, amplitude range in [-amp,+amp]
        sigma   : scale of events
        alpha   : Singularity order coressponding to each event
        phi     : Phase rotation component: [0,2)
        """

    
    if N%2 == 1:
	N+=1

    loc = __np.array([loc]).ravel()
    amp = __np.array([amp]).ravel()
    sigma = __np.array([sigma]).ravel()
    alpha = __np.array([alpha]).ravel()
    phi = __np.array([phi]).ravel()%2
    
    k=loc.shape[0]

    delta=N
    attrib_org=__np.array([loc,sigma,alpha,phi]).transpose()

    if k>1:
        perm=__np.argsort(attrib_org[:,0],axis=0)
        attrib_org=attrib_org[perm,:]
        loc=loc[perm]
        delta=__np.diff(loc)
        if k>2:
            delta=min(delta)
        if amp.shape[0]==1:
            amp=amp *__np.random.rand(k)
        else:
            amp=amp[perm]

            
    attrib_vect=_Null * __np.ones((4,N))
    attrib_vect[0,loc]=amp
    attrib_vect[1,loc]=attrib_org[:,1]
    attrib_vect[2,loc]=attrib_org[:,2]
    attrib_vect[3,loc]=attrib_org[:,3]
    
    trace=__np.zeros(N)
    events=__np.zeros((k,N))
    for i in range(k):
        # Number of generated events
        f=__manf.gmanf(N)
        f.update(attrib_org[i,:])
        events[i,:]=amp[i]*f.data(tnorm=1)[0] /max(abs(f.data(tnorm=1)[0] )) 
        # for Amp. same as spikes
        trace = trace+events[i,:]
    
    regul=None
    model_dict={'trace':trace,'events':events,'attrib_org':attrib_org,
                'attrib_vect':attrib_vect, 'k_act':k,'delta_act':delta,
                'regul':regul}
    return [trace,events,attrib_org,attrib_vect,k,delta,regul,model_dict]

# Add documentation of outputs to functions
det_trace.__doc__ += _doc_out
rand_trace.__doc__ += _doc_out
det_gausstrace.__doc__ += _doc_out
rand_gausstrace.__doc__ += _doc_out
#========================================
#  Reflectivity-Source Convolution Func.
#========================================
def src_conv(reflect,src=_df_src,scl=_df_scl,pad=15):
    """
    Generates seismic trace by convolving source signature 
    with reflectivity model
    Warning: Make distortion at edges (because of cicular effect of FFT), 
             No Possible Fix
 
    Input:
        reflect : reflectivity model (given by splinetrain function) in time
        src     : a wavelet class object to define source signature 
        scl     : scale of source signature to be used

    Output:
        trace   : Seismic Trace after source signature convolution
    """
    N=len(reflect)
    len_pad=round(N/pad)
    pad1=reflect[0]*__np.ones(len_pad)
    pad2=reflect[-1]*__np.ones(len_pad)
    r_pad=__np.concatenate((pad1,reflect,pad2))
    N2=len(r_pad)
    wav= src.wavf(N2,s=scl,tnorm=1,domain='freq',mode='complex')[0]
    trace_pad=__np.real( __ft.ifft(wav*__ft.fft(r_pad)) )
    trace=trace_pad[len_pad:N2-len_pad]
    #     wav= src.wavf(N,s=scl,tnorm=1,domain='time',mode='complex')[0]
    #     wav= __ft.ifftshift(wav)
    #     trace=__np.convolve(reflect,wav,mode='same')
    return trace/__misc.norm(trace)
