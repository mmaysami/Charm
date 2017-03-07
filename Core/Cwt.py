"""

      Continous Wavelet Tranform
    ==============================

    Functions
    ----------
    extrema   -- Find extrema of a 2-D matrix of wavelet coefficient
    mmlimage  -- Plot CWT coefficient overlayed with MML & their maxima point
    gauswavf  -- Gaussian Wavelet
    mexhwavf  -- Mexican Hat Wavelet
    morlwavf  -- Morlet Wavelet

    Class
    ----------
    wavelet 



AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""
#__all__=['extrema','mmlimage','gauswavf','mexhwavf','morlwavf','wavelet']
import numpy as _wavelet__np
import numpy as __np
import numpy.fft as _wavelet__ft
import numpy.fft as __ft
import Charm.Core.Misc as __misc
import Charm.Core.Misc as _wavelet__misc

from  platform import uname  as __uname
from   sys     import stderr as _wavelet__stderr
from   Charm import _eps,_df_xlabel
try:
    import pylab as __py
    import pylab as _wavelet__py 
except ImportError:
    from Charm import __Disp_Err
    __py=None


_df_scl=2**3
_df_sclrng=2**__np.linspace(-1,4,50)
_df_tnorm=1
_df_tflip=0
_df_rad=3
_df_acc=1e-4
_df_con=__np.sqrt(10)

_tsize=__py.rcParams['axes.titlesize']
_lsize=__py.rcParams['axes.labelsize']
#=====================================
#         Define Functions
#=====================================

def extrema(data,rad=_df_rad,acc=_df_acc):
    """
    find exterema of wavelet transform output (2-D signal)
    
    Input:
        data     : 2-D array  whose exterma is going to be found
	rad      : topological radius
	acc      : order of accuracy

    Output:
        Extrema points of data
    """
    nt = data.shape[0]
    if len(data.shape)==1:
        data=__misc.vector(data,'col')
    const,ind,flag=__misc.max2d(abs(data))
    ext = __np.zeros(data.shape,dtype=int)
    Ldif = data[rad:nt-rad,:] - data[0:nt-2*rad,:]
    Rdif = data[rad:nt-rad,:] - data[2*rad:nt,:]
    ext[rad:nt-rad,:] = __np.sign(Ldif)+__np.sign(Rdif) /2
    ext[rad:nt-rad,:] = ext[rad:nt-rad,:] * (abs(Ldif*data[rad:nt-rad,:]) > acc*const)
    return ext 


def mmlimage(max_points,mml,acw,scatter=1,cbar=0,label=1):
    """
    find exterema of wavelet transform output (2-D signal)
    
    Input:
        maxpoints : Array of location and scale index of maximum points 
                    along MMLs (Given by detect function)
	mml       : Array of modulus maxima lines
	acw       : Absolute value of CWT coefficient
        scales    : Scales used for CWT
        scatter   : Show maximum points as a scatter plot
        cbar      : Show colorbar
        label     : Show Title and axis labels

    Output:
        Plot image of MMLs on top of CWT coefficient and show maximums 
    """

    if __py == None:
        print >>__stderr, __Disp_Err
        pass
        
    data=__py.rot90(acw-mml*(mml>0)*acw,k=1)
    if __uname()[0].lower()[0:5]!='linux':
        data=__py.flipud(data)   
    rcdict0=__py.rcParams.copy()
    __py.rc('figure',figsize=(10,3))
    __py.figure()
    __py.imshow(data,aspect=4)
    if cbar:
        __py.colorbar()
    if scatter:
        __py.scatter(max_points[:,0],max_points[:,1],
                   s=25,c='w',hold=True)
    if label:
        __py.title('Modulus of CWT & MMLs',size=_tsize)
        __py.ylabel('Scale index',size=_lsize)
        __py.xlabel(_df_xlabel,size=_lsize)

    __py.rcParams.update(rcdict0)
    pass


def gauswavf(N,s=_df_scl,p=2,wtype='complex',tflip=_df_tflip,tnorm=_df_tnorm,
             domain=None,mode=None):
    """
    Zero-Phase Gaussian wavelet function (-1)**n d**n/dt**n (Theta)
    
    Input:
        N      : Number of samples
	s      : scale
	p      : order of derivative for Gaussian wavelet
	wtype  : type of wavelet  ['complex']/'real'
	tflip  : returns f(-t) if 1
        tnorm  : L2 normalized output if True
	domain : Output data domain ['time']/'freq'
	mode   : type of data as complex,real-imaginary,and absolute-phase ['complex']/'re-im'/'ab-ph'

    Output:
        Gaussian wavelet function in proper form 
    """
    w = 2*__np.pi/N * __np.arange(_eps,N/2+1)
    if tflip==1:
	w= -w
    psi = (w)**p * __np.exp(-(w*s)**2)
    psi=__misc.fillfreq(psi,tnorm=tnorm,stype=wtype,axis=0)[0:N/2+1]
    # Normalization for full frequnecy range is required if needed!
    return __misc.outtype(psi,domain=domain,mode=mode,tnorm=0,stype=wtype)


def mexhwavf(N,s=_df_scl,wtype='complex',tnorm=_df_tnorm,domain=None,mode=None):
    """
    Mexican Hat wavelet function (Symmetric)

    Input:
        N      : Number of samples
	s      : scale
	wtype  : type of wavelet  ['complex']/'real'
        tnorm  : L2 normalized output if True
	domain : Input data domain ['time']/'freq'
	mode   : type of data as complex,real-imaginary,and absolute-phase ['complex']/'re-im'/'ab-ph'
	
    Output:
        Mexican Hat wavelet function in proper form 
    """
    t   = __np.linspace(round(-N/2),0,N/2+1) # Effective Width : [-5,5]
    t   = __np.concatenate((t,-t[N/2-1:0:-1]))
    psi = (2/__np.sqrt(3)*__np.pi**-0.25) * (1-(t/s)**2)*__np.exp(-(t/s)**2/2)
    psi = __ft.fft(__ft.fftshift(psi))[0:N/2+1]
    psi=__misc.fillfreq(psi,tnorm=tnorm,stype=wtype,axis=0)[0:N/2+1]	 
    return __misc.outtype(psi,domain,mode=mode,tnorm=0,stype=wtype)


def morlwavf(N,s=_df_scl,w0=1.5*__np.pi,wtype='complex',tflip=_df_tflip,
             tnorm=_df_tnorm,domain=None,mode=None):
    """
    Morlet wavelet function (Symmetric)

    Input:
        N      : Number of samples
	s      : scale
	w0     : frequency shift parameter for Morlet wavelet
	wtype  : type of wavelet  ['complex']/'real'
	tflip  : returns f(-t) if 1
        tnorm  : L2 normalized output if True
	domain : Input data domain ['time']/'freq'
	mode   : type of data as complex,real-imaginary,and absolute-phase ['complex']/'re-im'/'ab-ph'
	
    Output:
        Morlet wavelet function in proper form 
    """
    w = 2*__np.pi/N * __np.arange(_eps,N/2+1)
    if tflip==1:
	w= -w
    psi =__np.sqrt(2*__np.pi)*__np.exp(-(w*s-w0)**2/2)
    psi=__misc.fillfreq(psi,tnorm=tnorm,stype=wtype,axis=0)[0:N/2+1]
    return __misc.outtype(psi,domain=domain,mode=mode,tnorm=0,stype=wtype)


#=====================================
#  Define Class
#=====================================
class wavelet(object):
    """
    wavelet class for doing continuous wavelet transfrom
    
      CWT class
    ===============================
    'gaus'  : Gaussian
    'mexh'  : Mexican Hat 
    'morl'  : Morlet

    wavf    : Returns wavelet function values in certain domains
    show    : shows wavelet function values in certain domains
    cwt     : Does Continuous Wavelet Transform on data
    showcwt : Show the result of CWT on data in proper form
    mml     : Computes the modulus maxima lines for CWT 
    mmle    : Extended version of mml method
    """  
    def __init__(self,name='gaus',wtype='complex',param={'p':2}):
	"""
	Initialize wavelet object

	Input:
	    name  : Name of wavelet
	    wtype : Type of wavelet ['complex']/'real'
	    param : Parameters specfic to this wavelet
	"""
	self.name=name.lower()
	self.param=param
	self.wtype=wtype
	    
    def __param__(self):
	"""
	Print parameters for wavelet object 
	"""
        print >>__stderr, 'Wavelet Object Parameters'
        print >>__stderr, '=========================' 
        print >>__stderr, 'Name  :',self.name
        print >>__stderr, 'Type  :',self.wtype
        print >>__stderr, 'Param :',self.param
        dict={'name':self.name,'wtype':self.wtype,'param':self.param}
        return dict

    def wavf(self,N,s=_df_scl,tflip=_df_tflip,tnorm=_df_tnorm,domain=None,
             mode=None,*args,**kwargs):
	"""
	Wavelet function
	
	Input:
            N     : Number of samples
	    s     : scale
            tflip : returns f(-t) if 1
            tnorm : L2 normalized output if True
	    domain: output data domain None (half frequency range)/['time']/'freq'
	    mode   : type of output: for frequency domain as complex,real-imaginary,and absolute-phase ['complex']/'re-im'/'ab-ph'
	
        Output:
                 Wavelet function in proper form 
	"""
	wavefun=eval(self.name+'wavf') # returns a pointer to function
	return wavefun(N,s=s,wtype=self.wtype,tflip=tflip,tnorm=tnorm,
                       domain=domain,mode=mode,*args,**self.param)


    def show(self,N,s=_df_scl,domain='time',mode=None,tshift=1,*args,**kwargs):
	"""
	Show Normalized & Centeralized Wavelet function
	
	Input:
            N      : Number of samples
	    s      : scale
	    domain : output data domain None 
                    (half frequency range)/['time']/'freq'
	    mode   : type of output: for frequency domain as complex,
                     real-imaginary,and absolute-phase ['complex']/'re-im'
                     /'ab-ph' and for time domain no shift 'orig', and shifted 
                     to center ['shift']
	
        Output:
                 Wavelet function in proper form and Generates proper figure 
	"""
        if __py == None:
            print >>__stderr, __Disp_Err
            pass
	wf=self.wavf(N,s=s,tnorm=1,domain=domain,mode=mode)
	__misc.plotdata(wf,domain=domain,mode=mode,stype=self.wtype,
                        tshift=tshift,label= True,*args,**kwargs)
	__py.title('Wavelet Function in '+domain.lower(),size=_tsize)
	#__py.show()
        pass

	
    def cwt(self,data,s=_df_sclrng,out='complex'):
	"""
	Computes the Continuous Wavelet Transform using 

	    F(x,sigma)  =  WT{f,psi(sigma^p)}
	       (sigma**p) /{x**p}(f*psi(sigma^p))

	Input:
	    data: the input data to be transformed (Assume no shift is needed)
	    s	: vector containing the scales
	    p	: the number of vanishing moments for Gaussian 
	    out : output type ['abs']/'complex'

	Output:	
	    The wavelet transform coefficient (Column index start from 
            fine scale [0] to coarse scale [-1])
	"""
	s = s.ravel()
	data=__np.array(data).ravel()
	Nd=len(data)
	Nd += Nd%2

        psi = self.wavf(Nd,s=__np.max(s),tnorm=1,domain='time',
                        mode='ab-ph',**self.param)[0]
        Np =  len(psi >= 10*__np.min(psi))
        Np += Np%2

        data=__np.concatenate([data,__np.zeros(Np)])
        N = Nd+ Np

	fdata = __ft.fft((data),N)

	wt = __np.zeros( (N,len(s)) ,dtype=__np.complex)
	for si in __np.arange(0,len(s)):
	    # psi data for half frequency range
	    psi = self.wavf(N,s=s[si],tflip=1,tnorm=1,domain='freq',
                            mode='complex',**self.param)[0]
	    wt[:,si] = fdata *psi
	   
	wt = __ft.ifft(wt,axis=0)           #fft along time x-axis=0


	#if self.name=='gaus':
       	# wt=__np.flipud(wt)
	if out == 'abs':
	    wt= abs(wt)
	elif out != 'complex':
	    raise ValueError,'Unsupported type for output. Please use "abs" or "complex"'

	return __np.array(wt[0:Nd,:])



    def showcwt(self,data,s=_df_sclrng,interpol='nearest',out=None,
                *args,**kwargs):
	"""
	Show real Wavelet Transformed 1D-data in proper form

	Input:
	    data     : the input data to be transformed 
                       (Assume no shift is needed)
	    s	     : vector containing the scales
	    interpol : interpolation of image 
                       (i.e. 'nearest', 'bilinear','gaussian',...)
	    out      : output type ['abs']/'complex'/None (for no output)

	Output:	
	    The Wavelet Transformed data and Generate proper figure

	"""
        if __py == None:
            print >>__stderr, __Disp_Err
            pass
	if out!=None:
	    cw = self.cwt(data,s=s,out=out)
	else:
	    cw = self.cwt(data,s=s,out='abs')
	ax1=__py.subplot(121)
        # Hide Tick Labels
        ylab=ax1.get_yticklabels()
        __py.setp(ylab,visible=0)

 	t=__np.arange(len(data))
	__py.plot(data,__np.flipud(t),lw=1)
	if min(data)>=0 :
	    __py.xlim(xmin=-1)
	if max(data)<=0 :
	    __py.xlim(xmax=1)

	__py.title('Data',size=_tsize)
	__py.xlabel('Amplitude',size=_lsize)
	__py.ylabel('Time',size=_lsize)

	ax2=__py.subplot(122)
        # Hide Tick Labels
        ylab=ax2.get_yticklabels()
        xlab=ax2.get_xticklabels()
        __py.setp(xlab,visible=0)
        __py.setp(ylab,visible=0)

	__py.imshow(abs(cw),interpolation=interpol,origin='upper',
                    *args,**kwargs)
	__py.title(self.wtype+' CWT of Data',size=_tsize)
	__py.xlabel('Scales',size=_lsize)
	__py.ylabel('Locations',size=_lsize)
	__py.colorbar()

	#__py.show()
	if out!=None:
	    return cw

    def mml(self,data,s=_df_sclrng,con=_df_con,out='complex',negline=0,
            rad=_df_rad,acc=_df_acc,*args,**kwargs):
	"""
	Computes the modulus maxima lines (Starting from finest scale) ) for 
        Continuous Wavelet Transform 

	Input:
	    data: the input data to be transformed (Assume no shift is needed)
	    s	    : Vector containing the scales
            con     : Connectivity factor (max. allowed distance between 
                      connected points)
	    out     : Output type ['abs']/'complex'
            negline : Computes negative mml lines
            rad     : topological radius (For finding extrema lines of cwt)
            acc     : order of accuracy (For finding extrema lines of cwt)

	Output:	
	    mml   : Modulus maxima lines (MML)
            cw    : Continuous wavelet transform (CWT)
            lnpos : Positive MML as list (for separately access to each line)
            ext   : Extrema points of CWT
	"""
        # NOTE : cw[:, fine[0] --> coarse [-1] ]
        # Form cwt and find extrema of coeffiicients
        cw = self.cwt(data,s=s,out=out)
        ext = extrema(abs(cw),rad=rad,acc=acc)

        con=con**2+1e-15
        if __np.isscalar(con):
            con = con*__np.ones((ext.shape[1]))

        # Find Pos. & Neg. location for the finest scale
        pos0=__np.nonzero(ext[:,0] >0)[0]
        neg0=__np.nonzero(ext[:,0] <0)[0]
        # Number of Pos. & Neg. lines in fine scale
        nbposl=len(pos0)
        nbnegl=len(neg0)

        # Structure for storing MML
        lnpos=list()
        lnneg=list()

        # Add start points to storage
        for i in range(nbposl):
            lnpos.append([[pos0[i],0]])

        for i in range(nbnegl):
            lnneg.append([[neg0[i],0]])

        # End points of lines at current scales
        endpp=__np.zeros((nbposl,2))
        endpn=__np.zeros((nbnegl,2))

        # Move along scales and find points to be added to lines                
        for k in __np.arange(1,ext.shape[1]):
            # Positive locations for current scale
            pos=__np.nonzero(ext[:,k] >0)[0]
            nbpos=len(pos)
            # If any pos. point is detected
            if nbpos:
                for ni in range(nbposl):
                    # Load end points for distance comparison
                    endpp[ni,:] = lnpos[ni][-1]
                    if __np.isscalar(lnpos[ni][-1]):
                        endpp[ni,:] = lnpos[ni]

                for ni in range(nbpos):                                       
                    # Find closest end point to pos[ni] to fit in one Pos.Line
                    newp = [pos[ni],k]
                    dist = __np.sum( (endpp-newp*__np.ones((nbposl,2)))**2
                                     ,axis=-1 )
                    # Check to hold connectivity condition
                    dist = dist* (dist < con[k] )
                    
                    #Index of nnz elements of dist
                    nnz=dist.nonzero()[0]
                    if len(nnz) != 0:
                        lm=__np.argmin(dist[nnz])
                        # Index of closest point among end points
                        lm=nnz[lm]

                        # Load End point of the line
                        prevp1= lnpos[lm][-1]
                        if __np.isscalar(lnpos[lm][-1]):
                            prevp1= lnpos[lm]
                        # Add point if the line has no point with this scale
                        if prevp1[1] != k:
                            lnpos[lm].append(newp)
                        # o.w. check which one is closer
                        else:
                            prevp2 = lnpos[lm][-2]
                            dist_prv = __np.sum((__np.array(prevp1) - 
                                                 __np.array(prevp2))**2)
                            dist_new = __np.sum((__np.array(newp) - 
                                                 __np.array(prevp2))**2)
                            # Replace end point if new point is closer
                            if dist_new < dist_prv :
                                lnpos[lm][-1]=newp


            # Adding Negative lines
            # If any neg. point is detected
            if negline:
                neg=__np.nonzero(ext[:,k] <0)[0]
                nbneg=len(neg)
                if nbneg:
                    for ni in range(nbnegl):
                        # Load end points for distance comparison 
                        endpn[ni,:] = lnneg[ni][-1]
                        if __np.isscalar(lnneg[ni][-1]):
                            endpn[ni,:] = lnneg[ni]

                    for ni in range(nbneg):                                       
                        # Find closes end point to neg[ni] to fit one Neg.Line
                        newp = [neg[ni],k]*__np.ones((nbnegl,2))
                        dist = __np.sum( (endpn-newp)**2,axis=-1 )
                        dist = dist* (dist < con[k] )
                        newp=[neg[ni],k]
                        nnz=dist.nonzero()[0]
                        if len(nnz) != 0:
                            lm=__np.argmin(dist[nnz])
                            lm=nnz[lm]
                            # Load End point of the line
                            prevp1= lnneg[lm][-1]
                            if __np.isscalar(lnneg[lm][-1]):
                                prevp1= lnneg[lm]

                            # Add point if the line has no point with this scale
                            if prevp1[1] != k:
                                lnneg[lm].append(newp)
                            # o.w. check which one is closer
                            else:
                                prevp2 = lnneg[lm][-2]
                                dist_prv = __np.sum((__np.array(prevp1) - 
                                                     __np.array(prevp2))**2)
                                dist_new = __np.sum((__np.array(newp) - 
                                                     __np.array(prevp2))**2)
                                # Replace end point if new point is closer
                                if dist_new < dist_prv :
                                    lnneg[lm][-1]=newp

        # Forming MML matrix from detected lines
        # Put 1 for Pos. lines and -1 for Neg. lines
        mml=__np.zeros(cw.shape,dtype=int)

        if nbposl==0:
            print >>__stderr, " ==> Warning : No event was detected! "
            return ( mml, cw, lnpos, ext )

        lengthpos=__np.zeros(nbposl)
        for i in range(nbposl):
            lengthpos[i]=len(lnpos[i])
            lnpos[i]=__np.array(lnpos[i])

        # Divide start points to chunk of adjacent ones
        stp=1+(__np.diff(pos0)>1).nonzero()[0]
        #print >>__stderr, "STP",stp
        stp=__np.concatenate(([0],stp,[len(pos0)]))
        endline=0
        shft=0
            
        newlnpos=list()     
        for k in range(len(stp)-1):
            # Find longest line among adjacent lines
            ind=__np.argmax(lengthpos[stp[k]:stp[k+1]])+stp[k]
            pline=lnpos[ind]
            # Shift of location to get central MML
            shft=(pos0[stp[k]]+pos0[stp[k+1]-1]) /2 -pline[0,0]
            # Form a new list of lines with only central lines 
            # (no repeat of adjacent lines)
            newlnpos.append(__np.array([pline[:,0]+shft,pline[:,1]]).transpose())

            # Add adjacent lines (make a thick MML with a width of rad)
            for r in range(-(rad/2),1+rad/2):
                mml[pline[:,0]+shft+r,pline[:,1]]=1

        # ======================================
        # Add negative lines to mml if required
        if negline:
            lengthneg=__np.zeros(nbnegl)       
            for i in range(nbnegl):
                lengthneg[i]=len(lnneg[i])
                lnneg[i]=__np.array(lnneg[i])

            # Divide start points to chunk of adjacent ones
            stp=1+(__np.diff(neg0)>1).nonzero()[0]
            stp=__np.concatenate(([0],stp,[len(neg0)]))
            endline=0
            shft=0

            for k in range(len(stp)-1):
                ind=__np.argmax(lengthneg[stp[k]:stp[k+1]])+stp[k]
                nline=lnneg[ind]
                # Shift of location of central MML
                shft=(neg0[stp[k]]+neg0[stp[k+1]-1]) /2 -nline[0,0]
            
                # Add adjacent lines (make a thick MML with a width of rad)
                for r in range(-(rad/2),1+rad/2):
                    mml[nline[:,0]+shft+r,nline[:,1]]=-1


        return ( mml, cw, newlnpos, ext )





    def mmle(self,data,s=_df_sclrng,con=_df_con,out='complex',negline=0,
             rad=_df_rad,acc=_df_acc,search=1./8,mindist=5,*args,**kwargs):
	"""
	Extended version of mml method
        Computes the modulus maxima lines (Starting from not only finest 
        scale but in a fine scale region) for Continuous Wavelet Transform 

	Input:
	    data: the input data to be transformed (Assume no shift is needed)
	    s	    : Vector containing the scales
            con     : Connectivity factor (max. allowed distance between 
                      connected points)
	    out     : Output type ['abs']/'complex'
            negline : Computes negative mml lines
            rad     : Topological radius (For finding extrema lines of cwt)
            acc     : Order of accuracy (For finding extrema lines of cwt)
            search  : Portion of scales to be searched for start point of MMLs
            mindist : Min. Distance of 2 MMLs to be deteced as separate lines

	Output:	
	    mml   : Modulus maxima lines (MML)
            cw    : Continuous wavelet transform (CWT)
            lnpos : Positive MML as list (To separately access to each line)
            ext   : Extrema points of CWT
	"""

        # Find normal lines starting from finest scale
        ( mml, cw, newlnpos, ext ) = self.mml(data,s=s,con=con,out=out,
                                              negline=negline,rad=rad,acc=acc)

        stscl=search*len(s)
        con=2*con**2+1e15
        # ================================
        # Check Positive lines if required
        flagp=0
        ext2=ext.copy()
        # Difference of Extrema & MML to catch other (minor) Lines
        ext2=(ext2-mml)
        # Check for new Postive MMLs (Patch)
        pos=(ext2>0).nonzero() 
        pos2=__np.array([__np.concatenate([[0],pos[0]])
                         ,__np.concatenate([[0],pos[1]])])
        # Separte points to lines with forming edge
        edgep=__np.logical_or((abs(__np.diff(pos2[1]))>con),
                         (abs(__np.diff(pos2[0]))>con)).nonzero()[0]

        lnpos=list()
        pos0=list()
        # Check All founded Positive MML
        for k in range(len(edgep)-1):
            # Take the line and sort it w.r.t scales
            pline=(pos[0][edgep[k]:edgep[k+1]],pos[1][edgep[k]:edgep[k+1]])
            perm=__np.argsort(pline[1])
            pline=__np.array((pline[0][perm],pline[1][perm])).transpose()
            # Check if the start point of line lies inside fine scale region
            if pline[0,1] < stscl:
                # Check if there is any line in neighborhood
                cond = mml[ pline[0,0]+range(-mindist,mindist+1),pline[0,1]]
                if not (cond>0).any():
                    # ======= Found a minor Pos. MML ========
                    flagp=1
                    lnpos.append(pline)
                    # Form start locations of new lines
                    pos0.append(pline[0,0])

        if flagp==1:
            # ======= Add Pos. MML(s) to output variables ========
            print >>__stderr, "NOTE : Positive minor MML(s) was detected "
            nbposl=len(lnpos)
            lengthpos=__np.zeros(nbposl)
            for i in range(nbposl):
                lengthpos[i]=len(lnpos[i])
                lnpos[i]=__np.array(lnpos[i])

            if nbposl !=0:
                # Divide start points to chunk of adjacent ones
                stp=1+(__np.diff(pos0)>1).nonzero()[0]
                stp=__np.concatenate(([0],stp,[len(pos0)]))

                for k in range(len(stp)-1):
                    ind=__np.argmax(lengthpos[stp[k]:stp[k+1]])+stp[k]
                    pline=lnpos[ind]
                    # Shift of location to get central central MML
                    shft=(pos0[stp[k]]+pos0[stp[k+1]-1]) /2 -pline[0,0]
                    # Form a new list of lines with only central lines
                    # (no repeat of adjacent lines)
                    newlnpos.append(__np.array([pline[:,0]+shft,
                                                pline[:,1]]).transpose())
                    # Add central line and adjacent lines 
                    # (make a thick MML with a width of rad)
                    for r in range(-(rad/2),1+rad/2):
                        mml[pline[:,0]+shft+r,pline[:,1]]=1

        # ================================
        # Check Negative lines if required
        if negline==1:
            # Check for new Negative MMLs (Patch)
            flagn=0
            lnneg=list()
            neg=(ext2<0).nonzero() 
            neg2=__np.array([__np.concatenate([[0],neg[0]]),
                             __np.concatenate([[0],neg[1]])])
            # Separte points to lines with forming edge
            edgen=__np.logical_or((__np.diff(neg2[1])>=con) ,
                                  (__np.diff(neg2[0])>=con)).nonzero()[0]

            lnneg=list()
            neg0=list()
            # Check for every founded Negative MML
            for k in range(len(edgen)-1):
                nline=(neg[0][edgen[k]:edgen[k+1]],neg[1][edgen[k]:edgen[k+1]])
                perm=__np.argsort(nline[1])
                nline=__np.array((nline[0][perm],nline[1][perm])).transpose()
                # Check if the start point of line lies inside fine scale region
                if nline[0,1] < stscl:
                    # Check if there is any line in neighborhood
                    cond = mml[ nline[0,0]+range(-mindist,mindist+1),nline[0,1]]
                    if not (cond>0).any():
                        # ======= Found a minor Pos. MML ========
                        flagn=1
                        lnneg.append(nline)
                        # Form start locations of new lines
                        neg0.append(nline[0,0])

            if flagn==1:
                # ======= Add Neg. MML(s) to output variables ========
                print >>__stderr, "NOTE : Negative minor MML(s) was detected "
                #                 mml[nline[:,0],nline[:,1]]=-1

                nbnegl=len(lnneg)
                lengthneg=__np.zeros(nbnegl)
                for i in range(nbnegl):
                    lengthneg[i]=len(lnneg[i])
                    lnneg[i]=__np.array(lnneg[i])
                    
                # Divide start points to chunk of adjacent ones
                stp=1+(__np.diff(neg0)>1).nonzero()[0]
                stp=__np.concatenate(([0],stp,[len(neg0)]))

                for k in range(len(stp)-1):
                    ind=__np.argmax(lengthneg[stp[k]:stp[k+1]])+stp[k]
                    nline=lnneg[ind]
                    # Shift of location to get central central MML
                    shft=(neg0[stp[k]]+neg0[stp[k+1]-1]) /2 -nline[0,0]
                    # Add central line and adjacent lines 
                    # (make thick MML with a width of rad)
                    for r in range(-(rad/2),1+rad/2):
                        mml[nline[:,0]+shft+r,nline[:,1]]=-1


        return ( mml, cw, newlnpos, ext)
