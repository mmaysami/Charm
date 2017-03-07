"""

           Manifold
    ======================
    Appearance Manifolds for fractionaly differentiated Gaussians

    Functions
    ----------
    movavg  
    smooth
    wgauss

    Class
    ----------
    gmanf  --  Gaussian Manifold with different attribute



AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""
# __all__=['movavg','smooth','wgauss','gmanf']



import copy  as _gmanf__copy
import numpy as _gmanf__np
import numpy as __np
import numpy.fft as _gmanf__ft
import numpy.fft as __ft
import Charm.Core.Misc as _gmanf__misc
import Charm.Core.Misc as __misc
from   sys     import stderr  as _gmanf__stderr
from   Charm   import _eps

try:
    import pylab as _gmanf__py 
    import pylab as __py
except ImportError:
    from Charm import __Disp_Err
    __py=_gmanf__py=None



_df_tshift=0
_df_tnorm=0
_df_dnorm=1

_df_stype='real'
_df_dtype='real'

_df_d3const=-1j

#=====================================
#  Define Class
#=====================================
class gmanf(object):
    """
    Gaussian Appearance Manifolds class
    ==================================
    help        : Describe how to create an instance of manifold
    update      : Update parameter values for manifold instance
    data        : Manifold values in certain domains
    param       : A summary of parameter values
    d0,d1,d2,d3 : Manifold partial derivatives w.r.t. is parameters
    grad        : Gradient of manifold 
    plot,show   : A figure handle or show the figure of manifold values
    smooth      : Uses Gaussian smoother on manifold data and 
                  may plot the results 
    norm        : Norm of Guassian waveform vector
    support     : Approximate support of data vector of manifold
    eval        : Attribute of the manifold oobject 
    copy        : Make a copy of manifold instance
    
    See gmanf.__init__  for more info on setting manifold instances.

    """
    def __init__(self,n=1024,tau=512,sigma=10,alpha=0,phi=None,s=0,fi=0):
	    """
	    Construct Smoothed Gaussian Pulse given all Parameteres  
	    Works for odd length as well as even length of signal

	    INPUT:
	        N     : Length of signal (Even/Odd?) [1024] 
		Param : [tau,sigma,alpha,phi]: values of Parameters [512,10,0,0]
                        if phi=-alpha/2 there is no phase shift 
		s     : Smoothing scale [0]
		fi    : Fractional integration order [0]

            OUTPUT:
		Manifold Object 
            """
	    if n%2==1:
		    n += 1
	    self.n = n
	    self.tau   = tau
	    self.sigma = sigma
	    self.alpha = alpha
            if phi==None:
                self.phi   = -0.5*alpha
            else:
                self.phi   = phi
	    self.s    = s
	    self.fi   = fi

    #def __call__(self,domain='time',mode='complex',*args,**kargs):
    #return self.data(domain='time',mode='complex',*args,**kargs)

    def param(self,*args,**kargs):
            """
            Returns all the parameters of manifold objects
            """
            dict={'n':self.n,'sigma':self.sigma,'tau':self.tau,
                  'alpha':self.alpha,'phi':self.phi,'s':self.s,'fi':self.fi}

            print >>__stderr,'\nGaussian Manifold Object Parameters'
            print >>__stderr,'===================================\n' 
            print >>__stderr,'Length :',self.n
            print >>__stderr,'Tau    :',self.tau
            print >>__stderr,'Sigma  :',self.sigma
            print >>__stderr,'Alpha  :',self.alpha
            print >>__stderr,'Phi    :',self.phi
            print >>__stderr,'Smoothing factor            :',self.s
            print >>__stderr,'Fraction Integration factor :',self.fi,'\n'
            pass

    # Makes self.param as a call to this function  
    #     param = property(__param__)


    def help(self):
            print >>__stderr, """
            Construct Smoothed Gaussian Pulse given all Parameteres"  
            Works for odd length as well as even length of signal\n"

            INPUT:
                N     : Length of signal (Even/Odd?) [1024] 
                Param : [tau,sigma,alpha,phi]: values of Parameters [512,10,0,0]
                        if phi=-alpha/2 there is no phase shift 
                s     : Smoothing scale [0]
                fi    : Fractional integration order [0]

            OUTPUT:
                Manifold Object 
                """


    def eval(self,att='tau'):
	    """
	    Returns attribute of the manifold oobject which is 
            pointed by string.

	    INPUT:
		str : String that determine attribute of manifold objecte [None]

            OUTPUT:
		Manifold Object's attribute
	    """
            return eval('self.'+str(att))


        
    def copy(self,s=None,fi=None):
	    """
	    Returns a copy of the manifold object with requested 
            change in parameters.

	    INPUT:
                if any input variable is not None, then the value 
		s     : Smoothing scale [None]
		fi    : Fractional integration scale [None]

            OUTPUT:
		Manifold Object copy with any change in parameters if needed
	    """
	    temp=__copy.copy(self)
            if s != None:
                temp.s=s
            if fi != None:
                temp.fi=fi

	    return temp


    def update(self,param):
            """
            Updates values of parameters for Gaussian manifold 
            with those of param

	    INPUT:
		param  : Tuple of values for (tau,sigma,alpha,phi)

            OUTPUT:
		Manifold Object copy with any change in parameters if needed
                """
            self.tau=param[0]
            self.sigma=param[1]
            self.alpha=param[2]
            self.phi=param[3]
            pass


        
    def support(self,eps=1e-5):
	"""
	Returns approximate support of Gaussian waveform in time domain
	"""
        tau=self.tau
        self.tau=self.n/2
	d=self.data(domain='time',mode='complex',tnorm=1)[0]
	cond=(abs(__np.gradient(d))>=eps).nonzero()[0]
	ind=range(len(d))
	mask=__np.zeros(len(d))
	mask[__np.logical_and(ind>=cond[0],ind<=cond[-1])]=1
        self.tau=tau
	return cond[-1]-cond[0]+1



    def data(self,domain='time',mode='complex',tnorm=_df_tnorm,tshift=_df_tshift):
	    """   
	    Returns data values of the manifold object in [time]/freq domain. 
            In frequency domain, returned data could be as [complex], 
            separated real-imaginary parts, or  absolute-phase values.
            Note1: Frequency range is similar to FFT e.g. :[0,2*pi] or [0,+f,-f]
	    Note2: Time domain data is shifted inorder to have 0 at center 
                   of X-axis

	    INPUT:
	        domain : domain of output data ['time']/'freq'/None 
                         (return input data without change)
	        mode   : type of output as complex,real-imaginary,and 
                         absolute-phase ['complex']/'re-im'/'ab-ph'
                tnorm  : L2 time normalized output if True
	        tshift : Time shift flag for time domain (centered in Time if 1)

            OUTPUT:
		Tuple of data parts.If data should have only one part then 
                the second part is None (Data is zero-mean)
	    """
            DC=0
	    N=self.n
	    (sigma,tau,alpha,phi)=(self.sigma,self.tau,self.alpha,self.phi)
	    (s,fi)=(self.s,self.fi)
            sigma=__np.sqrt(sigma**2+s**2)
            alpha+=fi

	    w=2*__np.pi/N * __np.arange(_eps,N/2+1)

            #Equivalent to: (-1j*w)**(-0.5*alpha-phi) * (1j*w)**(-0.5*alpha+phi)
	    F= (w)**(-alpha) *__np.exp(1j*__np.pi*phi) *__np.exp(-(w*sigma)**2 /2)* __np.exp(-1j*w* tau) 
            if alpha > 0:
                F[0] = DC

	    return __misc.outtype(F,domain,mode,tshift,tnorm=tnorm
                                  ,stype=_df_stype)


    def norm(self,domain='time',ord=2):
        """
        Returns norm of manifold values

	    INPUT:
	        domain : Domain of data to calculate norm ['time']/'freq'
                ord    : Norm order

           OUTPUT:
                Norm of manifold                
        """
        d=self.data(domain=domain,mode='complex')[0]
        return __np.linalg.norm(d,ord=ord)



    __dx_doc="""
        Construct Gradient of Smoothed Gaussian Pulse w.r.t. %s

        INPUT:
            domain : domain of output data ['time']/'freq'/None 
                     (return input data without change)
            mode   : type of output as complex,real-imaginary,and 
                     absolute-phase ['complex']/'re-im'/'ab-ph'
            tnorm  : L2 time normalized output if True
            tshift : Time shift flag for time domain (centered in Time if 1)


        OUTPUT:
             Estimated values of derivative w.r.t. %s
        """

    def d0(self,domain='time',mode='complex',tnorm=_df_dnorm,tshift=_df_tshift):
	    N = self.n
	    w = 2*__np.pi/N * __np.arange(_eps,N/2+1)
	    F = (-1j*w) * self.data(domain=None,tnorm=_df_tnorm)[0]
	    return __misc.outtype(F,domain,mode,tshift,tnorm=tnorm,
                                  stype=_df_dtype)

    def d1(self,domain='time',mode='complex',tnorm=_df_dnorm,tshift=_df_tshift):
	    N = self.n
	    w = 2*__np.pi/N * __np.arange(_eps,N/2+1)
	    F =-self.sigma*(w**2) * self.data(domain=None,tnorm=_df_tnorm)[0]
	    return __misc.outtype(F,domain,mode,tshift,tnorm=tnorm,
                                  stype=_df_dtype)

    def d2(self,domain='time',mode='complex',tnorm=_df_dnorm,tshift=_df_tshift):
            
	    N = self.n
	    w = 2*__np.pi/N * __np.arange(_eps,N/2+1)
	    F = - __np.log(w+0j) * self.data(domain=None,tnorm=_df_tnorm)[0]
            if self.alpha > 0:
                F[0]=0
	    return __misc.outtype(F,domain,mode,tshift,tnorm=tnorm,
                                  stype=_df_dtype)

    d0.__doc__=__dx_doc %("TAU","Tau")
    d1.__doc__=__dx_doc %("SIGMA","Sigma")
    d2.__doc__=__dx_doc %("ALPHA","Alpha")

    def d3(self,domain='time',mode='complex',tnorm=_df_dnorm,
           tshift=_df_tshift,dphi=_df_d3const):
	    """
	    Construct Gradient of Smoothed Gaussian Pulse w.r.t. PHI
	    
	    INPUT:
	        domain : domain of output data ['time']/'freq'/None 
                         (return input data without change)
	        mode   : type of output as complex,real-imaginary,and 
                         absolute-phase ['complex']/'re-im'/'ab-ph'
                tnorm  : L2 time normalized output if True
	        tshift : Time shift flag for time domain (centered in Time if 1)
                dphi   : Constant scaling number to be used for derivative 
                         to phase [see the code for more details]

	    OUTPUT:
		 Estimated values of derivative w.r.t. Phi
	    """
	    N = self.n

	    w = 2*__np.pi/N * __np.arange(_eps,N/2+1)
	    F = dphi* 1j*__np.pi* self.data(domain=None,tnorm=1)[0]

	    return __misc.outtype(F,domain,mode,tshift,tnorm=tnorm,
                                  stype=_df_dtype)



    def grad(self,domain='time',mode='complex',tnorm=_df_dnorm,
             tshift=_df_tshift,dphi=_df_d3const):
	    """
	    Construct Gradient of Smoothed Gaussian Pulse w.r.t. Parameters
	    
	    INPUT:
	        domain : domain of output data ['time']/'freq'/None 
                         (return input data without change)
	        mode   : type of output as complex,real-imaginary,and 
                         absolute-phase ['complex']/'re-im'/'ab-ph'
                tnorm  : L2 time normalized output if True
	        tshift : Time shift flag for time domain (centered in Time if 1)
                dphi   : Constant scaling number to be used for derivative 
                         to phase [see the code for more details]

	    OUTPUT:
		 Tuple of manifold's gradient with size of [N*4]
	    """
            d0=self.d0(domain,mode,tnorm,tshift)
            d1=self.d1(domain,mode,tnorm,tshift)
            d2=self.d2(domain,mode,tnorm,tshift)
            d3=self.d3(domain,mode,tnorm,tshift,dphi)
            G0=__np.array([d0[0],d1[0],d2[0],d3[0]]).transpose()
            if d0[1] != None:
                G1=__np.array([d0[1],d1[1],d2[1],d3[1]]).transpose()
            else:
                G1=None
            return (G0,G1)



    def smooth(self,s=0,fi=0,domain='time',mode='complex',tnorm=_df_tnorm,
               show=0,showparam={'domain':'time','mode':'re-im'}):
	    """
	    Uses Gaussian smoother and fractional integration to make 
            data smoother in time domain. 
            The smoothed data could be showed on a plot.

	    INPUT:
	        s      : smoothing factor [0]
		fi     : Fractional integration factor [0]
	        domain : output data domain ['time']/'freq'/None 
                         (half frequency range)
	        mode   : type of output as complex,real-imaginary,and 
                         absolute-phase ['complex']/'re-im'/'ab-ph'
                tnorm  : L2 time normalized output if True
		show   : flag to enable showing data on a plot
		showparam : parameters for plotting the data

             OUTPUT:
		Tuple of smooth data and figure if show=1.

	    """
            tmp=self.copy(s=s+self.s,fi=fi+self.fi)
	    if show==1:
		    tmp.show(**showparam)
	    return tmp.data(domain=domain,mode=mode,tnorm=tnorm)


    def plot(self,domain='time',mode='complex',label= True,tnorm=1,tamp=0,
             tshift=_df_tshift,*args,**kwargs):
	    """
	    Returns plot object of the Gaussian data in [time]/freq domain. 
            Data could be plotted in real-imagniray parts[mode='re-im'] or 
            absolute-phase parts (mode='ab-ph'). Use c='colourterm' to plot 
            in a certain colour. 

	    INPUT:
		domain : plot data domain ['time']/'freq'
		mode   : type of plot data as complex,real-imaginary,and 
                         absolute-phase ['complex']/'re-im'/'ab-ph'
		label  : Shows the labels and titles if True 
                tnorm  : L2 time normalized output if True
                tamp   : Max. of absolute amplitude in time domain
                         (Use one of these at any time and set the ther to 0, 
                         o.w. tamp will be dominant)
	        tshift : Time shift flag for time domain (centered in Time if 1)

            OUTPUT:
		plot object of data.
	    """
            _lsize=__py.rcParams['axes.labelsize']
            _tsize=__py.rcParams['axes.titlesize']

            if __py == None:
                print >>__stderr,__Disp_Err
                pass
	    if mode.lower()=='complex':
		    mode='re-im'
	    data = self.data(domain=domain,mode=mode,tshift=tshift,tnorm=tnorm)
            # Set amplitude of data if needed
            if tamp != 0:
                print >>__stderr,"Warning : tnorm argument will be dominated by tamp !"
                data0=tamp*data[0]/max(abs(data[0]))
                try:
                    data1=tamp*data[1]/max(abs(data[0]))
                except:
                    data1=None
                data=(data0,data1)

            __misc.plotdata(data,domain=domain,mode=mode,stype='real',
                            label= label,*args,**kwargs)
            if label:
                __py.title('Gaussian Manifold in '+domain.lower(),color='k',
                           style='italic',size=_tsize)    
            pass
	    


    def show(self,domain='time',mode='complex',label= True,tnorm=1,tamp=0,
             tshift=_df_tshift,*args,**kwargs):
	    """
	    Shows the generated figure by plot method of manifold object.

	    plots Gaussian data in [time]/freq domain. In frequency domain, 
            data could be plotted in real/imagniray parts[mode='re-im'] or 
            absolute/phase parts (mode='ab-ph'). Use c='colourterm' to plot 
            in a certain colour. 

	    INPUT:
		domain : plot data domain ['time']/'freq'
		mode   : type of plot data as complex,real-imaginary,and 
                         absolute-phase ['complex']/'re-im'/'ab-ph'
		label  : Shows the labels and titles if True 
                tnorm  : L2 time normalized output if True
                tamp   : Max. of absolute amplitude in time domain
	        tshift : Time shift flag for time domain (centered in Time if 1)

            OUTPUT:
		figure of plotted data.
	    """
            if __py == None:
                print >>__stderr,__Disp_Err
                pass
	    __py.figure()
	    self.plot(domain=domain,mode=mode,label= label,tnorm=tnorm,
                      tamp=tamp,tshift=tshift, *args,**kwargs)
	    __py.show()

#================================================================
#                         Define Functions
#================================================================

def smooth(I,s=0,fi=0,domain='time',mode='complex',tnorm=_df_tnorm,
           stype='real',image='time',show=0):
	"""
	Smoothing and fractional integration on input signal with Gaussian Pulse  	Note: Frequency range is similar to FFT e.g. :[0,2*pi] or [0,+f,-f]
	
	INPUT:
	     I  : Image signal in [Time(Not centeralized)]/freq Domain
	     s  : Smoothing factor
	     fi : Fractional integration order (Positive for integration)
	     domain : Output data domain ['time']/'freq'
	     mode   : type of output data as complex,real-imaginary,
                      and absolute-phase ['complex']/'re-im'/'ab-ph'
             tnorm  : L2 time normalized output if True
	     stype  : Determine wheter the data is real time signal or not 
                      ['real']/'complex'
	     image  : Domain of input data(Image)

	OUTPUT:
	     Smoothed Image (Zero-Mean)
	"""
	Nrm=__misc.norm(I)
        N=__np.size(I)
	if image.lower() == 'freq':
                print >>__stderr,"""Warning: This function needs the frequency 
                range of data to be in [0,2*pi] similar to fft output 
                without any shift !"""
		FI=I.ravel()[0:N/2+1]
	elif image.lower() == 'time':
		FI=__ft.fft(I)[0:N/2+1]        # fft=[0,+f,-f]
	else:
		raise ValueError,"""This Function can only handle data in 
                time/frequency domain."""	
        
        DC=0
	w = 2*__np.pi/N * __np.arange(_eps,N/2+1)
	fInt=1/(w)**fi
        if fi > 0:
            fInt[0]=DC
	OP=fInt*__np.exp(-(w*s)**2/ 2)
	Is=__misc.outtype(OP*FI,domain,mode,tnorm=tnorm,stype=stype)
	if show==1:
		__misc.plotdata(Is[0],domain=domain,mode=mode,stype=stype)
		__py.show()

	return Is


def movavg(I,n):
	"""
        Compute the len(n) moving average of I
	
	INPUT:
	     I  : Image signal in time
	     n  : Length of averaging window

	OUTPUT:
	     Moving Average with length N-n
	"""
        n = __np.int(n)
        N = len(I)
        assert(N>n)
        Is = __np.zeros(N-(n-1),Float)
        for i in range(n):
            Is += I[i:N-(n-1)+i]
        Is /= float(n)
        return Is



def wgauss(N,sigma,domain='freq',mode='complex'):
	"""
	Construct Gaussian wave with known sigma in frequency domain 

	INPUT:
	     N      : Length of signal (Even/Odd)
	     sigma  : scale parameter of gaussian wave
	     domain : output data domain ['time']/'freq'/None 
                      (half frequency range)
	     mode   : type of output: for frequency domain as complex,
                      real-imaginary,and absolute-phase ['complex']/'re-im'
                      /'ab-ph' and for time domain no shift 'orig', and 
                      shifted to center ['shift']


	OUTPUT:
	     Gussian Signal
	"""
	
	w= 2*__np.pi/N *__np.arange(0,N/2+1)
	F= __np.exp(- (w*sigma)**2 /2);
	return __misc.outtype(F,domain,mode)[0]

