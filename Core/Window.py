"""

            Window Functions
    ===================================
    wblackman
    wboxcar
    wbutterworth
    wgaussian


AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""
# __all__=['wblackman','wboxcar','wbutterworth','wgaussian']
import numpy     as __np
import numpy.fft as __ft


def wboxcar(data,loc,sigma,lhalf=5,rhalf=5):
	"""
	
	Input:
            data  : Signal to be windowed
            loc   : Location of event
            sigma : Scale of event
            lhalf : Size of left wing of boxcar window w.r.t sigma
            rhalf : Size of right wing of boxcar window w.r.t sigma
	    
	Output: tuple of masked data and mask itself
            wdata : windowed event with  zero elsewhere
            mask  : window function
	"""
        N=len(data)
        mask=__np.zeros(N)
        mask[max(loc-lhalf*sigma,0):min(N,loc+rhalf*sigma)]=1
        wdata=data*mask
        return wdata,mask



def wgaussian(data,loc,sigma,scl=3):
	"""
	Gaussian windowing for 1-D signals
	
	Input:
            data  : Signal to be windowed
            loc   : Location of event
            sigma : Scale of event
            scl   : Scale factor w.r.t. sigma (Sigma of Gaussian=sigma*scl) 
	    
	Output: tuple of masked data and mask itself
            wdata : windowed event with  zero elsewhere
            mask  : window function
	"""

        N=len(data)        
        t=__np.arange(0,N,dtype=int)
        scl=float(scl*sigma)
        mask= __np.exp(-0.5* ((t-loc)/scl)**2 ) 
        wdata=data*mask
        return wdata,mask


def wblackman(data,loc,sigma,left=3,right=3,taper=3):
	"""
	Blackman windowing for 1-D  signals
	
	Input:
            data  : signal to be windowed / length of signal as scalar
            loc   : Location of event
            sigma : scale of event
            left  : size of left wing of flat part w.r.t sigma
            right : size of right wing of flat part  w.r.t sigma
            taper : size of left & right taper w.r.t sigma
	    
	Output: tuple of masked data and mask itself
            wdata : windowed event with  zero elsewhere
            mask  : window function
            """
        if __np.isscalar(data):
                N=data
        else:
                N = len(data)
        rwbound=round(right*sigma)
        lwbound=round(left*sigma)
        eps=2*sigma
        twidth=eps*taper
        n=__np.arange(0,twidth)
        a0=0.42;a1=0.5;a2=0.08
        b0=a0-a1*__np.cos(2*__np.pi*n/(2*twidth))+a2*__np.cos(4*__np.pi*n/(2*twidth))
        bn=b0[::-1]

        mask=__np.zeros(__np.shape(data))
       
        mask[__np.maximum(0,loc-lwbound):__np.minimum(loc+rwbound,N)]=1
        # Check Left part of window 
        if __np.maximum(0,loc-lwbound-len(b0)) != 0:
                mask[__np.maximum(0,loc-lwbound-len(b0)):loc-lwbound]=b0
        else:
                mask[0:loc-lwbound]=1
        # Check Right part of window 
        if __np.minimum(loc+rwbound+len(bn),N) != N:
                mask[loc+rwbound:__np.minimum(loc+rwbound+len(bn),N)]=bn
        else:
                mask[loc+rwbound:N]=1
        if __np.isscalar(data):
                wdata=None
        else:
                wdata=data*mask
        return (wdata,mask)



def wbutterworth(data,loc,sigma,width=60,ordl=3,ordr=3):
	"""
	ButterWorth Wndowing for 1-D signals
	
	Input:
            data  : Signal to be windowed / length of signal as scalar
            loc   : Location of event
            sigma : Scale of event
            width : Factor for Width of window(Default width ~ 0.35*width*sigma)
            ordl  : Butterworth order of left wing 
            ordr  : Butterworth order of right wing 
 
	    
	Output: tuple of masked data and mask itself
            wdata : Windowed event with  zero elsewhere
            mask  : Window function
            """
        """
        # lenl  : Nonzero length of left half 
        # lenr  : Nonzero length of right half

        # left  : size of left wing of window w.r.t sigma
        # right : size of right wing of window  w.r.t sigma
        """
        left=0.1
        right=0.1
        
        if __np.isscalar(data):
                N = data
        else:
                N = len(data)
        width*=sigma
        t=__np.arange(0,N,dtype=int)
        loc=int(loc)
        # Left & Right half of the window function
        #ordl*=sigma
        #ordr*=sigma
        locl=left*width
        locr=right*width
        b0 = (1 +(t/locl)**(2*ordl)) ** -1 
        bn = (1 +(t/locr)**(2*ordr)) ** -1


        mask=__np.zeros(N)
        # Check Left & Right part of window        
        mask[loc:0:-1]=b0[0:loc]
        mask[loc:N]=bn[0:N-loc]
        lenr=N-loc
        lenl=loc

        if __np.isscalar(data):
                wdata=None
        else:
                wdata=data*mask
        return wdata,mask#,lenl,lenr

