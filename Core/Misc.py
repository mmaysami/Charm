"""
         Miscellaneous functions for Charm Package
    ==================================================

    get_path    --  Returns absolute path of the executing script
    abs_file    --  Returns absolute path of a file 
    search_file --  Given a search path, find file

    vector      --  Turn any 1-D vector to a column/row vector
    min2d       --  Find min. of a matrix and return its value and index
    max2d       --  Find Max. of a matrix and return its value and index  
    mse         --  Gives Mean Square Error with respect to nput data vectors 
    norm        --  Gives different order of norm of data  
     
    hilbert     --  Hilbert transform of a real data vector
    phase_shift --  Phase shift to data vector
    fillfreq    --  Convert pos. frequency spectrum to full frequency spectrum 
    outtype     --  Manages returning data in different domain and modes
    plotdata    --  Returns plot object of data in different domains and modes.



AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""
# __all__=['get_path','abs_file','search_file','vector','min2d','max2d',
#          'mse','norm','hilbert','phase_shift','fillfreq','outtype','plotdata']

from sys import stderr as __stderr
import numpy as __np
import numpy.fft as __ft

try:
        import pylab as __py
except ImportError:
        from Charm import Disp_Err as __Disp_Err
        __py=None
 
from Charm import __path__, _eps

#=============================
#   File Handling
#=============================
from sys import argv as __argv
import os.path as __path


def get_path():
	"""
        Returns absolute path of the executing script
	"""
        script=__argv[0]
        pathname = __path.dirname(script)
        dir=__path.abspath(pathname)
        return dir



def abs_file(file,folder):
	"""
        Returns absolute path to 'file' sitting inside 'folder' 
        relative to Charm folder     (inside executing script)
	"""
        if file in [None,'']:
                return None

        if not folder in __path__:
                dir =  __path__[0]
                pi   = __path.join(dir,folder)
        else:
                pi=folder

        if __path.exists(__path.join(pi, file)):
                return __path.abspath(__path.join(pi, file))                
        else:
                print "[Misc.py in abs_file] Warning: File not found!: "+file
                return None


def search_file(file, search_path=__path__):
	"""
        Extensive search of a file in a search path
	
	Input:
	    file        : Name of the file to be searched for
	    search_path : List of Locations to look for the file

	Output:
	    Absolute path to the file
	"""
        if isinstance(search_path,str):
                search_path=[search_path]
        file=__path.basename(file)
        if file in [None,'']:
                return None
        file_found=0
        for pi in search_path:
                if __path.exists(__path.join(pi, file)):
                        file_found = 1
                        break
        if file_found==1:
                return __path.abspath(__path.join(pi, file))
        else:
                print "Warning: File not found in the search path: "+file 
                return None
               

#=============================
#   Simple Vector Operators
#=============================
def vector(a,type='col'):
	"""
	Turn 1-D [un]structured python arrays to column/row vectors
	
	Input:
	    a    : 1-D array 
	    type : determine [column] or row vector is required

	Output:
	    Output structured vector
	"""
        v=__np.array(a).ravel()
        v.shape=(len(v),1)
        if type.lower() == 'row':
                v=v.transpose()     
        return v
        

def min2d(A):
	"""
	Gives Minimum value of a 2-D array
	
	Input:
	    A : 2-D array 
	    
	Output:
	    MIN  : Minimum value 
            Ind  : Index of min value as a tuple
            flag : False if it is the only minimum value 
	"""
        ind0=__np.argmin(A,axis=1)
        M=__np.zeros(len(ind0))
        for i in __np.arange(len(ind0)):
                M[i]=A[i,ind0[i]]

        ind1=__np.argmin(M)       
        MIN=M[ind1]
        Ind=(ind1,ind0[ind1])
        chk=(A==MIN)
        chk[Ind]=False
        flag=__np.any(chk)
        return (MIN,Ind,flag)

def max2d(A):
	"""
	Gives Maximum value of a 2-D array
	
	Input:
	    A : 2-D array 
	    
	Output:
	    MAX  : Maximum value 
            Ind  : Index of max value as a tuple
            flag : False if it is the only maximum value 
	"""
        ind0=__np.argmax(A,axis=1)
        M=__np.zeros(len(ind0))
        for i in __np.arange(len(ind0)):
                M[i]=A[i,ind0[i]]

        ind1=__np.argmax(M)       
        MAX=M[ind1]
        Ind=(ind1,ind0[ind1])
        chk=(A==MAX)
        chk[Ind]=False
        flag=__np.any(chk)     
        return (MAX,Ind,flag)

def mse(data1,data2=None):
	"""
	Gives Mean Square Error with respect to input data vector(s) 
	
	Input:
	    data1 : First data vector 
	    data2 : None/ Second data vector (same dimension as data1)  
	    
	Output:
	    Mean Square Error value
	"""
        data=__np.array(data1).ravel()
        if data2 !=None:
                data -= __np.array(data2).ravel()
        return sum( data**2 )/len(data)
        
def norm(data,mode='complex',ord=2):
	"""
	Gives different order of norm of data
	
	Input:
	    data : data that its norm is requried 
	    mode : type of input as complex,real-imaginary,and absolute-phase ['complex']/'re-im'/'ab-ph'
	    ord  : Norm order  
	    
	Output:
	    Norm of the data
	""" 
        if mode.lower()=='re-im':
                data=data[0]+1j*data[1]
        elif mode.lower()=='ab-ph':        
                data=data[0]
        elif mode.lower()=='complex':
                if isinstance(data,tuple):
                        data=data[0]
        else:
                raise ValueError,'Unsupported value for mode. Use complex, re-im, or ab-ph !'
        return __np.linalg.norm(data,ord=ord)

#=============================
#  Data Vector Manilpulations
#   Domains and plots
#=============================
def hilbert(d):
	"""
	Hilbert transform of a real data vector

	INPUT:
	    d  : real data in time/frequency domain

	OUTPUT:
	    Hilbert transform(real) of real input
        """
        Fd=__np.array(d)
        if __np.isreal(Fd).all():
                Fd=__ft.fft(Fd)
        w = __np.arange(len(d)/2+1)
        w[-1]=0
        w = __np.concatenate( ([w,-w[-2:0:-1]]))
        Hd = -1j * __np.sign(w)*Fd
        return __ft.ifft(Hd).real


def phase_shift(d,phi=0,tnorm=0,stype='real'):
        """
	Phase shift to data vector

	INPUT:
            d      : Data vector
            phi    : Phase shift
            tnorm  : L2 normalized output if True
            stype  : type of signal in time domain ['real']/'complex'

	OUTPUT:
             Phase shifted data vector   
        """
        N=len(d.ravel())                
        Fd=__ft.fft(d)[0:N/2+1]
        # Apply phase shift to ony postive frequency components
        Fd[1:N/2]*=__np.exp(1j*__np.pi*phi)
        dp=outtype(Fd,domain='time',tnorm=tnorm,stype=stype)[0]
        return dp



def fillfreq(FH,tnorm=0,stype='real',axis=0):
	"""
	Transform data with positive frequency to full frequency range for ifft 
	
	Input:
	    FH    : fourier domain data for zero and positive frequency range
            tnorm : L2 time normalized output if True, unchanged by default
	    stype : type of signal in time domain ['real']/'complex' 
	    axis  : axis of concatenation
	    
	Output:
	    Array of fourier domain data for whole (0,+,-) frequency range
	"""

	FH=__np.array(FH)
	N=2*( FH.shape[axis] -1)
	shap=__np.array(FH.shape)
	shap[axis]+=shap[axis]-2

    
	if axis==0:
		Nyqfreq=FH[N/2]
	else:
		Nyqfreq=FH[:,N/2]

	if stype.lower() == 'complex':
		negf=0
	elif stype.lower() == 'real':
		negf=1
		Nyqfreq=__np.real(Nyqfreq)
	else:
		raise ValueError,"Use 'real','complex' for signal type."			

	if axis==0:
	    if isinstance(FH,__np.matrix):
		Nyq=Nyqfreq
	    else:
		Nyq=[Nyqfreq]

	    F=__np.concatenate(( FH[0:N/2],Nyq ,negf*__np.conj(FH[N/2-1:0:-1] ) ),axis=axis)
	else:
	    F=__np.concatenate(( FH[:,0:N/2], Nyqfreq.reshape(shap[0],1), negf*__np.conj(FH[:,N/2-1:0:-1] )),axis=axis)

        if tnorm:
                nrm=norm(F,mode='complex',ord=2)
                if nrm  > _eps:
                        F/=norm(F,mode='complex',ord=2)/__np.sqrt(N)
                else :
                        pass
                        #print >>__stderr, "Warning: Norm of data is approximately zero"
	return F.reshape(shap)



def outtype(FH,domain='time',mode='complex',tshift=0,tnorm=0,stype='real',axis=0):
	"""
	Manages returning data in different domain and modes(for internal use).
	Note: Frequency range is similar to FFT e.g. :[0,2*pi] or [0,+f,-f]

	INPUT:
	    FH     : complex data in frequency domain for zero and positive frequencies only (f=0,+) 
	    domain : domain of output data ['time']/'freq'/None (return input data without change)
	    mode   : type of output as complex,real-imaginary,and absolute-phase ['complex']/'re-im'/'ab-ph'
	    tshift : Time shift flag for time domain (centered in time if 1)
            tnorm  : L2 normalized output if True
	    stype  : type of signal in time domain ['real']/'complex' 
	    axis   : axis to fill for negative frequencies

	OUTPUT:
  	     Tuple of data parts(full range).If data should have only oe part then the second part is None
	"""
        FH=__np.array(FH)
	if  domain == None :
		return (FH,None)

	F=fillfreq(FH,tnorm=tnorm,stype=stype,axis=axis)
	if mode==None:
		mode='complex'
	if domain.lower() == 'time' :
		# Recast to real if imaginary is close to zero
		data=__ft.ifft(F)
		if stype.lower() == 'real':
			data=__np.real(data)
		if tshift:
			data=__ft.ifftshift(data)
	elif domain.lower() == 'freq' :
		data=F
	else:
		raise ValueError,"Use 'time','freq', or None for domain."	


	if mode.lower()=='complex':
		return (data,None)
	elif mode.lower()=='re-im':
		return (__np.real(data),__np.imag(data))
	elif mode.lower()=='ab-ph':
		Rdata=__np.real(data)
		Idata=__np.imag(data)
                # 		PHdata=__np.arctan2(Idata,Rdata)
                PHdata=__np.angle(data)
		return (abs(data),PHdata)
	else:
		raise ValueError,"Use 'complex','re-im', or 'ab-ph' for mode."

		
def plotdata(data,domain='time',mode='complex',stype='real',label= True, tshift=False,*args,**kwargs):
    """
    Returns plot object  data in [time]/freq domain. Data could be plotted in 
    real-imagniray parts (mode='re-im') or absolute-phase parts (mode='ab-ph')

    INPUT:
    data   : tuple of data consisting of two elements 
             (last is None if  not required)
	domain : plot data domain ['time']/'freq'
	mode   : type of plot data as complex,real-imaginary,and 
	         absolute-phase ['complex']/'re-im'/'ab-ph'
	stype  : Determine wheter the data is real time signal or 
	         not ['real']/'complex'
	label  : Shows the labels and titles if True
	tshift : Shift origin of time data to the center if True

    OUTPUT:
	Generates a figre of data
    """
    if __py == None:
        print >>__stderr, __Disp_Err
        pass
    
    _tsize=__py.rcParams['axes.titlesize']
    _lsize=__py.rcParams['axes.labelsize']

    if not(isinstance(data,tuple)):
            data=(data,None)
    
    if not kwargs.has_key('lw'):
            kwargs['lw']=1
    if mode==None:
	    mode='complex'

    data0=data[0]
    data1=data[1]

    if mode.lower()=='complex' and stype.lower()=='complex':
	    print >>__stderr, "Warning: Complex data can not be plotted. Data recasted to real-imaginary parts !"
	    data0=__np.real(data[0])
	    data1=__np.imag(data[0])
	    mode='re-im'

    if tshift:
	    data0=__ft.ifftshift(data0)
	    if data1 != None:
		    data1=__ft.ifftshift(data1)

    #__py.clf()
    if domain.lower()=='time' and stype.lower()=='real':
	    __py.plot(__np.arange(len(data0.ravel())),data0,*args,**kwargs)
	    if label != False:
		    __py.xlabel('Time samples',color='k',
                                style='italic',size=_lsize)
		    __py.ylabel('Amplitude',color='k',
                                style='italic',size=_lsize)
		    __py.title('Data in time',color='k',
                               style='italic',size=_tsize)
	    return

    if __np.logical_and(domain.lower()=='freq', label):
	    x=__np.linspace(0,2*__np.pi,len(data0.ravel()))
    else:
	    x=__np.arange(len(data0.ravel()))

    if mode.lower() == 're-im':
	    __py.subplot(212)
	    __py.plot(x,data1,'r',*args,**kwargs)
	    if label != False :
		    __py.xlabel(domain.title()+' samples',size=_lsize)
		    __py.ylabel('Imaginary',size=_lsize)

	    __py.subplot(211)
	    __py.plot(x,data0,*args,**kwargs)
	    if label != False :
		    __py.ylabel('Real',size=_lsize)
		    __py.title('Data in '+domain.lower(),size=_tsize)

    elif mode.lower() == 'ab-ph':	    
	    __py.subplot(212)
	    __py.plot(x,data1,'r',*args,**kwargs)
	    if label != False :
		    __py.xlabel(domain.title()+' samples',size=_lsize)
		    __py.ylabel('Phase',size=_lsize)

	    __py.subplot(211)
	    __py.plot(x,data0,*args,**kwargs)
	    if label != False :
		    __py.ylabel('Absolute value',size=_lsize)
		    __py.title('Data in '+domain.lower(),size=_tsize)
