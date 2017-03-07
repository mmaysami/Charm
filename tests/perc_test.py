import Charm     as __chr
import numpy     as __np
import scipy     as __sci
import numpy.fft as __ft
import pylab     as __py
from   Charm  import __path__, __Disp_Err, _eps
__misc=__chr.Misc

dt=4e-3
Nt=512
t=__np.arange(Nt)*dt
dz=0.1524
Z1,Z2,Z3=100,200,400
VP1,VP2,VP3=1000,2000,3000
DEN1,DEN2,DEN3=1000,1000,1000
RC= (VP2*DEN2-VP1*DEN1)*1./(VP2*DEN2+VP1*DEN1)
print "Reflection Coeff. : ",RC


#===================================
#   CASE 1  - 2 POINTS
#===================================
z0=[Z1,Z2]
vp0=[VP1,VP2]
den0=[DEN1,DEN2]
print "Travel Time1: ",__np.sum(__np.diff(z0)*1./__np.diff(vp0))





# Rd0,Ru0,T0=__chr.layers(z0,vp0,den0,Nt=Nt,dt=1e0*dt) 
# def layers(z,vp,rho,Nt=2001,dt=0.004):

z,vp,rho=z0,vp0,den0
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
    r    = (rho[i+1]*s1 - rho[i]*s2) / (rho[i+1]*s1+rho[i] *s2)
    r    = __np.ones((Nf,1)) * r
    # Calculate the phase shift operator
    s   = __np.ones((Nf,1))*s1
    phi  = __np.exp(1j * s * w * dz[i])
    # Calculate the R downgoing & Upgoing + T downgoing
    Rd  += (T**2)*(phi**2) *r
    Ru   = -r + (phi**2)*Ru
    T   *= phi

    # Calculate the inverse fft's and apply the taper (laplace)
    TT=__ft.irfft(__np.ravel(T).conj())
    RR=__ft.irfft(__np.ravel(Rd).conj())
    Ru=__ft.irfft(__np.ravel(Ru).conj())

#     TT1=__misc.fillfreq(T,tnorm=0,stype='real')
#     TT=    __np.real(__sci.ifft(__np.conj(TT1)))
#     RR1=__misc.fillfreq(Rd,tnorm=0,stype='real')
#     RR=    __np.real(__sci.ifft(__np.conj(RR1)))
    
    
#     nt=Nt
#     TT = __np.concatenate([T[0:nt/2,:],[__np.real(T[nt/2,:])],__np.conj(T[nt/2:1:-1,:]) ])
#     TT = __np.real(__np.fft.ifft(__np.conj(TT)))

#     RR = __np.concatenate([Rd[0:nt/2,:],[__np.real(Rd[nt/2,:])],__np.conj(Rd[nt/2:1:-1,:])])
#     RR = __np.real(__np.fft.ifft(__np.conj(RR)))




Rd0,Ru0,T0=RR,Ru,TT
__py.figure()
__py.subplot(411);__py.plot(t,Rd0);__py.ylabel('Rd')    
__py.title("COARSE  SCALE")
__py.subplot(412);__py.plot(t,T0);__py.ylabel('T')
__py.subplot(413);__py.plot(z0,vp0);__py.ylabel('Vp')
__py.subplot(414);__py.plot(z0[0:-1],__np.diff(vp0),'.');__py.ylabel("Diff Vp")







#===================================
#   CASE 2 - RANGE
#===================================
# z=__np.arange(0,512,dz)
# vp=__np.zeros(z.shape)
# den=__np.zeros(z.shape)
# ind1=__np.argmin(abs(z-300))
# vp[0:ind1],vp[ind1:]=VP1,VP2
# den[0:ind1],den[ind1:]=DEN1,DEN2
# # vp=__np.linspace(VP1,VP2,len(z))
# # den=__np.linspace(DEN1,DEN2,len(z))
# print "Travel Time2: ",__np.sum(dz/vp)

# Rd,Ru,T=__chr.layers(z,vp,den,Nt=512,dt=1e0*dt) 
# __py.figure()
# __py.subplot(411);__py.plot(t,Rd);__py.ylabel('Rd')    
# __py.title("FINE SCALE")
# __py.subplot(412);__py.plot(t,T);__py.ylabel('T')
# __py.subplot(413);__py.plot(z,vp);__py.ylabel('Vp')
# __py.subplot(414);__py.plot(z[0:-1],__np.diff(vp));__py.ylabel("Diff Vp")



