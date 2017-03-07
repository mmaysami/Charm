#! /usr/bin/env python 
"""
M: constrain modulus
E: Young Modulus
K: Bulk Modulus
G: Shear Modulus
K=E/3/(1-2NU)
G=E/2/(1+2NU)
M=E(1-Nu)/(1+Nu)/(1-2Nu)

NU = (3K-2G)/(6K+2G)

Vs=sqrt(G/Rho)
Vp=sqrt(M/Rho)
cof = sqrt(0.5*(1-2*NU)/(1-NU))
Vs = cof*Vp

AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""

import numpy as __np
# NU =  Poisson Ratio(0..0.5) : Liq=0.5,Rocks:0.1-0.3
_df_dic=dict(DEN1 =1400,
             DEN2 =1700,
             VP1  =1690.3,
             VP2  =1798.7,
             K1   =2.8e9,
             K2   =3.5e9,
             G1   =0.9e9,
             G2   =1.5e9,
             NU1  =0.355,
             NU2  =0.417,
             PC   =0.3116,
             DENC =1493.5,
             BETA =0.41,
             LP   =1,
             HP   =2,
             VS1  =801.784,
             VS2  =939.336)

#==============================
#   HB Model and Velocities 
#  from overall density only
#==============================
def  __HBd_vs(den,dic=_df_dic):
    """
    Generates S-wave velocity profile from density profile

    Input:
         den  : Density profile array
         dic  : Dictionary of end-member properties and binary mixture factors
                G1,G2,DEN1,DEN2,BETA,DENC are required
    Output:
         vs   : S-wave Velocity profile array
    """
    for (name,value) in dic.items():
        exec('%s = value' %name) 

    if __np.isscalar(den):
        den=list([den])
    N = len(den)
    vs=__np.zeros(N)
    for k in range(N):
        p = (den[k]-DEN1) / (DEN2-DEN1) 
        q = 1-p

        if den[k] < DENC:   # p < Pc
            Gp = 1. / ( q/G1 + p/G2 )
            vs[k] = __np.sqrt(Gp/den[k])
        elif den[k]==DEN2:
            vs[k]= __np.sqrt(G2/DEN2)
        else:               # p > Pc
            ps  = p*((den[k]-DENC)/(DEN2-DENC))**BETA
            Gp  = (1-ps)**2 / (q/G1+(p-ps)/G2) + ps*G2
            vs[k] = __np.sqrt(Gp/den[k])
    return vs



def __HBd_vp(den,dic=_df_dic):
    """
    Generates P-wave velocity profile from density profile

    Input:
         den  : Density profile array
         dic  : Dictionary of end-member properties and binary mixture factors
                K1,K2,G1,G2,DEN1,DEN2,BETA,DENC are required
    Output:
         vp   : P-wave Velocity profile array
    """

    for (name,value) in dic.items():
        exec('%s = value' %name) 

    if __np.isscalar(den):
        den=list([den])
    N = len(den)
    vp=__np.zeros(N)

    for k in range(N):
        p = (den[k]-DEN1) / (DEN2-DEN1) 
        q = 1-p
        vs = __HBd_vs(den[k], dic)
        if den[k] < DENC:
            Kp = 1. /( q/K1+p/K2 )
            vp[k] = __np.sqrt(Kp/den[k]+4*(vs**2)/3.)
        else:
            ps  = p*((den[k]-DENC)/(DEN2-DENC))**BETA
            Gp  = (1-ps)**2 / ( q/K1+(p-ps)/K2 ) + ps*K2 
            vp[k] = __np.sqrt(Gp/den[k]+4*(vs**2)/3.)
    return vp



def vHB_den(den,mode='VP',dic=_df_dic):
    """
    Generates P or S-wave velocity profile from density profile

    Input:
         den  : Density profile array
         mode : 'VP' for P-wave and 'VS' for S-wave velocity 
         dic  : Dictionary of end-member properties and binary mixture factors
                [K1,K2],G1,G2,DEN1,DEN2,BETA,DENC are required
    Output:
         V    : Velocity profile array
    """

    # mode = VP/VS
    if mode.upper() == 'VP':
        V = __HBd_vp(den,dic)
    else:
        V = __HBd_vs(den,dic)
    return V

#================================
#       HB Velocity 
# from volume fraction & Density
#================================
def __Xe(p,flag='K',dic=_df_dic):
    """
    Compute effective modulus values Ke, Ge
    Note : It is only defined and used for P > Pc 

    Input:
         p    : Volume fraction of HP element in mixture
         flag : 'K' for Bulk Modulus and 'G' for Shear Modulus  
         dic  : Dictionary of end-member properties and binary mixture factors
                [K1,K2] or [G1,G2], PC, BETA are required
    Output:
         out  : Either Ke or Ge
    """

    for (name,value) in dic.items():
        exec('%s = value' %name) 
    p=__np.array(p)
    X1=eval(flag.upper()+'1')  # K1 or G1
    X2=eval(flag.upper()+'2')  # K2 or G2
    if p==1:
        out=X2
    else:
        ps=p*((p-PC)/(1-PC))**BETA
        out = (1-ps)/((1-p)/X1+(p-ps)/X2)
    return out


def __Xp(p,flag='K',dic=_df_dic):
    """
    Compute overal modulus values Kp, Gp

    Input:
         p    : Volume fraction of HP element in mixture
         flag : 'K' for Bulk Modulus and 'G' for Shear Modulus  
         dic  : Dictionary of end-member properties and binary mixture factors
                [K1,K2] or [G1,G2], PC, BETA are required
    Output:
         out  : Either Kp or Gp
    """
    for (name,value) in dic.items():
        exec('%s = value' %name) 
    p=__np.array(p)
    X1=eval(flag.upper()+'1')  # K1 or G1
    X2=eval(flag.upper()+'2')  # K2 or G2
    if p < PC:
        out = 1./ ((1-p)/X1+p/X2)
    else:
        ps=p*((p-PC)/(1-PC))**BETA
        out = __Xe(p,flag,dic)*(1-ps)+X2*ps
    return out


def vHB_p(p,mode='VP',dic=_df_dic):
    """
    Compute velocity and density profile based on HB percolation model 
    for binary mixtures

    Input:
         p    : Volume fraction of HP element in mixture
         mode : 'VP' for P-wave and 'VS' for S-wave velocity 
         dic  : Dictionary of end-member properties and binary mixture factors
                K1,K2,G1,G2,DEN1,DEN2,BETA,DENC are required
    Output:
         v    : Velociy profile array (S-wave or P-wave velocity)
         den  : Density profile array
    """
    for (name,value) in dic.items():
        exec('%s = value' %name) 
    p = __np.array(p)
    den=(1-p)*DEN1+p*DEN2
    N=len(p)
    v=__np.zeros(N)
    for k in range(N):
        if mode.upper() == 'VS':
            v[k] = __np.sqrt(__Xp(p[k],'G',dic)/den[k])
        if mode.upper() == 'VP':
            v[k] = __np.sqrt((__Xp(p[k],'K',dic)+4*__Xp(p[k],'G',dic)/3.)/den[k])
    return v,den



#==============================
#      Velocity Models
#==============================
def vp_model(p,model='HB',dic=_df_dic):
    """
    Compute velocity profile based on different models or density profile
    for binary mixtures verus volume fraction

    NOTE:   K1, G1 > K2,G2 (1=HP,2=LP)
    
    Input:
         p    : Volume fraction of HP element in mixture
         model: model to be used for velocity 
               'DEN' ==> out : density
               'RT'  ==> out : Ray Theory Velocity
               'EMT' ==> out : Equivalent Medium Theory Velocity
               'HB'  ==> out : Herrmann-Bernabe Model Velocity

         dic  : Dictionary of end-member properties and binary mixture factors
                K1,K2,G1,G2,DEN1,DEN2,BETA,DENC are required
    Output:
         out  : Velociy/Density profile array (S-wave or P-wave velocity)
    """

    for (name,value) in dic.items():
        exec('%s = value' %name) 

    v1 = __np.sqrt((K1+4*G1/3.)/DEN1)
    v2 = __np.sqrt((K2+4*G2/3.)/DEN2)
    if __np.isscalar(p):
        p=list([p])
    N = len(p)
    den = __np.zeros(N)
    out = __np.zeros(N)

    if model.upper() == 'HB':     # HB Vel
        out,den = vHB_p(p,'VP',dic)
    else:            
        for k in range(N):
            pk = p[k]
            den[k] = (1-pk)*DEN1+pk*DEN2
            if model.upper() == 'DEN':            
                out[k] = den[k]
            elif model.upper() == 'RT':           
                out[k] = 1/((1-pk)/v1+pk/v2)
            elif model.upper() == 'EMT':          
                out[k] = __np.sqrt(1/( ((1-pk)/(DEN1*v1**2)+pk/
                                         (DEN2*v2**2)) * den[k] ))
                   
    return out


#==============================
#         HS Bounds
#==============================
def bounds(p,mode='RV',dic=_df_dic):

    """
    Compute 

    
    Input:
         p    : Volume fraction of HP element in mixture
         mode : type of bounds to be computed
                'HS' for Hashin-Shtrikman
                'RV' for Ruess & Voigt
         dic  : Dictionary of end-member properties and binary mixture factors
                K1,K2,G1,G2,DEN1,DEN2,BETA,DENC are required
    Output:
         VP_L : Lower bound for Vp
         VP_U : Upper bound for Vp
         VS_L : Lower bound for Vs
         VS_U : Upper bound for Vs
         RHO  : Density
    """
    for (name,value) in dic.items():
        exec('%s = value' %name) 

    if __np.isscalar(p):
        p=list([p])
    p=__np.array(p)
    q=1-p

    RHO = q*DEN1+p*DEN2

    if mode.upper()=='HS':
        # Hashin-shtrikman Bounds
        a1=3*K1/(3*K1+4*G1)
        b1=6*(K1+2*G1)/ (5*(3*K1+G1))
        a2=3*K2/(3*K2+4*G2)
        b2=6*(K2+2*G2)/ (5*(3*K2+G2))

        KL=K1*( 1+ p*(K2-K1) / (q*(K2-K1)*a1+K1) )
        KU=K2*( 1+ q*(K1-K2) / (p*(K1-K2)*a2+K2) )
        GL=G1*( 1+ p*(G2-G1) / (q*(G2-G1)*b1+G1) )
        GU=G2*( 1+ q*(G1-G2) / (p*(G1-G2)*b2+G2) )
    else:
        # Reuss Voigt averages
        KL=1/ (q/K1+ p/K2)
        KU=q*K1+p*K2
        GL=1/ (q/G1+ p/G2)
        GU=q*G1+p*G2

    VS_L=__np.sqrt(GL/RHO)
    VS_U=__np.sqrt(GU / RHO)
    VP_L=__np.sqrt( (KL+4*GL/3.) / RHO)
    VP_U=__np.sqrt( (KU+4*GU/3.) / RHO)


    return (VP_L,VP_U,VS_L,VS_U,RHO)

