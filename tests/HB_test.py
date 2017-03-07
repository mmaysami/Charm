from Charm.Extra import Percolation as pr
import numpy as __np
import pylab as py

"""
dic=dict(DEN1=1400,
     DEN2=1700,
     VP1=1690.3,
     VP2=1798.7,
     K1=2.8e9,
     K2=3.5e9,
     G1=0.9e9,
     G2=1.5e9,
     PC=0.3116,
     NU=0.22,
     DENC=1493.5,
     BETA=0.41,
     LP=1,
     HP=2)
"""
dic=pr._df_dic
for (name,value) in pr._df_dic.items():
    exec('%s = value' %name) 
        

p1=0.0
p2=1.0

den1=(1-p1)*DEN1+p1*DEN2        
den2=(1-p2)*DEN1+p2*DEN2
p=__np.arange(p1,p2,0.001)
den=__np.linspace(den1,den2,len(p))


vs=pr.vHB_den(den,'VS',dic)
VS,den0=pr.vHB_p(p,'VS',dic)

vp=pr.vHB_den(den,'VP',dic)
VP,den0=pr.vHB_p(p,'VP',dic)

(VP_LHS,VP_UHS,VS_LHS,VS_UHS,RHO)=pr.bounds(p,'RV',dic)

py.figure();py.plot(den,vs);py.title ("vs  - Density")
py.plot(den,VS_LHS);py.plot(den,VS_UHS)
py.figure();py.plot(den,vp);py.title ("vp  - Density")
py.plot(den,VP_LHS);py.plot(den,VP_UHS)

py.figure();py.plot(p,VS);py.title ("diff VS  - p")
py.plot(p,VS_LHS);py.plot(p,VS_UHS)
py.figure();py.plot(p,VP);py.title ("diff VP  - p")
py.plot(p,VP_LHS);py.plot(p,VP_UHS)
py.show()
