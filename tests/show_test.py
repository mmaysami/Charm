#! /usr/bin/env python 

# Author: M. Maysami
#         Seismic Laboratory for Imaging and Modeling
#         Department of Earch & Ocean Sciences
#         The University of British Columbia
#         
# Date  : April, 07
 
import os
import numpy as np
import Charm as chr
import rsf

try:
    import pylab as py
    MA  = py.cm.colors.ma
    MA2 = py.matplotlib.numerix.ma
except ImportError:
    raise ImportError,chr.Pylab_Error
    

"""
Show characterization results of sfchar script. It will show initial data, reconstructed data, and data overlayed with singularity orders.
"""
saveflg=0

# Set path to data file
df_file=chr.search_file('LM01.pyd')
prefix=os.path.join( chr.__path__[4],'Seismic' )
if not os.access(prefix,os.F_OK):
    os.makedirs(prefix)



f=raw_input("Name of data file?['LM01.pyd']")
if f=='':
    f=df_file
else:
    f=chr.search_file(f)

# Lod data & attributes
data=np.load(f)
n3,n2,n1 = data.shape

var=['seismic','recons','amp','sigma','alpha','phi']
for ind in range(len(var)):
    exec(var[ind]+'=data['+str(ind)+',:,:].transpose()')


#================================
#         Plot Figures 
#================================

# ORIGINAL Data
# py.matshow(seismic,cmap=py.cm.gray,fignum=1);py.colorbar()
# py.title('Real Seismic Data (Migrated)',size=18)
# py.xlabel('Offset',size=18)
# py.ylabel('Depth / Time',size=18)
# if saveflg:
#     py.savefig(os.path.join(prefix,"section_real.eps"))


# RECONSTRUCTED Data
# py.matshow(recons,cmap=py.cm.gray,fignum=2);py.colorbar()
# py.title('Reconstructed Seismic Data (Estimated)',size=18)
# py.xlabel('Offset',size=18)
# py.ylabel('Depth / Time',size=18)
# if saveflg:
#     py.savefig(os.path.join(prefix,"section_reconst.eps"))



# Estiamated ATTRIBUTES
alpha2=alpha+np.roll(alpha,+1,axis=0)+np.roll(alpha,-1,axis=0)
alpha2[-1,-1],alpha2[-1,-2]=0.5,-5
alpha2=MA.array(alpha2)
MSK_alpha=MA.masked_where(np.logical_or((alpha2==0),(alpha2 <-5)),alpha2)
# MSK_alpha=MA.masked_values(alpha2,0)


py.figure(3)
py.imshow(seismic,cmap=py.cm.gray,interpolation='nearest');#py.colorbar()
py.imshow(MSK_alpha,cmap=py.cm.jet,interpolation='nearest');py.colorbar()
# py.matshow(seismic,cmap=py.cm.gray,alpha=1,fignum=3);#py.colorbar()
# py.matshow(MSK_alpha,cmap=py.cm.spectral,alpha=1,fignum=3);py.colorbar()
py.title('  Estimated Singularity Orders',size=18)
py.xlabel('Offset',size=18)
py.ylabel('Depth / Time',size=18)
py.show()
if saveflg:
    py.savefig(os.path.join(prefix,"section_alpS10X.eps"))


















#===============================================================

# py.figure()
# py.matshow(seismic[520:580,:],cmap=py.cm.gray,alpha=1,fignum=4);py.colorbar()
# py.matshow(MSK_alpha[520:580,:],cmap=py.cm.spectral,alpha=1,fignum=4);py.colorbar()
# py.title('  Estimated Singularity Orders',size=18)
# py.xlabel('Offset',size=18)
# py.ylabel('Depth / Time',size=18)
# if saveflg:
#     py.savefig(prefix+"section_alpZ1.eps")

# py.figure()
# py.matshow(seismic[675:745,:],cmap=py.cm.gray,alpha=1,fignum=5);py.colorbar()
# py.matshow(MSK_alpha[675:745,:],cmap=py.cm.spectral,alpha=1,fignum=5);py.colorbar()
# py.title('  Estimated Singularity Orders',size=18)
# py.xlabel('Offset',size=18)
# py.ylabel('Depth / Time',size=18)
# if saveflg:
#     py.savefig(prefix+"section_alpZ2.eps")


# py.figure()
# py.matshow(seismic[1060:1140,:],cmap=py.cm.gray,alpha=1,fignum=6);py.colorbar()
# py.matshow(MSK_alpha[1060:1140,:],cmap=py.cm.spectral,alpha=1,fignum=6);py.colorbar()
# py.title('  Estimated Singularity Orders',size=18)
# py.xlabel('Offset',size=18)
# py.ylabel('Depth / Time',size=18)
# if saveflg:
#     py.savefig(prefix+"section_alpZ3.eps")
# py.show()











#===============================================================

# py.figure()
# py.imshow(seismic[520:580,:],cmap=py.cm.gray,alpha=1);#py.colorbar()
# py.imshow(MSK_alpha[520:580,:],cmap=py.cm.spectral,alpha=1);#py.colorbar()
# py.title('   Singularity Orders (Estimated)',size=18)
# py.axis('off')
# if saveflg:
#     py.savefig(prefix+"section_alpZ1.eps")

# py.figure()
# py.imshow(seismic[675:745,:],cmap=py.cm.gray,alpha=1);#py.colorbar()
# py.imshow(MSK_alpha[675:745,:],cmap=py.cm.spectral,alpha=1);#py.colorbar()
# py.axis('off')
# if saveflg:
#     py.savefig(prefix+"section_alpZ2.eps")


# py.figure()
# py.imshow(seismic[1060:1140,:],cmap=py.cm.gray,alpha=1);#py.colorbar()
# py.imshow(MSK_alpha[1060:1140,:],cmap=py.cm.spectral,alpha=1);#py.colorbar()
# py.xlabel('Offset',size=18)
# py.axis('off')
# if saveflg:
#     py.savefig(prefix+"section_alpZ3.eps")
# py.show()

