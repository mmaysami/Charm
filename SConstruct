# This scons script run the characterization analysis on standard input.
# It is calling 'sfchar.py <IN.rsf >OUT.rsf in scons format

# AUTHOR:   Mohammad Maysami
#           Seismic Laboratory for Imaging and Modeling (SLIM)
#           Department of Earth & Ocean Sciences (EOSC)
#           The University of British Columbia (UBC)

# LICENSE:  You may use this code only under the conditions and terms 
#           of the license contained in the file LICENSE provided with
#           this source code. If you do not agree to these terms you may 
#           not use this software.
# DATE:     Jul , 07

from rsfproj import *
from os.path import join
import Charm as chr

Data   = chr.__path__[2]    #'Data/'
Results= chr.__path__[3]    #'Results/'
Mkdir(Results)

# IN.rsf & OUT.rsf are sample names and have to be changed.
stdin=join(Data,'IN.rsf')
stdout=join(Results,'OUT.rsf')
command='./sfchar.py smooth=0 major=0.0 window=0 solver=LBFGS'
# replace IN.rsf with 2-D seismic section 
# and OUT.rsf with output file name that contains seimic, reconstruction of seismic 
# with estimated events, amplitudes of estimated events, scales, singularity orders, 
# and phase components of events as 2-D slices respectively in third dimension. 


Flow(stdout,stdin,command)
Default(Results)
End()

