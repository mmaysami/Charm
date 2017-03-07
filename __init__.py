"""

       ChaRM Package for characterizing reflectivity models
    ==========================================================
    
    Misc        --  Miscellaneous Function for data and file handling
    API         --  RSF Programming interface (Read & Write RSF files)
    Synthetize  --  Generate Synthetic traces with fractional splines and spikes
    Manifold    --  Gaussian Image Appearance Manifolds module
    Cwt         --  Continious wavelet transform module
    Window      --  Window functions to be used for segmentation of events
    Steps       --  Steps of Characterization with Detection-Estimation method
    Main        --  Actuall Characterization submodules to work with real or    
                    synthetic data


AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""


from os.path import join    as __join,abspath as __abspath
from os      import environ as __env,makedirs as __makedirs
from os      import access  as __access,F_OK   as __F_OK
from sys     import version_info as __vinfo

#=========================================
#    Set Package global variables
#=========================================
__author__ = """
   Mohammad Maysami  
   Seismic Laboratory for Imaging and Modeling (SLIM)
   Department of Earth & Ocean Sciences (EOSC)
   The University of British Columbia (UBC)
"""
__license__ =  """ 
LICENSE INFORMATION

You may use this code only under the conditions and terms of the
license contained in the file LICENSE provided with this source
code. If you do not agree to these terms you may not use this
software."""
__Disp_Err = """ Pylab(Matplotlib) Module is neither installed nor set properly on your machine."""

_df_input  = "INPUT_XXXX.rsf"  # Default input rsf file to read seismic data
_df_rx     = 790               # Trace number in data file to be analyzied
_df_xlabel = 'Time (Sample)'   # Default xLabel for figure Time/Depth/Sampl
_Null      =-100               # Null value for vectors
_eps       = 1e-16             # pre-defined epsilon


#    Python Version Check
assert __vinfo > (2,4) , ( "Please use a version"
                           " of Python greater than 2.4" )
#=========================================
#     Set different path arguments
#=========================================
__path__[0] = __abspath(__path__[0])
__path__.append(__join(__path__[0],'Core'))              # __path__[1]
for key in ['Data','Results','pydata','Demos','Extra']:  # __path__[i]
    if __env.has_key('CHARM_'+key):
        __path__.append(__env['CHARM_'+key])  
    else:
        __path__.append(__join(__path__[0],key)) 
    if not __access(__path__[-1],__F_OK):
        __makedirs(__path__[-1])

#=========================================
#            Import modules 
#=========================================

from Core.Misc         import *
from Core.API          import *
from Core.Manifold     import *
from Core.Cwt          import *
from Core.Synthesize   import *
from Core.Window       import *
from Core.Steps        import *
from Core.Main         import *
from Core.Show         import *

from Extra.Logs        import *
from Extra.Percolation import *
from Extra.Model       import *
#=========================================
#     Delete un neccessary variables 
#           from namespace
#=========================================
try:
    del __ft,__np,__py,__os,__sys,__cpk,__copy,__path,__chr,__os,__sci
except:
    pass



