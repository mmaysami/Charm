#!/usr/bin/env python
__doc__="""
DESCERIPTION:

    Executable script working on stdin & stdout RSF Files
    Find attributes of 2-D input seismic section and will generate a 3-D 
    section with original datam reconstructed data, locations, scales, 
    singularity orders, and phases. (I also saves the data to a python file 
    names 'tmp#.pyd') 
    Note that original seismic is stored in stdout without any smoothing !


SYNTAX:    sfchar <input.rsf >output.rsf solver=[FBGS]/LM/LS 
                    smooth=0  major=0.1 window=0 log=__stderr

ARGUMENTS: solver -- solver to be used in estimation of attributes
           smooth -- smoothing traces to remove high frequency content
           major  -- lower bound  of norm percentage to pick up major events 
                     in detection part
           window -- Set if analysing a part of data  ***         
           log    -- file or stderr to write logs

IO:        input.rsf  -- 2-D seismic section 
           output.rsf -- Contains seimic, reconstruction of seismic 
                         with estimated events, amplitudes of 
                         estimated events, scales, singularity orders, 
                         and phase components of events as 2-D slices 
                         respectively in third dimension. 


NOTE:      Use sfshow.py script to image the results saved to out.rsf


AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""
 
import os    as __os
import numpy as __np
import rsf   as __rsf
import Charm.Core.Misc as __misc
from sys import stderr   as __stderr
from Charm import _Null,__path__, Main as __main


# _df_dfile=__misc.abs_file(__misc._df_datafile,'data')
# par = __rsf.Par(['a=1','b=2','smooth=1.0','major=0.0','solver=BFGS'])

par    = __rsf.Par()
help     = par.int('help')

if help:
    print >>__stderr,__doc__
    raise SystemExit


stdout = __rsf.Output()
stdin  = __rsf.Input() 

solver=par.string('solver')
smooth=par.float('smooth',0)
major=par.float('major',0)
window=par.float('window',0)
# log1=par.string('log')
if solver == None :
    solver='BFGS'
for name in ['smooth','major','window']:
    if eval(name+' ==None') :
        exec(name+"=0")
log1=open('job.log','w+')
log2=log1 #__stderr

n1 = stdin.int("n1")
n2 = stdin.int("n2")
data=__np.zeros((n2,n1),dtype='f')
stdin.read(data)
data = data.transpose()


"""
PUT YOUR CODE HERE !
"""
# " Init Attribute Arrays and Estimated Data " 
att_dict={'amp':0,'sigma':1,'alpha':2,'phi':3}
for name in att_dict.keys():
    exec(name+'=_Null*__np.ones((n1,n2))')

#" Set time range for analyze"
if window:
    t_rng=__np.arange(2/0.004,4.6/0.004,dtype=__np.int32)
else:
    t_rng=__np.arange(0,data.shape[0],dtype=__np.int32)

# "Analyze each trace and store estimated attributes"
data_est=__np.zeros((n1,n2))
for rx in range(n2):
    print >>log2,"\n Trace #"+str(rx)
    input=data[t_rng,rx]
    DICT = __main.char(type='trace',input=input,user=0,solver=solver,
                        smooth=smooth,major=major,log=log1)
    
    attrib_vect=DICT['attrib_vect_est'][:,0:len(t_rng)]
    trace_est=DICT['trace_est'][0:len(t_rng)]


    for (name,value) in att_dict.items():
        exec(name+'[t_rng,rx]=attrib_vect['+str(value)+',:]')

    data_est[t_rng,rx]=trace_est

# Assign final array of data and all attributes 
OUT=_Null *__np.ones((len(att_dict)+2,n2,n1),dtype='f')
OUT[0,:,:]=data.transpose()
OUT[1,:,:]=data_est.transpose()

for (name,value) in att_dict.items():
    exec('OUT['+str(value+2)+',:,:]='+name+'.transpose()')

dir = __misc.get_path()
prefix=__path__[3]
if not __os.access(prefix,__os.F_OK):
    __os.makedirs(prefix)
OUT.dump(__os.path.join(prefix,'tmp'+str(__np.int32(smooth))+'.pyd'))
print >>log2,"\n Characterization is done."
log1.close()
"""
End of your Code .
"""

# Note that original seismic is stored in stdout without any smoothing !
n3=OUT.shape[0]
stdout.put("n1",n1)
stdout.put("n2",n2)
stdout.put("n3",n3)

stdout.put("Null",_Null)
stdout.put("smooth",smooth)
stdout.put("df_major",major)
stdout.put("df_Alpha",DICT['_df_Ai'])
stdout.put("estimate",DICT['estimate_type'])
stdout.put("Niter1",DICT['est1_dict']['niter'])
stdout.put("Niter2",DICT['est2_dict']['niter'])
stdout.write(OUT)
stdout.fileclose()
stdin.fileclose()
