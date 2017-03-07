"""

       API
    ==========
    Python API for RSF (Madagascar) to transfer data betwee RSF & Python

    sftrace2d
    sfsection2d    
    sfread2d
    sfread3d
    sfstore2d
    sfstore3d

AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software.
"""
#__all__=['header','sfwrite','sftrace2d','sfsection2d','sfread2d','sfread3d']

import numpy as __np
import rsf   as __rsf
import os    as __os
from   Charm.Core import Misc as __misc
from   Charm      import __path__,_df_input

_df_IOerr="The input file not found in the path: "
_df_whdr={'o1':0,'o2':0,'d1':None,'d2':None}
_df_rhdr  ={'o1':0,'o2':0,'label1':'Time','label2':'Lateral direction',
         'unit1':'','unit2':''}

def header(stdin):
    """
    (Internal use only)
    Reads header words from input RSF input (2-D or 3-D)
    
    Input:
        stdin    : RSF standard input object
    Output:
        hdr      : RSF standard headers

    """
    hdr=_df_rhdr.copy()
    params=[['stdin.int','n1','n2','n3'],
            ['stdin.float','o1','o2','o3','d1','d2','d3','Null'],
            ['stdin.string','label1','label2','unit1','unit2']]
    for ind in params:
        for prm in ind[1::]:
            stdset=ind[0]
            exec("%s = %s(prm)" %(prm,stdset))
            if eval('%s != None' %prm):
                hdr.update({prm:eval(prm)})
                
    return hdr


def sfwrite(infile,outfile,hdr_dic=None): 
        """
        Convert 2-D and 3-D python data (file) to RSF data file 
        
        Input:
             infile  : String(with ".pyd") of the python input file name (trace)
                       or a numpy array
             outfile : String(with ".rsf") of the output rsf file name
             hdr_dic : Dictionary of header info for RSF output file 
                       which includes:
                         o1,o2,[o3] : Origins
                         d1,d2,[d3] : Sampling rate
        Output:
             data    : dataset 
             hdr     : Dictionary of header information
        """
        # Load data array to python
        if isinstance(infile,__np.ndarray):
                data=__np.array(infile)
        else:
                if  not __os.path.isfile(infile):
                        infile=__misc.search_file(infile,__path__)
                        if infile==None:
                                raise IOError, _df_IOerr+_df_input
                data=__np.load(infile)

        data=__np.float32(data)
        DIM=len(data.shape)
        if DIM==2:
                n1,n2 = data.shape
                dim_dict=dict(n1=n1,n2=n2)
                hdr=_df_whdr.copy()
        else:
                n1,n2,n3 = data.shape
                dim_dict=dict(n1=n1,n2=n2,n3=n3)
                hdr=_df_whdr.copy()
                hdr.update({'o3':_df_whdr['o2'],'d3':_df_whdr['d2']})
        stdout  = __rsf.Output(outfile)

        if hdr_dic !=None:
                hdr=hdr_dic.copy()
        hdr.update(dim_dict)

        for name,value in hdr.items():
                if value != None:
                        stdout.put(name,value)

        stdout.write(data.transpose())
	stdout.fileclose()
        return data,hdr



def sfread2d(infile=_df_input,outfile=None): 
        """
        Convert 2-D RSF data file to python data (file)
        
        Input:
             infile  : String(with ".rsf") of the rsf input file name (trace)
             outfile : String(with extension) of the output data file name/ None
        Output:
             data    : RSF dataset 
             hdr     : Dictionary of header information
        """
        if not __os.path.isfile(infile):
                infile=__misc.search_file(infile,__path__)
                if infile==None:
                        raise IOError, _df_IOerr+_df_input
        par=__rsf.Par(['a=1','b=2'])
        stdin  = __rsf.Input(infile)

        hdr=header(stdin)
        n1,n2=hdr['n1'],hdr['n2']

        data=__np.zeros((n2,n1),dtype='f')
        stdin.read(data)
        stdin.fileclose()

        data=data.transpose()
        if outfile != None:
            data.dump(outfile)
        return data,hdr

def sfread3d(infile=_df_input,outfile=None): 
        """
        Convert 3-D RSF data file to python data (file)
        
        Input:
             infile  : String(with ".rsf") of the rsf input file name (trace)
             outfile : String(with extension) of the output data file name/ None
        Output:
             data    : RSF dataset 
             hdr     : Dictionary of header information
        """
        if not __os.path.isfile(infile):
                infile=__misc.search_file(infile,__path__)
                if infile==None:
                        raise IOError, _df_IOerr+_df_input
        par=__rsf.Par(['a=1','b=2'])
        stdin  = __rsf.Input(infile)

        hdr=header(stdin)
        n1,n2,n3=hdr['n1'],hdr['n2'],hdr['n3']

        data=__np.zeros((n3,n2,n1),dtype='f')
        stdin.read(data)
        stdin.fileclose()

        if outfile != None:
            data.dump(outfile)
        return data,hdr


def sfsection2d(infile=_df_input,outfile=None,f1=None,f2=None,n1=None,n2=None): 
        """
        Convert windowed section of 2-D dataset from RSF to python data (file)
        
        Input:
             infile  : String(with ".rsf") of the rsf input file name (trace)
             outfile : String(with extension) of the output data file name/ None
             f1      : Start index for windowing the data section in 1st axis
             f1      : Start index for windowing the data section in 2nd axis   
             n1      : Number of sample to be included in window in 1st axis
             n2      : Number of sample to be included in window in 2nd axis
        Output:
             data    : RSF dataset 
             hdr     : Dictionary of header information
        """
        if not __os.path.isfile(infile):
                infile=__misc.search_file(infile,__path__)
                if infile==None:
                        raise IOError, _df_IOerr+_df_input
        par=__rsf.Par(['a=1','b=2'])
        stdin  = __rsf.Input(infile)

        hdr=header(stdin)
        n1,n2,n3=hdr['n1'],hdr['n2'],hdr['n3']

        data=__np.zeros((n2,n1),dtype='f')
        stdin.read(data)
        stdin.fileclose()

        data=data.transpose()
        if f1 != None:
                data=data[f1:f1+n1,:]
        if f2 != None:
                data=data[:,f2:f2+n2]

        if outfile != None:
            data.dump(outfile)
        return data,hdr


def sftrace2d(infile=_df_input,outfile=None,tracenum=0): 
        """
        Convert windowed trace of 2-D dataset from RSF to python data(file)
        
        Input:
             infile   : String(with ".rsf") of the rsf input file name (trace)
             outfile  : String(with extension) of the output data file name/None
             tracenum : Number of trace to get from section

        Output:
             trace  : Python array for the trace
             hdr    : Dictionary of header information

        """
        if not __os.path.isfile(infile):
                infile=__misc.search_file(infile,__path__)
                if infile==None:
                        raise IOError, _df_IOerr+_df_input
        par=__rsf.Par(['a=1','b=2'])
        stdin  = __rsf.Input(infile)

        hdr=header(stdin)
        n1,n2=hdr['n1'],hdr['n2']

        section=__np.zeros((n2,n1),dtype='f')
        stdin.read(section)
        stdin.fileclose()

        section=section.transpose()
        trace = section[:,tracenum]
        trace.shape=(trace.shape[0],)

        if outfile != None:
            trace.dump(outfile)

        return trace,hdr
