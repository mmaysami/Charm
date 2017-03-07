#! /usr/bin/env python 
"""

          Show
    ======================
    It Handles plotting characterization results from data files. It Shows 
    characterization results of sfchar script (stdin). It will show initial 
    data, reconstructed data, and data overlayed with singularity orders.

    pydshow  -- Show results loaded from python data file
    rsfshow  -- Show results loaded from rsf data file  


AUTHOR:    Mohammad Maysami
           Seismic Laboratory for Imaging and Modeling (SLIM)
           Department of Earth & Ocean Sciences (EOSC)
           The University of British Columbia (UBC)

LICENSE:   You may use this code only under the conditions and terms 
           of the license contained in the file LICENSE provided with
           this source code. If you do not agree to these terms you may 
           not use this software."""

 
import os    as __os
import rsf   as __rsf
import numpy as __np
from os.path import join   as __join, abspath as __abspath
from sys     import stderr as __stderr
from numpy   import ma     as __ma
from Charm   import _Null, __Disp_Err, __path__,header
from Charm.Core.API import _df_rhdr as _tags
try:
    import pylab as __py
except ImportError:
    raise ImportError,__Disp_Err


_im_dict ={'interpolation':'nearest','aspect':'auto'}
_par_dict={'saveflag':0,'savename':'section','Null':_Null,'seisflag':0,
           'reconflag':0,'shift':0,'attrib':'alpha'}

_tsize=__py.rcParamsDefault['axes.titlesize']
_lsize=__py.rcParamsDefault['axes.labelsize']

def show_rsf(stdin,strin=0,par_dict=_par_dict.copy(),hdrlabel=True,
            im_dict=_im_dict.copy()):
    """
    (Internal use only)
    sfshow <input.rsf save=0 att=alpha Null=_Null save=0 seismic=0 recons=0
    Show characterization results of sfchar script stored in RSF form. 
    It can show original and reconstructed data, the difference of them, 
    and data overlayed with any of attributes.

    Input:
        stdin    : RSF standard input object
        strin    : Set to 1 if stdin is string, and 0 if stdin is actual stdin  
        par_dict : Dictionary of RSF par object
        hdrlabel    : Set labels from header of stdin
        im_dict  : Dictionary of plotting preferences
    """

    showflg=1
    if strin:
        stdin=__rsf.Input(stdin)
        showflg=0


    if hdrlabel:
        tags=header(stdin)
    else:
        tags=_tags

    n1 = stdin.int("n1")
    n2 = stdin.int("n2")
    n3 = stdin.int("n3")
    if n3 == None:
        n3=1        
        par_dict.update({'attrib':'0','seisflag':1})   
    else:
        print >>__stderr, "\n  ===== RSF File Loading Summary =====\n"
        print >>__stderr,"Estimation Algorithm: "+str(stdin.string("estimate"))
        print >>__stderr,"Smoothing factor used for traces: "+str(
            stdin.float("smooth"))
        print >>__stderr,"Threshold percent to pick major events: "+str(
            stdin.string("df_major"))
        print >>__stderr,"Default initial guess for singularity order: "+str(
            stdin.float("df_Alpha"))
        print >>__stderr,"Number of Initial Iteration: "+str(
            stdin.int("Niter1"))
        print >>__stderr,"Number of Final Iteration: "+str(stdin.int("Niter2"))
        print >>__stderr,"Data File Null Value: "+str(stdin.float("Null"))
        print >>__stderr,"User Defined Null Value: "+str(par_dict["Null"])

    data=__np.zeros((n3,n2,n1),dtype='f')
    stdin.read(data)
    stdin.fileclose()
    show(data,par_dict=par_dict,tags=tags,im_dict=im_dict,showflg=showflg)
    return data



def show_pyd(f,label=True,tags=_tags.copy(),im_dict=_im_dict.copy()):
    """
    (Internal use only)
    Show characterization results stored in PYD form . It can show original and 
    reconstructed data, the difference of them, 
    and data overlayed with any of attributes.

    Input:
        f       : String pointing to the pyd file
        label   : Show the Labels
        tags    : Dictionary to be used for labeling and axes limits
        im_dict : Dictionary of plotting preferences
    """
    print >>__stderr, "\n  ===== PYD File Loading Summary =====\n"
    data=__np.load(f)
    show(data,im_dict=im_dict,label=label,tags=tags)
    return data


                    
def show(data,par_dict=_par_dict.copy(),label=1,tags=_tags.copy(),
         im_dict=_im_dict.copy(),showflg=1,**kwargs):
    """
    (Internal use only)
    Main script to show images while data is loaded to python array format. 
    It can show original and reconstructed data, the difference of them, 
    and data overlayed with any of attributes.

    Input:
        data     : Loaded 3-D data in array format
        par_dict : Dictionary of RSF par object incuding 
             attrib    -- Attribute to be plotted on top of seismic data
             tags      -- Dictionary to be used for labeling and axes limits
             seisflag  -- Flag to plot original seismic data
             reconflag -- Flag to plot reconstructed data
             saveflag  -- Flag to save all the figures
             savename  -- Initial Name to be used for saved figures
             shift     -- Shift valuesof attributes with this
             Null      -- Null value in the data
        label    : Show the Labels
        tags     : Dictionary to be used for labeling and axes limits
        im_dict  : Dictionary of plotting preferences
        showflg  : Call pylab show if set (for calling from terminal)
                   for python environment set it to zero
    """
    for (name, value) in tags.items():
        exec('%s = value' %name)
        
    tmp=_par_dict.copy()
    tmp.update(par_dict)
    par_dict=tmp.copy()
    for (name, value) in par_dict.items():
        exec('%s = value' %name)

    (n3,n2,n1)=data.shape
    vars={'seismic':0,'recons':1,'amplitude':2,'scale':3,'alpha':4,'phase':5}
    seismic=data[0,:,:].transpose()
    if attrib not in ["0",""]:
        exec('%s = data[%s,:,:].transpose()' %(attrib,vars[attrib]))
    print "\n  ==> Loading data is done !\n"


    #-------------------------------
    #         Plot Figures 
    #-------------------------------
    try:
        Xextent=o2,o2+d2*n2
    except:
        Xextent=0,seismic.shape[1]
    try:
        Yextent=o1+d1*n1,o1
    except:
        Yextent=0,seismic.shape[0]

    im_dict.update({'extent':__np.concatenate([Xextent,Yextent]) })
    # Set path for saving figures
    if saveflag:
        dir=__path__[3]  # [5] for demos
        prefix=__join(dir,'Fig_res')
        fname=savename
        ext='pdf'
    
        if not __os.access(prefix,__os.F_OK):
                __os.makedirs(prefix)


    if seisflag:
        __py.figure()     # SEISMIC Data
        __py.imshow(seismic,cmap=__py.cm.gray,**im_dict);__py.colorbar()
        __py.title('Real Seismic Data',size=_tsize)
        __py.xlabel('%s (%s)' %(label2,unit2),size=_lsize)
        __py.ylabel('%s (%s)' %(label1,unit1),size=_lsize)
        if saveflag:
            __py.savefig(__join(prefix,"%s_real.%s" %(fname,ext)),dpi=100)
        #__py.grid()

    if reconflag:
        recons =data[1,:,:].transpose()
        __py.figure()     # RECONSTRUCTED Data
        __py.imshow(recons,cmap=__py.cm.gray,**im_dict);__py.colorbar()
        __py.title('Reconstructed Seismic Data (Estimated)',size=_tsize)
        __py.xlabel('%s (%s)' %(label2,unit2),size=_lsize)
        __py.ylabel('%s (%s)' %(label1,unit1),size=_lsize)
        if saveflag:
            __py.savefig(__join(prefix,"%s_reconst.%s"%(fname,ext)),dpi=100)


    if seisflag and reconflag:
        __py.figure()     # Difference of Seismic & Reconstruction
        __py.imshow(recons-seismic,cmap=__py.cm.jet,**im_dict);__py.colorbar()
        __py.title('Difference of Reconstructed and Original Seismic Data',
                   size=_tsize)
        __py.xlabel('%s (%s)' %(label2,unit2),size=_lsize)
        __py.ylabel('%s (%s)' %(label1,unit1),size=_lsize)
        if saveflag:
            __py.savefig(__join(prefix,"%s_diff.%s"%(fname,ext)),dpi=100)


    if attrib not in ["0",""]:
        limits={'alpha':(-8,0),'scale':(0,5)}
        att=eval(attrib)
        attW=att+__np.roll(att,+1,axis=0)+__np.roll(att,-1,axis=0)
        MSK_att=__ma.masked_values(att,Null)
        if attrib.lower() in ['alpha','scale']:
            limit=limits[attrib.lower()]
            MSK_att=__ma.masked_where(__np.logical_or(att<=limit[0],
                                                      att>=limit[1]),att)
            MSK_att[-1,-1],MSK_att[-1,0]=limit[0],limit[1]

        elif attrib.lower()=='phase':
            MSK_att = MSK_att%2
        MSK_att[-1,-2]=0
        __py.figure()
        __py.imshow(seismic,cmap=__py.cm.gray,**im_dict);#__py.colorbar()
        __py.imshow(MSK_att+shift,cmap=__py.cm.jet,**im_dict);__py.colorbar()
        __py.title('  Estimated %s' %attrib,size=_tsize)
        __py.xlabel('%s (%s)' %(label2,unit2),size=_lsize)
        __py.ylabel('%s (%s)' %(label1,unit1),size=_lsize)
        if saveflag:
            __py.savefig(__join(prefix,"%s_%s.%s" %(fname,attrib,ext)),dpi=100)

    if showflg:
        try:
            __py.show()
        except:
            pass
    pass



def showgray(data,label=False,tags=_tags.copy(),im_dict=None):
    """
    Show a gray scale image of 2D data
    
    Input:
        data    : 2-D data array 
        label   : Show the Labels
        tags    : Dictionary to be used for labeling and axes limits
        im_dict : Dictionary of image properties
    """
    if im_dict==None:
        im_dict=_im_dict.copy()
    if label:
        for (name, value) in tags.items():
            exec('%s = value' %name)
    try:
        Xextent=o2,o2+d2*n2
    except:
        Xextent=0,data.shape[1]
    try:
        Yextent=o1+d1*n1,o1
    except:
        Yextent=0,data.shape[0]

    im_dict.update({'extent':__np.concatenate([Xextent,Yextent]) })


    __py.imshow(data,cmap=__py.cm.gray,**im_dict)
    if label:
        __py.xlabel('%s (%s)' %(label2,unit2),size=_lsize)
        __py.ylabel('%s (%s)' %(label1,unit1),size=_lsize)
    pass



def showmasked(data,Null=_Null,limit=None,shift=0,
               im_dict=_im_dict.copy(),out=0,**kwargs):
    """
    Shows a images of 2-D data array.

    Input:
        data     : 2-D data array 
        Null     : Null value in the data
        limit    : 
        shift    : Shift valuesof attributes with this
        im_dict  : Dictionary of image properties
    """
    MSK_data=__ma.masked_values(data,Null)
    if limit !=None:
        MSK_data=__ma.masked_where(
            __np.logical_or(data<=min(limit),data>=max(limit)),data)

    #     MSK_data[-1,-1],MSK_data[-1,0]=limit[0],limit[1]
    __py.figure()
    __py.imshow(MSK_data+shift,cmap=__py.cm.jet,**im_dict);__py.colorbar()
    __py.show()
    if out:
        return MSK_data
    else:
        pass

# ************************************************************
#                Main Body (Running Script)
# ************************************************************
# if __name__ == '__main__':

#     # Set path to data file
#     _df_str="pLM01.pyd"
#     _df_file=misc.abs_file(_df_str, __path__[2])

#     f=raw_input("Name of data file?["+_df_str+"]  ")
#     if f=='':
#         if _df_file==None:
#             raise ValueError, """The default file does not exist!. 
#             Please check the path to default file as well as 'show.py' """
#         f=_df_file
#     else:
#         f=misc.search_file(f,search_path=__path__)
#         #         f=misc.abs_file(f,_df_folder)

#     print >>__stderr,"\n"
#     print >>__stderr,__abspath(f)
#     # Check if the file exists
#     if f==None:
#         raise IOError, "The input file does not exist."
#     elif not __os.path.isfile(f):
#         raise IOError, "The input file does not exist."

#     if f[-3::].lower()=='pyd':
#         data=pydshow(f)
#     else:
#         par = __rsf.Par(['a=1','b=1'])
#         stdin  = __rsf.Input(f)
#         data=rsfshow(stdin)
