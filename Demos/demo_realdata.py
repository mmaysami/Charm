#!/usr/bin/env python
import os,sys,rsf
import numpy as np
import scipy as sci
import pylab as py
import Charm as chr
tsize=py.rcParams['axes.titlesize']
lsize=py.rcParams['axes.labelsize']

realfile1=chr.search_file('veritas08w1.rsf',search_path=chr.__path__)
realfile2=chr.search_file('veritas08.rsf',search_path=chr.__path__)
dir = chr.__path__[5]
prefix=os.path.join(dir,'Intro')
if not os.access(prefix,os.F_OK):
    os.makedirs(prefix)

f=raw_input('Name of file?')
if f=='':
    f=realfile1

d,h=chr.sfread2d(f)
chr.show_rsf(rsf.Input(f))
py.title('Real Seismic Section',size=tsize)
py.savefig(os.path.join(prefix,"realsection.eps"))

py.figure()
py.plot(d[:,150],lw=1)
py.title('Trace # 850',size=tsize)
py.ylabel('Amplitude',size=lsize)
py.xlabel('Depth / Time',size=lsize)
py.savefig(os.path.join(prefix,"realtrace.eps"))



chr.show_rsf(rsf.Input(realfile2))
py.title('Real Seismic Data',size=tsize)
py.savefig(os.path.join(prefix,"realdata.eps"))
py.show()
