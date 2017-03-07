#!/usr/bin/env python
import numpy as __np
import Charm as __chr

# infile=__np.array([[88,1.5,6,7],[10,12,6.3,3.1],[2.0,8.7,33,12.4]],'f')
a=__np.random.randn(6,47)
infile='ddd.pyd'
a.dump(infile)
outfile='ddd.rsf'
hdr_dic=None


__chr.sfwrite(infile,outfile,hdr_dic)
