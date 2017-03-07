#!/bin/bash

# AUTHOR:   Mohammad Maysami
#           Seismic Laboratory for Imaging and Modeling (SLIM)
#           Department of Earth & Ocean Sciences (EOSC)
#           The University of British Columbia (UBC)

# LICENSE:  You may use this code only under the conditions and terms 
#           of the license contained in the file LICENSE provided with
#           this source code. If you do not agree to these terms you may 
#           not use this software.
# DATE:    Jul , 07

# SYNTAX:  type "./full_char.sh in command line to submit the job to cluster


# Fetch and load input data into files
cd `dirname $0`
cd ./data
echo
echo Fetch data from FTP server 
pwd
scons

# Call Job from Cluster
cd ..
echo
echo Call job from cluster
pwd
qsub char.sh
