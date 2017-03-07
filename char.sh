#! /bin/bash
#PBS -q sngl
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -j oe
#PBS -N Job_Char100
#PBS -M mmaysami@eos.ubc.ca


# AUTHOR:   Mohammad Maysami
#           Seismic Laboratory for Imaging and Modeling (SLIM)
#           Department of Earth & Ocean Sciences (EOSC)
#           The University of British Columbia (UBC)

# LICENSE:  You may use this code only under the conditions and terms 
#           of the license contained in the file LICENSE provided with
#           this source code. If you do not agree to these terms you may 
#           not use this software.
# DATE:    Jul , 07

# SYNTAX:  type "qsub char.sh" on command line to submit the job to cluster

EXEC="scons"
ARGS=""
WORKDIR="$HOME/PyTools/Charm"

#### jump to working directory
if test -d "$WORKDIR"; then
	cd $WORKDIR || exit 1
else
	cd $PBS_O_WORKDIR || exit 1
fi

echo ========== Executing ==========
$EXEC $ARGS


