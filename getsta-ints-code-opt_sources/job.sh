#!/bin/bash

#PBS -N HCORE
#PBS -l nodes=1:ppn=28
#PBS -j oe
#PBS -l walltime=1:00:00
#PBS -q debug


cd $PBS_O_WORKDIR

NP=`cat $PBS_NODEFILE|wc -l`

#ulimit -s unlimited


#  export OMP_NUM_THREADS=4
 export OMP_STACKSIZE="1Gb"

source /public/software/profile.d/compiler_intel-compiler-2017.5.239.sh
source /public/software/profile.d/mpi_intelmpi-2017.4.239.sh

./Get2eStaMoints he -rep 
# ./Coll_2e  ap-h2y -rep
# ./Coll_2e  ap-h2z -rep

