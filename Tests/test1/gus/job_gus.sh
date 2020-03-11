#!/bin/sh -f

#PBS -N JobName
#PBS -l     cput=6:00:00
#PBS -l walltime=6:00:00
#PBS -m ae
#PBS -l nodes=node05:ppn=1

RUNGMS=/mnt/nfs/home/common/gamess/myrungms
input="h2"
output="h2.out"

#list of possible output files
#comment what you don't need if you don't want to copy large files
#or add what you need
#hcore="hcore_energy.txt"
ehcore="hcore_energy.txt"
hcore="hcore_eigenvectors.txt"
orb="orbital.txt"
dat=$input".dat"
hforb="HForbenergy.txt"
mocoef="mocoef.txt"
dipao="dipAO.txt"
aobasis="AObasis.txt"
h1eao="h1eAO.txt"
ovao="overlapAO.txt"
xyz="xyz.txt"
h1emo="h1eMO.txt"
kinmo="kinMO.txt"
potmo="potMO.txt"
i2emo="fort.223"


##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

##############################################################
#                                                            #
#   The prologue script automatically makes a directory      #
#   on the local disks for you.  The name of this directory  #
#   depends on the job id, but you need only refer to it     #
#   using ${WORKDIR}.                                        #
#                                                            #
##############################################################

SERVER=$PBS_O_HOST
WORKDIR=/home/scratch/PBS_$PBS_JOBID
SCP=/usr/bin/scp
SSH=/usr/bin/ssh

######################################################################
#                                                                    #
#   To minimize communications traffic, it is best for your job      #
#   to work with files on the local disk of the compute node.        #
#   Hence, one needs to transfer files from your permanent home      #
#   directory tree to the directory ${WORKDIR} automatically         #
#   created by PBS on the local disk before program execution,       #
#   and to transfer any important output files from the local        #
#   disk back to the permanent home directory tree after program     #
#   execution is completed.                                          #
#                                                                    #
#   There are essentially two ways to achieve this: (1) to use the   #
#   PBS stagein and stageout utilities, or (2) to manually copy the  #
#   files by commands in this script.  The stagein and stageout      #
#   features of OpenPBS are somewhat awkward, especially since       #
#   wildcards and macros in the file lists cannot be used.  This     #
#   method also has some timing issues.  Hence, we ask you to use    #
#   the second method, and to use secure copy (scp) to do the file   #
#   transfers to avoid NSF bottlenecks.                              #
#                                                                    #
######################################################################

#####################################################
#                                                   #
#    Specify the permanent directory(ies) on the    #
#    server host.  Note that when the job begins    #
#    execution, the current working directory at    #
#    the time the qsub command was issued becomes   #
#    the current working directory of the job.      #
#                                                   #
#####################################################


PERMDIR=$PBS_O_WORKDIR
SERVPERMDIR=${PBS_O_HOST}:${PERMDIR}

echo server is $SERVER
echo workdir is $WORKDIR
echo permdir is $PERMDIR
echo servpermdir is $SERVPERMDIR
echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo ' '
echo ' '

###############################################################
#                                                             #
#    Transfer files from server to local disk.                #
#                                                             #
###############################################################

stagein()
{
 echo ' '
 echo Transferring files from server to compute node
 echo Writing files in node directory  ${WORKDIR}
 mkdir ${WORKDIR}
 cd ${WORKDIR}

 ${SCP} ${SERVPERMDIR}/* .

# echo Files in node work directory are as follows:
# ls -l
}

############################################################
#                                                          #
#    Execute the run.  Do not run in the background.       #
#                                                          #
############################################################
runprogram()
{
# program_executable < input_file > output_file
  mkdir scr
  SCR=${WORKDIR}/scr
  $RUNGMS $input nico 1 ${SCR}  > $output
  cd $SCR
   mv fort.77 soint.txt
   mv fort.223 $directory/moint.txt
   head -7 orbital1.txt >>orbital.tmp
   head -1 orbital2.txt >>orbital.tmp
   cat orbital.tmp HForbenergy.txt > orbital.txt
}

###########################################################
#                                                         #
#   Copy necessary files back to permanent directory.     #
#                                                         #
###########################################################

stageout()
{
 echo ' '
 echo Transferring files from compute nodes to server
 echo Writing files in permanent directory  ${PERMDIR}
 cd ${WORKDIR}
 SCR=${WORKDIR}/scr

 ${SCP} $output  ${SERVPERMDIR}
 ${SCP} $SCR/$ehcore  ${SERVPERMDIR}
 ${SCP} $SCR/$hcore  ${SERVPERMDIR}
 ${SCP} $SCR/$orb  ${SERVPERMDIR}
 ${SCP} $SCR/$dat  ${SERVPERMDIR}
 ${SCP} $SCR/$hforb  ${SERVPERMDIR}
 ${SCP} $SCR/$mocoef  ${SERVPERMDIR}
 ${SCP} $SCR/$dipao  ${SERVPERMDIR}
 ${SCP} $SCR/$aobasis  ${SERVPERMDIR}
 ${SCP} $SCR/$h1eao  ${SERVPERMDIR}
 ${SCP} $SCR/$ovao  ${SERVPERMDIR}
 ${SCP} $SCR/$xyz  ${SERVPERMDIR}
 ${SCP} $SCR/$h1emo  ${SERVPERMDIR}
 ${SCP} $SCR/$kinmo  ${SERVPERMDIR}
 ${SCP} $SCR/$potmo  ${SERVPERMDIR}
 ${SCP} $SCR/$i2emo  ${SERVPERMDIR}/moint.txt

# echo Final files in permanent data directory:
# ${SSH} ${SERVER} "cd ${PERMDIR}; ls -l"
 }

#####################################################################
#                                                                   #
#  The "qdel" command is used to kill a running job.  It first      #
#  sends a SIGTERM signal, then after a delay (specified by the     #
#  "kill_delay" queue attribute (set to 60 seconds), unless         #
#  overridden by the -W option of "qdel"), it sends a SIGKILL       #
#  signal which eradicates the job.  During the time between the    #
#  SIGTERM and SIGKILL signals, the "cleanup" function below is     #
#  run. You should include in this function commands to copy files  #
#  from the local disk back to your home directory.  Note: if you   #
#  need to transfer very large files which make take longer than    #
#  60 seconds, be sure to use the -W option of qdel.                #
#                                                                   #
#####################################################################

early()
{
 echo ' '
 echo ' ############ WARNING:  EARLY TERMINATION #############'
 echo ' '
 }

trap 'early; stageout' 2 9 15


##################################################
#                                                #
#   Staging in, running the job, and staging out #
#   were specified above as functions.  Now      #
#   call these functions to perform the actual   #
#   file transfers and program execution.        #
#                                                #
##################################################

stagein
runprogram
stageout 

###############################################################
#                                                             #
#   The epilogue script automatically deletes the directory   #
#   created on the local disk (including all files contained  #
#   therein.                                                  #
#                                                             #
###############################################################

exit

