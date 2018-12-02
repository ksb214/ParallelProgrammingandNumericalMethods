#!/bin/sh
#PBS -q lionxo-nuce530
#
#   Request 1 processors on 1 node
#
#PBS -l nodes=4:ppn=1
#
#   Request 4 hours of walltime
#
#PBS -l walltime=2:00:00
#
#   Request that regular output and terminal output go to the same file
#
#PBS -j oe
# change the current working directory to the directory where
# the executable file 'foo' can be found
#
cd $PBS_O_WORKDIR
#
# Writing the date and time of execution
echo " "
echo " "
echo "Job started on `hostname` at `date`"
#
#   Now we want to run the program "hello".  "hello" is in
#   the directory that this script is being submitted from,
#   $PBS_O_WORKDIR.
#
# run the executable file 'foo' using the qmpirun script
/usr/global/bin/icmpirun ./p1
#
#
#
echo " "
echo "Job Ended at `date`"
echo " "
