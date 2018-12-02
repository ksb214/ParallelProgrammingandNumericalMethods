# This is a sample PBS script. It will request 1 processor on 1 node
# for 4 hours.
#
#PBS -q unisys-nuce530
#
#   Request 1 processors on 1 node
#
#PBS -l nodes=2:ppn=2
#
#   Request 4 hours of walltime
#
#PBS -l walltime=8:00:00
#
#   Request that regular output and terminal output go to the same file
#
#PBS -j oe
#
#   The following is the body of the script. By default, 
#   PBS scripts execute in your home directory, not the 
#   directory from which they were submitted. The following
#   line places you in the directory from which the job
#   was submitted.
# 
cd $PBS_O_WORKDIR
#
#   Now we want to run the program "hello".  "hello" is in
#   the directory that this script is being submitted from,
#   $PBS_O_WORKDIR.
#
echo " "
echo " "
echo "Job started on `hostname` at `date`"
time ./a.out
echo " "
echo "Job Ended at `date`"
echo " "

