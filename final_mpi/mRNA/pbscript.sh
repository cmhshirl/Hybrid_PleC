#PBS -l nodes=2:ppn=24
#PBS -l walltime=50:00:00
#PBS -o output.txt
#PBS -e error.txt

#PBS -M cmhshirl@vt.edu
#PBS -m abe

## nodes = total number of nodes you need
## ppn = processors per node that you will need
## walltime = amount of time your job will be allowed before being forcefully removed
#####################################################################################
cd $PBS_O_WORKDIR
## How many cores total do we have?
NO_OF_CORES=`cat $PBS_NODEFILE | egrep -v '^#'\|'^$' | wc -l | awk '{print $1}'`
##
## Main execution
echo "Job Started at: `date`"
mpiexec -np $NO_OF_CORES -machinefile $PBS_NODEFILE cell_lsodar
echo "Job Ended at: `date`"
