#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=4.5G
#$ -pe mpi 40
#$ -N pso_calculation
#$ -A KCL_Cedric 
#$ -P Gold
#$ -cwd 

module purge
module load gcc-libs/4.9.2
module load compilers/intel/2017/update1
module load mpi/intel/2017/update1/intel
module load vasp/5.4.4-18apr2017/intel-2017-update1
module load openblas/0.3.7-serial/gnu-4.9.2
module load python3/3.8


source $HOME/work_pso/kcl_pso/bin/activate
export OMP_NUM_THREADS=1

start=`date +%s`
python test_1.py  > output.out
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
