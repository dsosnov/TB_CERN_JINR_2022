#!/bin/bash

#temporary
set echo on

date

module load mpi/ompi-4.1
source /home/sosnov/straw/root/bin/thisroot.sh

# echo SLURM_JOB_ID=$SLURM_JOB_ID
# echo SLURM_ARRAY_JOB_ID=$SLURM_ARRAY_JOB_ID
echo SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID
echo SLURM_ARRAY_TASK_COUNT=$SLURM_ARRAY_TASK_COUNT

echo FILE=$FILE

root -b -q -n -e 'gROOT->ProcessLine(".L link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"'${FILE}'\"))->Loop(0, '${SLURM_ARRAY_TASK_ID}', '${SLURM_ARRAY_TASK_COUNT}')")'

# echo "PWD: `pwd -P`"
# ls -l

# date

# sbatch --mem=4Gb --array=0-27 -n 1 slurm_run_apv.sh --export=FILE="run45" --job_name="run_apv"
