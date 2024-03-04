#!/bin/bash
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -l walltime=00:08:00

cd $PBS_O_WORKDIR

module load anaconda3/personal

source activate nbs

python qochas_model.py