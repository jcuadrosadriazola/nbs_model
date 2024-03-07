#!/bin/bash
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -l walltime=24:00:00

module load anaconda3/personal
source activate nbs

cd $PBS_O_WORKDIR

python qochas_model.py
