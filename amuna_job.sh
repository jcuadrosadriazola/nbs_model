#!/bin/bash
#PBS -l select=1:ncpus=4:mem=32gb
#PBS -l walltime=24:00:00

module load anaconda3/personal
source activate nbs

cd $PBS_O_WORKDIR

python amunas_model.py
