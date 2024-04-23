#!/bin/bash
#PBS -N HCP100_feather_to_mat
#PBS -q defaultQ
#PBS -l walltime=6:00:00
#PBS -j oe
#PBS -o /headnode1/abry4213/github/fMRI_FeaturesDisorders/cluster_output/HCP100_feather_to_mat.out
#PBS -V
#PBS -l select=1:ncpus=1:mem=20GB:mpiprocs=1

#module load Anaconda3-5.1.0
/usr/physics/Modules/3.2.8/bin/modulecmd bash load Anaconda3-5.1.0 --silent
source /usr/physics/python/anaconda3/etc/profile.d/conda.sh 

# Activate the pyspi environment
conda activate pyspi

cd /headnode1/abry4213/github/fMRI_FeaturesDisorders/code/feature_extraction/

# First, we need to convert our time-series feather file to a Matlab .mat file to be read in properly
HCP100_time_series_file_base=/headnode1/abry4213/data/HCP100/raw_data/functional_MRI/HCP100_fMRI_TS

# Run the feather_to_mat.py script with the file base as the input argument, indicating that the output file should be a mat file
python3 feather_to_mat.py ${HCP100_time_series_file_base} mat