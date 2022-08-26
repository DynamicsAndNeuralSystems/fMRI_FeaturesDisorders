#!/usr/bin/env bash

# Define paths
export github_dir="/headnode1/abry4213/github/"
export fmri_github_dir="/headnode1/abry4213/github/fMRI_FeaturesDisorders/"

###############################################################################
# UCLA Schizophrenia
###############################################################################
export dataset_ID="UCLA_Schizophrenia"
export data_path="/headnode1/abry4213/data/UCLA_Schizophrenia/"
export subject_metadata_file="UCLA_Schizophrenia_sample_metadata.Rds"
export noise_procs="AROMA+2P AROMA+2P+GMR AROMA+2P+DiCER"

# Split data into numpy files
python3 ${fmri_github_dir}/helper_functions/data_prep_and_QC/split_MTS_into_npy.py \
--github_dir ${fmri_github_dir} \
--data_path ${data_path} \
--noise_procs $noise_procs \
--dataset_ID ${dataset_ID}

# Create sample.yaml file for this dataset
for noise_proc in $noise_procs
do
    noise_label=$(echo $noise_proc | sed "s/\+/_/g")
    Rscript $github_dir/pyspi-distribute/create_yaml_for_samples.R \
    --data_dir ${data_path}/raw_data/numpy_files/${noise_label}/ \
    --sample_metadata_file ${data_path}/${subject_metadata_file} \
    --ID_var Sample_ID \
    --label_vars Diagnosis \
    --dim_order ps \
    --overwrite
done

###############################################################################
# ABIDE ASD
###############################################################################
export dataset_ID="ABIDE_ASD"
export data_path="/headnode1/abry4213/data/ABIDE_ASD/"
export sample_metadata_file="ABIDE_ASD_sample_metadata.Rds"
export noise_procs="FC1000"

# Split data into numpy files
cmd="python3 ${fmri_github_dir}/helper_functions/data_prep_and_QC/split_MTS_into_npy.py \
--github_dir ${fmri_github_dir} \
--data_path ${data_path} \
--noise_procs $noise_procs \
--dataset_ID ${dataset_ID}"
echo $cmd
$cmd

for noise_proc in $noise_procs
do
    noise_label=$(echo $noise_proc | sed "s/\+/_/g")
    cmd="Rscript $github_dir/pyspi-distribute/create_yaml_for_samples.R \
    --data_dir ${data_path}/raw_data/numpy_files/${noise_label}/ \
    --sample_metadata_file ${data_path}/${sample_metadata_file} \
    --ID_var Sample_ID \
    --label_vars Diagnosis \
    --dim_order ps \
    --overwrite"
    echo $cmd
    $cmd
done
