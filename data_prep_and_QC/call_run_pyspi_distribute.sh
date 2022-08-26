#!/usr/bin/env bash
export github_dir=/headnode1/abry4213/github/
export config_file=${github_dir}/fMRI_FeaturesDisorders/data_prep_and_QC/pyspi14_config.yaml
export pyspi_script_dir=${github_dir}/pyspi-distribute/
export email="abry4213@uni.sydney.edu.au"

# UCLA Schizophrenia
dataset_ID="UCLA_Schizophrenia"
export data_path=/headnode1/abry4213/data/${dataset_ID}
for noise_label in AROMA_2P AROMA_2P_GMR AROMA_2P_DiCER
do
    python $pyspi_script_dir/distribute_jobs.py \
    --data_dir ${data_path}/raw_data/numpy_files/${noise_label}/ \
    --compute_file $pyspi_script_dir/pyspi_compute.py \
    --template_pbs_file $pyspi_script_dir/template.pbs \
    --pyspi_config $config_file \
    --sample_yaml ${data_path}/raw_data/numpy_files/${noise_label}/sample.yaml \
    --pbs_notify a \
    --email $email --walltime_hrs 8 --cpu 1 --mem 20 --table_only
done


# ABIDE ASD
dataset_ID="ABIDE_ASD"
export data_path=/headnode1/abry4213/data/${dataset_ID}
for noise_label in FC1000
do
    cmd="python $pyspi_script_dir/distribute_jobs.py \
    --data_dir ${data_path}/raw_data/numpy_files/${noise_label}/ \
    --compute_file $pyspi_script_dir/pyspi_compute.py \
    --template_pbs_file $pyspi_script_dir/template.pbs \
    --pyspi_config $config_file \
    --sample_yaml ${data_path}/raw_data/numpy_files/${noise_label}/sample.yaml \
    --pbs_notify a \
    --email $email --walltime_hrs 6 --cpu 1 --mem 20 --table_only"
    echo $cmd
    $cmd
done