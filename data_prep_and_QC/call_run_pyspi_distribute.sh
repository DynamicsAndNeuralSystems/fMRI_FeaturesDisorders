#!/usr/bin/env bash
export github_dir="$1"
export config_file="$2"
export email="$3"
export dataset_ID="$4"
export data_path="$5"
export noise_procs="$6"
export walltime_hrs="$7"

noise_procs_list=$(echo $noise_procs | sed "s/;/ /g")
export pyspi_script_dir=${github_dir}/pyspi-distribute/

for noise_label in $noise_procs_list
do
    python $pyspi_script_dir/distribute_jobs.py \
    --data_dir ${data_path}/raw_data/numpy_files/${noise_label}/ \
    --compute_file $pyspi_script_dir/pyspi_compute.py \
    --template_pbs_file $pyspi_script_dir/template.pbs \
    --pyspi_config $config_file \
    --sample_yaml ${data_path}/raw_data/numpy_files/${noise_label}/sample.yaml \
    --pbs_notify a \
    --email $email --walltime_hrs $walltime_hrs --cpu 1 --mem 20 --table_only
done