export project_dir="$1"
export github_dir="$2"
export config_file="$3"
export pyspi_script_dir="$4"
export email="$5"

# export project_dir=/headnode1/abry4213
# export github_dir=${project_dir}/github/
# export config_file=${github_dir}/fMRI_FeaturesDisorders/data_prep_and_QC/pyspi19_config.yaml
# export pyspi_script_dir=${project_dir}/github/pyspi-distribute/

# UCLA Schizophrenia
dataset_ID="UCLA_Schizophrenia"
for noise_label in AROMA_2P AROMA_2P_GMR AROMA_2P_DiCER
do
    python $pyspi_script_dir/distribute_jobs.py --data_dir ${project_dir}/data/${dataset_ID}/pydata/${noise_label}/ \
    --compute_file $pyspi_script_dir/pyspi_compute.py \
    --template_pbs_file $pyspi_script_dir/template.pbs \
    --pyspi_config $config_file \
    --sample_yaml ${project_dir}/data/${dataset_ID}/pydata/${noise_label}/sample.yaml \
    --pbs_notify a \
    --email $email --walltime_hrs 2 --cpu 2 --mem 8 --table_only
done


# ABIDE ASD
dataset_ID="ABIDE_ASD"
for noise_label in FC1000
do
    python $pyspi_script_dir/distribute_jobs.py --data_dir ${project_dir}/data/${dataset_ID}/pydata/${noise_label}/ \
    --compute_file $pyspi_script_dir/pyspi_compute.py \
    --template_pbs_file $pyspi_script_dir/template.pbs \
    --pyspi_config $config_file \
    --sample_yaml ${project_dir}/data/${dataset_ID}/pydata/${noise_label}/sample.yaml \
    --pbs_notify a \
    --email $email --walltime_hrs 1 --cpu 2 --mem 8 --table_only
done