##########################################################################################
export github_dir=/headnode1/abry4213/github/
export univariate_feature_set="catch22"
export pairwise_feature_set="pyspi14_mod_di_gaussian"
export email="abry4213@uni.sydney.edu.au"
export python_to_use=/headnode1/abry4213/.conda/envs/pyspi/bin/python3
export pkl_file="calc_di_gaussian.pkl"
export sample_yaml="sample.yaml"
export pyspi_ncpus=2
export pyspi_mem=40

cd $github_dir/fMRI_FeaturesDisorders/data_prep_and_QC/pyspi_QC_analysis/

# UCLA Schizophrenia
export dataset_ID="UCLA_Schizophrenia"
export data_path=/headnode1/abry4213/data/${dataset_ID}/
export sample_metadata_file=${dataset_ID}_sample_metadata.Rds
export brain_region_lookup="Brain_Region_info.csv"
export noise_procs="AROMA+2P+GMR"
export main_noise_proc="AROMA+2P+GMR"
export label_vars="Diagnosis"
export pyspi_walltime_hrs=8

# # ABIDE ASD
# export dataset_ID="ABIDE_ASD"
# export data_path=/headnode1/abry4213/data/${dataset_ID}/
# export sample_metadata_file=${dataset_ID}_sample_metadata.Rds
# export brain_region_lookup="Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv"
# export noise_procs="FC1000"
# export main_noise_proc="FC1000"
# export label_vars="Diagnosis"
# export pyspi_walltime_hrs=6

# # Run pyspi-distribute
# bash call_run_pyspi_distribute.sh \
# $github_dir \
# ${github_dir}/fMRI_FeaturesDisorders/data_prep_and_QC/pyspi14_mod_config.yaml \
# $email \
# $dataset_ID \
# $data_path \
# $noise_procs \
# $pyspi_walltime_hrs \
# $pyspi_mem \
# $pyspi_ncpus \
# $pkl_file \
# $sample_yaml

# Integrate results from pyspi-distribute
qsub -v github_dir=$github_dir,data_path=$data_path,pkl_file=$pkl_file,python_to_use=$python_to_use,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,sample_metadata_file=$sample_metadata_file,brain_region_lookup=$brain_region_lookup,noise_procs=$noise_procs,main_noise_proc=$main_noise_proc,dataset_ID=$dataset_ID \
-N clean_pairwise_di_gaussian_data_${dataset_ID} \
-o ${github_dir}/fMRI_FeaturesDisorders/cluster_output/clean_pairwise_di_gaussian_data_${dataset_ID}_out.txt \
-m a -M $email \
call_clean_pairwise_di_gaussian_data.pbs
