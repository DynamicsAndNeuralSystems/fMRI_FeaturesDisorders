##########################################################################################
export github_dir=/headnode1/abry4213/github/
export pairwise_feature_set="pyspi14"
export email="abry4213@uni.sydney.edu.au"
export conda_env="pyspi"
export pyspi_ncpus=1
export pyspi_mem=40
export num_jobs=16
export cluster_queue="defaultQ"
export pyspi_config=${github_dir}/fMRI_FeaturesDisorders/prep_data_and_QC/pyspi14_config.yaml
export sample_yaml="sample.yaml"
export python_to_use=/headnode1/abry4213/.conda/envs/${conda_env}/bin/python3
export SPI_directionality_file=${github_dir}/fMRI_FeaturesDisorders/classification_analysis/SPI_Direction_Info.csv
export pkl_file="calc_pyspi14.pkl"
export label_vars="Diagnosis"
export num_null_iters=1000

#module load Anaconda3-5.1.0
/usr/physics/Modules/3.2.8/bin/modulecmd bash load Anaconda3-5.1.0 --silent
source /usr/physics/python/anaconda3/etc/profile.d/conda.sh 

# Activate the given conda environment
conda activate $conda_env

##########################################################################################
# UCLA CNP

# Define data paths
export dataset_ID="UCLA_CNP"
export data_path=/headnode1/abry4213/data/UCLA_CNP/
export sample_metadata_file="UCLA_CNP_sample_metadata.feather"
export noise_proc="AROMA+2P+GMR"
export noise_label="AROMA_2P_GMR"
export brain_region_lookup="UCLA_CNP_Brain_Region_Lookup.feather"
export pyspi_walltime_hrs=2
# Run univariate and pairwise linear SVM
pairwise_feature_file=$data_path/processed_data/${dataset_ID}_${noise_label}_${pairwise_feature_set}_filtered.feather

for scaling_type in standard robustsigmoid; do 
    for comparison_group in Schizophrenia Bipolar ADHD; do
        for univariate_feature_set in catch2 catch22 catch24; do
            univariate_feature_file=$data_path/processed_data/${dataset_ID}_${noise_label}_${univariate_feature_set}_filtered.feather

            # Run univariate SVM
            job_memory=40
            qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_file=$univariate_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_jobs=$num_jobs \
            -N ${dataset_ID}_${comparison_group}_${univariate_feature_set} \
            -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_classification_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${scaling_type}_scaler_out.txt \
            -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB -q $cluster_queue \
            $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_univariate_classification.pbs

            # # Run pairwise SVM
            # job_memory=180
            # qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,pairwise_feature_file=$pairwise_feature_file,SPI_directionality_file=$SPI_directionality_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_jobs=$num_jobs \
            # -N ${dataset_ID}_${comparison_group}_${pairwise_feature_set} \
            # -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_pairwise_classification_${dataset_ID}_${comparison_group}_${pairwise_feature_set}_${scaling_type}_scaler_out.txt \
            # -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB -q $cluster_queue \
            # $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_pairwise_classification.pbs
            
            # # Run combined univariate+pairwise SVM
            # job_memory=180
            # qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,univariate_feature_file=$univariate_feature_file,pairwise_feature_set=$pairwise_feature_set,pairwise_feature_file=$pairwise_feature_file,SPI_directionality_file=$SPI_directionality_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iter,num_jobs=$num_jobs \
            # -N ${dataset_ID}_${comparison_group}_${univariate_feature_set}_${pairwise_feature_set} \
            # -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_combo_classification_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${pairwise_feature_set}_${scaling_type}_scaler_out.txt \
            # -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GBB -q $cluster_queue  \
            # $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_combined_univariate_pairwise_classification.pbs
        done
    done
done

##########################################################################################
# ABIDE ASD

# Define data paths
export dataset_ID="ABIDE_ASD"
export data_path=/headnode1/abry4213/data/ABIDE_ASD/
export sample_metadata_file=ABIDE_ASD_sample_metadata.feather
export noise_proc="FC1000"
export noise_label="FC1000"
export brain_region_lookup="ABIDE_ASD_Brain_Region_Lookup.feather"
export pyspi_walltime_hrs=3

# Run univariate and pairwise linear SVM
pairwise_feature_file=$data_path/processed_data/${dataset_ID}_${noise_label}_${pairwise_feature_set}_filtered.feather

# for scaling_type in standard robust; do 
#     for comparison_group in ASD; do
#         for univariate_feature_set in catch2 catch22 catch24; do
#             univariate_feature_file=$data_path/processed_data/ABIDE_ASD_FC1000_${univariate_feature_set}_filtered.feather

#             # Run univariate SVM
#             qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_file=$univariate_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_jobs=$num_jobs \
#             -N ${dataset_ID}_${comparison_group}_${univariate_feature_set} \
#             -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_classification_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${scaling_type}_scaler_out.txt \
#             -m a -M $email \
#             $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_univariate_classification.pbs

#             # Run pairwise SVM
#             qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,pairwise_feature_file=$pairwise_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_jobs=$num_jobs \
#             -N ${dataset_ID}_${comparison_group}_${pairwise_feature_set} \
#             -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_pairwise_classification_${dataset_ID}_${comparison_group}_${pairwise_feature_set}_${scaling_type}_scaler_out.txt \
#             -m a -M $email \
#             $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_pairwise_classification.pbs

#             # Run combined univariate+pairwise SVM
#             qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,univariate_feature_file=$univariate_feature_file,pairwise_feature_set=$pairwise_feature_set,pairwise_feature_file=$pairwise_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_jobs=$num_jobs \
#             -N ${dataset_ID}_${comparison_group}_${univariate_feature_set}_${pairwise_feature_set} \
#             -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_combo_classification_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${pairwise_feature_set}_${scaling_type}_scaler_out.txt \
#             -m a -M $email \
#             $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_combined_univariate_pairwise_classification.pbs
#         done
#     done
# done