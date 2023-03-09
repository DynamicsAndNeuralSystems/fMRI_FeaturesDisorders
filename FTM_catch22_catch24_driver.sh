##########################################################################################
export github_dir=/headnode1/abry4213/github/
export email="abry4213@uni.sydney.edu.au"
export conda_env="pyspi"
export cluster_queue="yossarian"
export python_to_use=/headnode1/abry4213/.conda/envs/${conda_env}/bin/python3
export pairwise_feature_set="pyspi14"
export scaling_type="robustsigmoid"
export num_folds=10
export num_repeats=10
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

for comparison_group in Schizophrenia; do #Bipolar ADHD; do
    # Iterate over catch2, catch22, and catch24 for the rest of the analyses
    for univariate_feature_set in catch2 catch22 catch24; do
        univariate_feature_file=$data_path/processed_data/${dataset_ID}_${noise_label}_${univariate_feature_set}_filtered.feather

        # # Run univariate SVM
        # # Main SVM
        # num_jobs=10
        # job_memory=20
        # num_hours=3
        # qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_file=$univariate_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
        # -N ${dataset_ID}_${comparison_group}_${univariate_feature_set} \
        # -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_classification_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${scaling_type}_scaler_out.txt \
        # -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
        # $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_univariate_classification.pbs

        # Nulls
        job_memory=60
        num_hours=48
        num_jobs=10
        qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_file=$univariate_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
        -N ${dataset_ID}_${comparison_group}_${univariate_feature_set} \
        -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_nulls_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${scaling_type}_scaler_out.txt \
        -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
        $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_univariate_nulls.pbs

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
export comparison_group="ASD"

# # Iterate over catch2, catch22, and catch24 for the rest of the analyses
for univariate_feature_set in catch2 catch22 catch24; do
    univariate_feature_file=$data_path/processed_data/${dataset_ID}_${noise_label}_${univariate_feature_set}_filtered.feather

    # Run univariate SVM

    # # Main SVM
    # num_jobs=10
    # job_memory=80
    # num_hours=3
    # qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_file=$univariate_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
    # -N ${dataset_ID}_${comparison_group}_${univariate_feature_set} \
    # -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_classification_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${scaling_type}_scaler_out.txt \
    # -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
    # $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_univariate_classification.pbs

    # Nulls
    num_jobs=10
    job_memory=180
    num_hours=48
    qsub -v dataset_ID=$dataset_ID,data_path=$data_path,comparison_group=$comparison_group,univariate_feature_set=$univariate_feature_set,pairwise_feature_set=$pairwise_feature_set,univariate_feature_file=$univariate_feature_file,sample_metadata_file=$sample_metadata_file,noise_proc=$noise_proc,scaling_type=$scaling_type,num_null_iters=$num_null_iters,num_folds=$num_folds,num_jobs=$num_jobs,num_repeats=$num_repeats \
    -N ${dataset_ID}_${comparison_group}_${univariate_feature_set} \
    -o $github_dir/fMRI_FeaturesDisorders/cluster_output/run_univariate_nulls_${dataset_ID}_${comparison_group}_${univariate_feature_set}_${scaling_type}_scaler_out.txt \
    -m a -M $email -l select=1:ncpus=$num_jobs:mem=${job_memory}GB,walltime=${num_hours}:00:00 -q $cluster_queue \
    $github_dir/fMRI_FeaturesDisorders/classification_analysis/call_univariate_nulls.pbs

done